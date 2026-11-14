#!/usr/bin/env python3
#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2019 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

"""
Simulate a diploid sample carrying a known set of structural variants and write
paired-end reads which are already aligned to the reference.

No read mapper is involved. Reads are sampled from two simulated haplotypes, and
each read's alignment is *projected* from the known haplotype-to-reference
coordinate map, imitating the conventions of bwa-mem: breakpoint-spanning reads
become a soft-clipped primary alignment plus hard-clipped supplementary
alignments linked by SA tags, small indels are placed inline in the CIGAR, and
reads falling entirely inside novel inserted sequence are unmapped.

This makes the output an optimistic test case: there is no mismapping of the kind
a real mapper produces in repetitive sequence, other than a mapping quality
penalty in sequence which is not unique within the simulated region. It is meant
for regression testing and for checking that a method change does what it is
designed to do, not as a substitute for benchmarking on real data.

Three sources of noise are simulated so that a loss of precision can be seen:
chimeric read pairs, reads ending in sequence unrelated to the genome, and
"artifact" variants which are present in a small fraction of reads only (as if
from a third, low-coverage haplotype). Artifact variants are not part of the truth
set, so a call at one of them counts as a false positive.

Outputs (using the --out-prefix value as a file name prefix):
  PREFIX.ref.fa      reference sequence of the simulated region
  PREFIX.sam         unsorted alignments
  PREFIX.truth.tsv   one line per simulated variant allele
  PREFIX.artifacts.tsv  one line per low-fraction artifact variant
"""

import argparse
import bisect
import os
import sys

import numpy as np


BASES = b"ACGT"
COMP_TABLE = bytes.maketrans(b"ACGTN", b"TGCAN")


def revcomp(seq):
    return seq.translate(COMP_TABLE)[::-1]


def log(msg):
    sys.stderr.write("[simulateSVReads] %s\n" % (msg))


#
# reference access
#

class Fasta(object):
    """
    Minimal indexed fasta reader (requires a samtools style .fai file)
    """

    def __init__(self, path):
        self.path = path
        self.index = {}
        faiPath = path + ".fai"
        if not os.path.isfile(faiPath):
            raise Exception("Can't find fasta index: '%s'" % (faiPath))
        for line in open(faiPath):
            w = line.rstrip("\n").split("\t")
            self.index[w[0]] = tuple(int(x) for x in w[1:5])
        self.fp = open(path, "rb")

    def fetch(self, chrom, start, end):
        (length, offset, lineBases, lineWidth) = self.index[chrom]
        if start < 0 or end > length or start >= end:
            raise Exception("Invalid fetch region %s:%i-%i" % (chrom, start, end))
        self.fp.seek(offset + (start // lineBases) * lineWidth + (start % lineBases))
        readSize = (end - start) + ((end - start) // lineBases + 2) * (lineWidth - lineBases)
        return self.fp.read(readSize).replace(b"\n", b"").replace(b"\r", b"")[: end - start].upper()


def parseRegion(regionStr):
    (chrom, span) = regionStr.rsplit(":", 1)
    (start, end) = span.replace(",", "").split("-")
    return (chrom, int(start) - 1, int(end))


#
# variant model
#

class Edit(object):
    """
    A single change applied to one haplotype

    kind is one of:
    DEL     reference [pos,end) is removed
    INS     seq is inserted before reference position pos (pos == end)
    DELINS  reference [pos,end) is replaced by seq
    DUP     reference [pos,end) is duplicated in tandem
    INV     reference [pos,end) is inverted
    MEI     seq is inserted after reference position (end-1) and the target site [pos,end) is duplicated
            on both sides of it
    """

    def __init__(self, kind, pos, end, seq=b""):
        self.kind = kind
        self.pos = pos
        self.end = end
        self.seq = seq

    def span(self):
        return (self.pos, self.end)


class TruthRecord(object):
    def __init__(self, svType, pos, end, svLen, hapMask, group, note):
        self.svType = svType
        self.pos = pos
        self.end = end
        self.svLen = svLen
        self.hapMask = hapMask
        self.group = group
        self.note = note


def editToTruth(edit, hapMask, group, note):
    k = edit.kind
    if k == "DEL":
        return TruthRecord("DEL", edit.pos, edit.end, edit.end - edit.pos, hapMask, group, note)
    if k == "INS":
        return TruthRecord("INS", edit.pos, edit.pos, len(edit.seq), hapMask, group, note)
    if k == "DELINS":
        return TruthRecord("DELINS", edit.pos, edit.end, len(edit.seq) - (edit.end - edit.pos), hapMask, group, note)
    if k == "DUP":
        return TruthRecord("DUP", edit.pos, edit.end, edit.end - edit.pos, hapMask, group, note)
    if k == "INV":
        return TruthRecord("INV", edit.pos, edit.end, edit.end - edit.pos, hapMask, group, note)
    if k == "MEI":
        return TruthRecord("INS", edit.end, edit.end, len(edit.seq) + (edit.end - edit.pos), hapMask, group, note)
    raise Exception("Unknown edit kind: " + k)


def logUniformInt(rng, low, high):
    return int(round(np.exp(rng.uniform(np.log(low), np.log(high)))))


class InsertSequenceSource(object):
    """
    Provides insertion sequences which are real genomic sequence, but not present in the simulated region
    """

    def __init__(self, fasta, chrom, simStart, simEnd, rng, teGtf):
        self.fasta = fasta
        self.chrom = chrom
        self.simStart = simStart
        self.simEnd = simEnd
        self.rng = rng
        self.chromLen = fasta.index[chrom][0]
        self.alus = []
        if teGtf is not None:
            for line in open(teGtf):
                w = line.split("\t")
                if len(w) < 9 or w[0] != chrom: continue
                if 'gene_id "AluY' not in w[8]: continue
                (s, e) = (int(w[3]) - 1, int(w[4]))
                if (e - s) < 280 or (e - s) > 320: continue
                if not (e < simStart - 100000 or s > simEnd + 100000): continue
                self.alus.append((s, e, w[6]))

    def novel(self, size):
        for _ in range(1000):
            s = int(self.rng.integers(1000000, self.chromLen - size - 1000000))
            if not (s + size < self.simStart - 100000 or s > self.simEnd + 100000): continue
            seq = self.fasta.fetch(self.chrom, s, s + size)
            if b"N" in seq: continue
            if self.rng.random() < 0.5: seq = revcomp(seq)
            return seq
        raise Exception("Can't find a source for novel insertion sequence")

    def mobileElement(self):
        """
        return (sequence, isReversed) for an Alu-like insertion with a fresh poly-A tail
        """
        if not self.alus:
            return (None, False)
        (s, e, strand) = self.alus[int(self.rng.integers(0, len(self.alus)))]
        seq = self.fasta.fetch(self.chrom, s, e)
        if b"N" in seq: return (None, False)
        if strand == "-": seq = revcomp(seq)
        seq = seq.rstrip(b"A") + b"A" * int(self.rng.integers(15, 31))
        isReversed = (self.rng.random() < 0.5)
        if isReversed: seq = revcomp(seq)
        return (seq, isReversed)


def makeVariants(args, rng, regionLen, insertSource):
    """
    return (edits per haplotype, truth records)
    """
    hapEdits = ([], [])
    truth = []

    siteSpacing = args.site_spacing
    siteCount = (regionLen - 2 * args.edge_margin) // siteSpacing

    # (label, weight)
    classes = [
        ("del_small", 22), ("del_medium", 12), ("del_large", 4),
        ("ins_small", 14), ("ins_medium", 6), ("ins_large", 5), ("ins_mei", 8),
        ("dup_small", 6), ("dup_large", 5),
        ("inv", 6),
        ("tiny_indel", 4),
        ("compound_del_del", 5), ("compound_del_nested", 3), ("compound_del_ins", 3), ("compound_ins_ins", 2),
    ]
    if args.class_weights is not None:
        # replace the default class mixture, e.g. "dup_small=1,dup_large=2"
        knownLabels = set(c[0] for c in classes)
        classes = []
        for word in args.class_weights.split(","):
            (label, weight) = word.split("=")
            if label not in knownLabels:
                raise Exception("Unknown variant class '%s', expected one of: %s" %
                                (label, ",".join(sorted(knownLabels))))
            classes.append((label, float(weight)))
    labels = [c[0] for c in classes]
    weights = np.array([c[1] for c in classes], dtype=float)
    weights /= weights.sum()

    def zygosity():
        r = rng.random()
        if r < 0.35: return (0, 1)
        if r < 0.675: return (0,)
        return (1,)

    def addSimple(edit, group, note, haps=None):
        if haps is None: haps = zygosity()
        for h in haps: hapEdits[h].append(edit)
        hapMask = sum(1 << h for h in haps)
        truth.append(editToTruth(edit, hapMask, group, note))

    for siteIndex in range(siteCount):
        center = args.edge_margin + siteIndex * siteSpacing + int(rng.integers(0, siteSpacing // 4))
        label = labels[int(rng.choice(len(labels), p=weights))]
        group = "site%04i" % (siteIndex)

        if label == "del_small":
            size = logUniformInt(rng, 50, 300)
            addSimple(Edit("DEL", center, center + size), group, label)
        elif label == "del_medium":
            size = logUniformInt(rng, 300, 2000)
            addSimple(Edit("DEL", center, center + size), group, label)
        elif label == "del_large":
            size = logUniformInt(rng, 2000, 8000)
            addSimple(Edit("DEL", center, center + size), group, label)
        elif label == "ins_small":
            size = logUniformInt(rng, 50, 150)
            addSimple(Edit("INS", center, center, insertSource.novel(size)), group, label)
        elif label == "ins_medium":
            size = logUniformInt(rng, 150, 400)
            addSimple(Edit("INS", center, center, insertSource.novel(size)), group, label)
        elif label == "ins_large":
            size = logUniformInt(rng, 500, 3000)
            addSimple(Edit("INS", center, center, insertSource.novel(size)), group, label)
        elif label == "ins_mei":
            (seq, isReversed) = insertSource.mobileElement()
            if seq is None:
                addSimple(Edit("INS", center, center, insertSource.novel(300)), group, "ins_medium")
            else:
                tsd = int(rng.integers(10, 17))
                addSimple(Edit("MEI", center, center + tsd, seq), group, label + ("_rev" if isReversed else "_fwd"))
        elif label == "dup_small":
            size = logUniformInt(rng, 100, 600)
            addSimple(Edit("DUP", center, center + size), group, label)
        elif label == "dup_large":
            size = logUniformInt(rng, 600, 6000)
            addSimple(Edit("DUP", center, center + size), group, label)
        elif label == "inv":
            size = logUniformInt(rng, 500, 6000)
            addSimple(Edit("INV", center, center + size), group, label)
        elif label == "tiny_indel":
            size = int(rng.integers(10, 50))
            if rng.random() < 0.5:
                addSimple(Edit("DEL", center, center + size), group, label)
            else:
                addSimple(Edit("INS", center, center, insertSource.novel(size)), group, label)
        elif label == "compound_del_del":
            # two partially overlapping deletions, one per haplotype
            sizeA = logUniformInt(rng, 80, 600)
            sizeB = logUniformInt(rng, 80, 600)
            shift = int(rng.integers(max(20, sizeA // 4), max(21, (3 * sizeA) // 4)))
            addSimple(Edit("DEL", center, center + sizeA), group, label, (0,))
            addSimple(Edit("DEL", center + shift, center + shift + sizeB), group, label, (1,))
        elif label == "compound_del_nested":
            # two deletions sharing the same start, as seen at tandem repeats
            sizeA = logUniformInt(rng, 60, 300)
            sizeB = sizeA + logUniformInt(rng, 50, 300)
            addSimple(Edit("DEL", center, center + sizeA), group, label, (0,))
            addSimple(Edit("DEL", center, center + sizeB), group, label, (1,))
        elif label == "compound_del_ins":
            sizeA = logUniformInt(rng, 80, 500)
            sizeB = logUniformInt(rng, 60, 250)
            insPos = center + int(rng.integers(10, sizeA - 10))
            addSimple(Edit("DEL", center, center + sizeA), group, label, (0,))
            addSimple(Edit("INS", insPos, insPos, insertSource.novel(sizeB)), group, label, (1,))
        elif label == "compound_ins_ins":
            sizeA = logUniformInt(rng, 60, 200)
            sizeB = logUniformInt(rng, 60, 200)
            addSimple(Edit("INS", center, center, insertSource.novel(sizeA)), group, label, (0,))
            addSimple(Edit("INS", center, center, insertSource.novel(sizeB)), group, label, (1,))
        else:
            raise Exception("unhandled variant class: " + label)

    return (hapEdits, truth)


#
# haplotype construction
#

class Haplotype(object):
    """
    A haplotype is stored as an ordered list of segments, each one of:
    ('R', refStart, refEnd)  forward copy of reference [refStart,refEnd)
    ('V', refStart, refEnd)  reverse-complement copy of reference [refStart,refEnd)
    ('I', sequence)          sequence which is not in the reference
    """

    def __init__(self, refSeq, edits, snpPositions, snpBases):
        hapRef = bytearray(refSeq)
        for (pos, base) in zip(snpPositions, snpBases):
            hapRef[pos] = base
        hapRef = bytes(hapRef)

        segs = []
        cur = 0
        for edit in sorted(edits, key=lambda x: (x.pos, x.end)):
            if edit.pos < cur:
                raise Exception("Overlapping edits on one haplotype at reference position %i" % (edit.pos))
            k = edit.kind
            if k == "DEL":
                segs.append(("R", cur, edit.pos)); cur = edit.end
            elif k == "INS":
                segs.append(("R", cur, edit.pos)); segs.append(("I", edit.seq)); cur = edit.pos
            elif k == "DELINS":
                segs.append(("R", cur, edit.pos)); segs.append(("I", edit.seq)); cur = edit.end
            elif k == "DUP":
                segs.append(("R", cur, edit.end)); segs.append(("R", edit.pos, edit.end)); cur = edit.end
            elif k == "INV":
                segs.append(("R", cur, edit.pos)); segs.append(("V", edit.pos, edit.end)); cur = edit.end
            elif k == "MEI":
                segs.append(("R", cur, edit.end)); segs.append(("I", edit.seq)); cur = edit.pos
            else:
                raise Exception("Unknown edit kind: " + k)
        segs.append(("R", cur, len(refSeq)))
        self.segs = [s for s in segs if (s[0] == "I" and len(s[1]) > 0) or (s[0] != "I" and s[2] > s[1])]

        parts = []
        self.segHapStart = []
        offset = 0
        for seg in self.segs:
            self.segHapStart.append(offset)
            if seg[0] == "R":
                part = hapRef[seg[1]:seg[2]]
            elif seg[0] == "V":
                part = revcomp(hapRef[seg[1]:seg[2]])
            else:
                part = seg[1]
            parts.append(part)
            offset += len(part)
        self.seq = b"".join(parts)
        self.length = offset
        self.segHapStartArray = np.array(self.segHapStart, dtype=np.int64)
        self.segIsForwardRef = np.array([s[0] == "R" for s in self.segs], dtype=bool)
        self.segRefStart = np.array([(s[1] if s[0] != "I" else -1) for s in self.segs], dtype=np.int64)

    def pieces(self, hapStart, hapEnd):
        """
        Describe the haplotype interval [hapStart,hapEnd) as a list of pieces:
        (readOffset, length, refStart, refEnd, strand) with strand in ('+','-') for reference derived pieces,
        and (readOffset, length, None, None, None) for inserted sequence
        """
        ret = []
        segIndex = bisect.bisect_right(self.segHapStart, hapStart) - 1
        while segIndex < len(self.segs):
            segStart = self.segHapStart[segIndex]
            if segStart >= hapEnd: break
            seg = self.segs[segIndex]
            segLen = (len(seg[1]) if seg[0] == "I" else (seg[2] - seg[1]))
            a = max(hapStart, segStart)
            b = min(hapEnd, segStart + segLen)
            if b > a:
                if seg[0] == "R":
                    ret.append((a - hapStart, b - a, seg[1] + (a - segStart), seg[1] + (b - segStart), "+"))
                elif seg[0] == "V":
                    ret.append((a - hapStart, b - a, seg[2] - (b - segStart), seg[2] - (a - segStart), "-"))
                else:
                    ret.append((a - hapStart, b - a, None, None, None))
            segIndex += 1
        return ret


#
# alignment projection
#

class Aln(object):
    """
    One alignment of a read, in reference orientation
    """
    __slots__ = ("pos", "end", "cigar", "isReverse", "isChainReverse", "queryStart", "queryEnd", "aligned")

    def __init__(self, pos, end, cigar, isReverse, isChainReverse, queryStart, queryEnd, aligned):
        self.pos = pos
        self.end = end
        self.cigar = cigar                    # list of (length, op) in reference order, without clipping
        self.isReverse = isReverse            # alignment strand as written to the SAM flag
        self.isChainReverse = isChainReverse  # true if the haplotype is inverted relative to the reference here
        self.queryStart = queryStart          # aligned span in haplotype-forward read coordinates
        self.queryEnd = queryEnd
        self.aligned = aligned


def cigarString(aln, readLen, clipOp):
    lead = aln.queryStart
    trail = readLen - aln.queryEnd
    if aln.isChainReverse: (lead, trail) = (trail, lead)
    out = []
    if lead: out.append("%i%s" % (lead, clipOp))
    for (length, op) in aln.cigar: out.append("%i%s" % (length, op))
    if trail: out.append("%i%s" % (trail, clipOp))
    return "".join(out)


def maskJunkPieces(pieces, readLen, junkLeft, junkRight):
    """
    Mark the first junkLeft and last junkRight bases of a read (in haplotype-forward read coordinates) as
    sequence which does not come from the genome
    """
    validStart = junkLeft
    validEnd = readLen - junkRight
    ret = []
    if junkLeft: ret.append((0, junkLeft, None, None, None))
    for (qOff, length, refStart, refEnd, strand) in pieces:
        a = max(qOff, validStart)
        b = min(qOff + length, validEnd)
        if b <= a: continue
        if strand is None:
            ret.append((a, b - a, None, None, None))
        elif strand == "+":
            ret.append((a, b - a, refStart + (a - qOff), refStart + (b - qOff), strand))
        else:
            ret.append((a, b - a, refEnd - (b - qOff), refEnd - (a - qOff), strand))
    if junkRight: ret.append((validEnd, junkRight, None, None, None))
    return ret


def projectRead(hap, hapStart, hapEnd, isReadReverse, opt, junkLeft=0, junkRight=0):
    """
    Project one read, covering haplotype interval [hapStart,hapEnd), onto the reference

    Returns a list of Aln objects, primary alignment first. The list is empty if the read is unmapped.
    """
    readLen = hapEnd - hapStart
    pieces = hap.pieces(hapStart, hapEnd)
    if junkLeft or junkRight:
        pieces = maskJunkPieces(pieces, readLen, junkLeft, junkRight)

    # group pieces into chains which a mapper would report as a single gapped alignment
    chains = []
    current = None
    pendingInsert = None
    for piece in pieces:
        (qOff, length, refStart, refEnd, strand) = piece
        if strand is None:
            if pendingInsert is not None:
                piece = (pendingInsert[0], pendingInsert[1] + length, None, None, None)
            pendingInsert = piece
            continue
        isJoined = False
        if current is not None and current["strand"] == strand:
            last = current["pieces"][-1]
            insLen = (pendingInsert[1] if pendingInsert is not None else 0)
            if strand == "+":
                delLen = refStart - last[3]
            else:
                delLen = last[2] - refEnd
            if (0 <= delLen <= opt.max_inline_indel) and (insLen <= opt.max_inline_indel) and \
                    (qOff == last[0] + last[1] + insLen):
                current["pieces"].append(piece)
                current["gaps"].append((delLen, insLen))
                isJoined = True
        if not isJoined:
            current = {"strand": strand, "pieces": [piece], "gaps": []}
            chains.append(current)
        pendingInsert = None

    alns = []
    for chain in chains:
        cp = chain["pieces"]
        gaps = chain["gaps"]
        # a gapped alignment is only extended through an indel when enough sequence is aligned beyond it
        while len(cp) > 1 and cp[0][1] < opt.min_indel_flank:
            cp = cp[1:]; gaps = gaps[1:]
        while len(cp) > 1 and cp[-1][1] < opt.min_indel_flank:
            cp = cp[:-1]; gaps = gaps[:-1]

        aligned = sum(p[1] for p in cp)
        if aligned < opt.min_aligned: continue

        ops = [(cp[0][1], "M")]
        for (piece, (delLen, insLen)) in zip(cp[1:], gaps):
            if insLen: ops.append((insLen, "I"))
            if delLen: ops.append((delLen, "D"))
            if (not insLen) and (not delLen) and ops[-1][1] == "M":
                ops[-1] = (ops[-1][0] + piece[1], "M")
            else:
                ops.append((piece[1], "M"))

        isChainReverse = (chain["strand"] == "-")
        if isChainReverse:
            ops = ops[::-1]
            pos = cp[-1][2]
            end = cp[0][3]
        else:
            pos = cp[0][2]
            end = cp[-1][3]
        alns.append(Aln(pos, end, ops, (isReadReverse != isChainReverse), isChainReverse,
                        cp[0][0], cp[-1][0] + cp[-1][1], aligned))

    if not alns: return alns

    # mappers do not clip a handful of bases off a read end, they align through them instead
    refLen = opt.refLen
    for aln in alns:
        lead = aln.queryStart
        trail = readLen - aln.queryEnd
        if aln.isChainReverse: (refLead, refTrail) = (trail, lead)
        else:                  (refLead, refTrail) = (lead, trail)
        if 0 < refLead <= opt.max_absorbed_clip and (aln.pos - refLead) >= 0 and aln.cigar[0][1] == "M":
            aln.pos -= refLead
            aln.cigar[0] = (aln.cigar[0][0] + refLead, "M")
            if aln.isChainReverse: aln.queryEnd += refLead
            else:                  aln.queryStart -= refLead
        if 0 < refTrail <= opt.max_absorbed_clip and (aln.end + refTrail) <= refLen and aln.cigar[-1][1] == "M":
            aln.end += refTrail
            aln.cigar[-1] = (aln.cigar[-1][0] + refTrail, "M")
            if aln.isChainReverse: aln.queryStart -= refTrail
            else:                  aln.queryEnd += refTrail

    primaryIndex = max(range(len(alns)), key=lambda i: alns[i].aligned)
    return [alns[primaryIndex]] + [a for (i, a) in enumerate(alns) if i != primaryIndex]


def getNonUniqueKmerPrefixSum(refSeq, k):
    """
    return an array P with P[i] = number of reference positions < i which start a k-mer occurring more than once
    in the region (on either strand)
    """
    codes = np.frombuffer(refSeq.translate(bytes.maketrans(b"ACGTN", b"\x00\x01\x02\x03\x00")), dtype=np.uint8)
    n = len(codes) - k + 1
    fwd = np.zeros(n, dtype=np.uint64)
    rev = np.zeros(n, dtype=np.uint64)
    for j in range(k):
        b = codes[j:j + n].astype(np.uint64)
        fwd = (fwd << np.uint64(2)) | b
        rev = rev | ((np.uint64(3) - b) << np.uint64(2 * j))
    canonical = np.minimum(fwd, rev)
    (_, inverse, counts) = np.unique(canonical, return_inverse=True, return_counts=True)
    isRepeat = (counts[inverse] > 1)
    prefix = np.zeros(len(codes) + 1, dtype=np.int64)
    prefix[1:n + 1] = np.cumsum(isRepeat)
    prefix[n + 1:] = prefix[n]
    return prefix


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ref", required=True, help="indexed reference fasta")
    parser.add_argument("--region", required=True, help="region to simulate, e.g. chr21:32000001-35000000")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--te-gtf", default=None,
                        help="optional repeatmasker style GTF, used to pick real Alu copies as mobile element insertions")
    parser.add_argument("--sample", default="SIM1")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--coverage", type=float, default=30.)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--fragment-mean", type=float, default=450.)
    parser.add_argument("--fragment-sd", type=float, default=60.)
    parser.add_argument("--site-spacing", type=int, default=15000, help="distance between simulated variant sites")
    parser.add_argument("--edge-margin", type=int, default=20000)
    parser.add_argument("--class-weights", default=None,
                        help="replace the default mixture of variant classes, format: CLASS=WEIGHT[,CLASS=WEIGHT...]")
    parser.add_argument("--snp-rate", type=float, default=0.0008, help="heterozygous SNP rate per haplotype")
    parser.add_argument("--max-inline-indel", type=int, default=25,
                        help="largest indel kept inside a single gapped alignment")
    parser.add_argument("--min-indel-flank", type=int, default=12)
    parser.add_argument("--min-aligned", type=int, default=30, help="shortest alignment which is reported")
    parser.add_argument("--max-absorbed-clip", type=int, default=4)
    parser.add_argument("--kmer-size", type=int, default=32)
    parser.add_argument("--chimera-rate", type=float, default=0.004,
                        help="fraction of fragments whose two reads come from unrelated positions")
    parser.add_argument("--junk-tail-rate", type=float, default=0.01,
                        help="fraction of reads which end in sequence unrelated to the genome")
    parser.add_argument("--artifact-fraction", type=float, default=0.1,
                        help="depth of the artifact haplotype as a fraction of --coverage (0 to disable)")
    parser.add_argument("--artifact-spacing", type=int, default=45000, help="distance between artifact variants")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    fasta = Fasta(args.ref)
    (chrom, regionStart, regionEnd) = parseRegion(args.region)
    refSeq = fasta.fetch(chrom, regionStart, regionEnd)
    if b"N" in refSeq:
        raise Exception("The simulated region must not contain 'N' bases")
    refLen = len(refSeq)
    args.refLen = refLen
    simChrom = chrom
    log("simulated region %s:%i-%i (%i bases)" % (chrom, regionStart + 1, regionEnd, refLen))

    # the simulated region becomes a stand-alone reference, coordinates are relative to its start
    with open(args.out_prefix + ".ref.fa", "w") as fp:
        fp.write(">%s\n" % (simChrom))
        text = refSeq.decode()
        for i in range(0, refLen, 60): fp.write(text[i:i + 60] + "\n")

    insertSource = InsertSequenceSource(fasta, chrom, regionStart, regionEnd, rng, args.te_gtf)
    (hapEdits, truth) = makeVariants(args, rng, refLen, insertSource)
    log("simulated %i variant alleles, %i candidate mobile element source copies" % (len(truth), len(insertSource.alus)))

    # Artifact variants sit halfway between real variant sites, and are only seen in the reads of a third
    # low-coverage haplotype
    artifactEdits = []
    artifactTruth = []
    if args.artifact_fraction > 0:
        # real variants occupy at most the first ~80% of each site interval
        pos = args.edge_margin + int(0.85 * args.site_spacing)
        index = 0
        while pos < (refLen - args.edge_margin):
            group = "artifact%04i" % (index)
            r = rng.random()
            if r < 0.5:
                edit = Edit("INS", pos, pos, insertSource.novel(logUniformInt(rng, 60, 600)))
            elif r < 0.85:
                edit = Edit("DEL", pos, pos + logUniformInt(rng, 60, 2000))
            else:
                edit = Edit("DUP", pos, pos + logUniformInt(rng, 100, 2000))
            artifactEdits.append(edit)
            artifactTruth.append(editToTruth(edit, 4, group, "artifact"))
            index += 1
            # keep artifacts on the grid of positions between real variant sites
            pos += max(1, int(round(args.artifact_spacing / float(args.site_spacing)))) * args.site_spacing
    with open(args.out_prefix + ".artifacts.tsv", "w") as fp:
        fp.write("#chrom\tpos0\tend0\ttype\tsvlen\tgenotype\thapMask\tgroup\tclass\n")
        for t in artifactTruth:
            fp.write("%s\t%i\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" %
                     (simChrom, t.pos, t.end, t.svType, t.svLen, "artifact", t.hapMask, t.group, t.note))
    log("simulated %i artifact variants at %.0f%% of sample depth" % (len(artifactTruth), 100 * args.artifact_fraction))

    with open(args.out_prefix + ".truth.tsv", "w") as fp:
        fp.write("#chrom\tpos0\tend0\ttype\tsvlen\tgenotype\thapMask\tgroup\tclass\n")
        for t in sorted(truth, key=lambda x: (x.pos, x.end)):
            gt = {1: "1|0", 2: "0|1", 3: "1|1"}[t.hapMask]
            fp.write("%s\t%i\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" %
                     (simChrom, t.pos, t.end, t.svType, t.svLen, gt, t.hapMask, t.group, t.note))

    # heterozygous SNPs, kept away from the simulated breakpoints so that the truth coordinates stay exact
    blocked = np.zeros(refLen, dtype=bool)
    for edits in (hapEdits[0], hapEdits[1], artifactEdits):
        for e in edits: blocked[max(0, e.pos - 5):e.end + 5] = True
    haps = []
    for hapIndex in range(2):
        snpCount = rng.binomial(refLen, args.snp_rate)
        snpPos = np.unique(rng.integers(0, refLen, snpCount))
        snpPos = snpPos[~blocked[snpPos]]
        snpBases = []
        for p in snpPos:
            refBase = refSeq[p:p + 1]
            alts = [b for b in (b"A", b"C", b"G", b"T") if b != refBase]
            snpBases.append(alts[int(rng.integers(0, 3))][0])
        haps.append(Haplotype(refSeq, hapEdits[hapIndex], snpPos.tolist(), snpBases))
        log("haplotype %i: %i bases, %i segments, %i SNPs" % (hapIndex, haps[-1].length, len(haps[-1].segs), len(snpPos)))

    # (haplotype, depth) pairs to sample reads from
    hapDepths = [(haps[0], args.coverage / 2.), (haps[1], args.coverage / 2.)]
    if artifactEdits:
        hapDepths.append((Haplotype(refSeq, artifactEdits, [], []), args.coverage * args.artifact_fraction))

    nonUniquePrefix = getNonUniqueKmerPrefixSum(refSeq, args.kmer_size)
    log("%.2f%% of reference %i-mers are not unique within the region" %
        (100. * nonUniquePrefix[-1] / refLen, args.kmer_size))

    def blockMapq(aln):
        span = aln.end - aln.pos
        if span < args.kmer_size + 8: return 30
        kmerCount = span - args.kmer_size + 1
        uniqueCount = kmerCount - int(nonUniquePrefix[aln.pos + kmerCount] - nonUniquePrefix[aln.pos])
        if uniqueCount == 0: return 0
        if uniqueCount < 10: return 20
        return 60

    readLen = args.read_length
    qualChoices = np.array([37, 25, 12], dtype=np.uint8)
    qualProbs = np.array([0.86, 0.10, 0.04])
    qualErrorProb = np.power(10., -qualChoices.astype(float) / 10.)
    baseCode = np.zeros(256, dtype=np.uint8)
    for (i, b) in enumerate(BASES): baseCode[b] = i
    codeBase = np.frombuffer(BASES, dtype=np.uint8)

    properLow = args.fragment_mean - 4.5 * args.fragment_sd
    properHigh = args.fragment_mean + 4.5 * args.fragment_sd

    samFp = open(args.out_prefix + ".sam", "w")
    samFp.write("@HD\tVN:1.6\tSO:unsorted\n")
    samFp.write("@SQ\tSN:%s\tLN:%i\n" % (simChrom, refLen))
    samFp.write("@RG\tID:sim\tSM:%s\tPL:ILLUMINA\n" % (args.sample))
    samFp.write("@PG\tID:simulateSVReads\tPN:simulateSVReads\tCL:%s\n" % (" ".join(sys.argv)))

    stats = {"fragments": 0, "split": 0, "unmapped": 0, "complex": 0, "chimeras": 0, "junk": 0}

    for (hapIndex, (hap, hapDepth)) in enumerate(hapDepths):
        fragCount = int(hapDepth * hap.length / (2. * readLen))
        fragLen = np.clip(np.rint(rng.normal(args.fragment_mean, args.fragment_sd, fragCount)).astype(np.int64),
                          readLen, int(args.fragment_mean + 6 * args.fragment_sd))
        fragStart = rng.integers(0, hap.length - fragLen)
        isFragReverse = (rng.random(fragCount) < 0.5)

        hapArray = np.frombuffer(hap.seq, dtype=np.uint8)

        # haplotype-forward coordinates of the upstream and downstream read of each fragment
        readStart = [fragStart, fragStart + fragLen - readLen]

        # the downstream read is sequenced from the reverse strand, except in a chimeric fragment where the
        # second read comes from an unrelated position in either orientation
        isChimera = (rng.random(fragCount) < args.chimera_rate)
        readStart[1] = np.where(isChimera, rng.integers(0, hap.length - readLen, fragCount), readStart[1])
        isMateReverse = [np.zeros(fragCount, dtype=bool), np.where(isChimera, rng.random(fragCount) < 0.5, True)]
        stats["chimeras"] += int(isChimera.sum())

        # number of bases at the 3' end of each read which are replaced by unrelated sequence
        junkLen = []
        for mate in range(2):
            isJunk = (rng.random(fragCount) < args.junk_tail_rate)
            junkLen.append(np.where(isJunk, rng.integers(8, 70, fragCount), 0))
            stats["junk"] += int(isJunk.sum())

        # sequence and base qualities for all reads, both in haplotype-forward orientation
        readSeqs = []
        readQuals = []
        for mate in range(2):
            window = baseCode[hapArray[readStart[mate][:, None] + np.arange(readLen)[None, :]]]
            qualIndex = rng.choice(len(qualChoices), size=window.shape, p=qualProbs)
            isError = (rng.random(window.shape) < qualErrorProb[qualIndex])
            shift = rng.integers(1, 4, size=window.shape, dtype=np.uint8)
            window = np.where(isError, (window + shift) % 4, window).astype(np.uint8)
            quals = qualChoices[qualIndex].copy()

            # A junk tail is at the 3' end of the read as sequenced, which is the high-coordinate end of a
            # forward read and the low-coordinate end of a reverse read in haplotype-forward orientation. Half of
            # the junk tails are given a very low base quality.
            column = np.arange(readLen)[None, :]
            junk = junkLen[mate][:, None]
            isJunkBase = np.where(isMateReverse[mate][:, None], column < junk, column >= (readLen - junk))
            window = np.where(isJunkBase, rng.integers(0, 4, size=window.shape, dtype=np.uint8), window).astype(np.uint8)
            isLowQualTail = (rng.random(fragCount) < 0.5)[:, None]
            quals = np.where(isJunkBase & isLowQualTail, 2, quals).astype(np.uint8)

            readSeqs.append(codeBase[window])
            readQuals.append((quals + 33).astype(np.uint8))

        # fast path: reads contained in a single forward reference segment
        simple = []
        simplePos = []
        for mate in range(2):
            s0 = np.searchsorted(hap.segHapStartArray, readStart[mate], side="right") - 1
            s1 = np.searchsorted(hap.segHapStartArray, readStart[mate] + readLen - 1, side="right") - 1
            isSimple = (s0 == s1) & hap.segIsForwardRef[s0] & (junkLen[mate] == 0)
            simple.append(isSimple)
            simplePos.append(hap.segRefStart[s0] + (readStart[mate] - hap.segHapStartArray[s0]))

        simpleCigar = [(readLen, "M")]

        for fragIndex in range(fragCount):
            qname = "sim:h%i:%i" % (hapIndex, fragIndex)
            isRev = bool(isFragReverse[fragIndex])

            # mate 0 is the haplotype-upstream read. It is sequenced from the forward strand of the haplotype,
            # and the downstream read from the reverse strand. A reversed fragment only swaps which of the two is
            # labeled as read 1.
            mateAlns = []
            for mate in range(2):
                isReadReverse = bool(isMateReverse[mate][fragIndex])
                if simple[mate][fragIndex]:
                    pos = int(simplePos[mate][fragIndex])
                    alns = [Aln(pos, pos + readLen, simpleCigar, isReadReverse, False, 0, readLen, readLen)]
                else:
                    hs = int(readStart[mate][fragIndex])
                    junk = int(junkLen[mate][fragIndex])
                    (junkLeft, junkRight) = ((junk, 0) if isReadReverse else (0, junk))
                    alns = projectRead(hap, hs, hs + readLen, isReadReverse, args, junkLeft, junkRight)
                    stats["complex"] += 1
                mateAlns.append(alns)

            stats["fragments"] += 1
            mapqs = [[blockMapq(a) for a in alns] for alns in mateAlns]

            # a read in non-unique sequence is still placed with some confidence when its mate is anchored
            for mate in range(2):
                other = 1 - mate
                if mateAlns[mate] and mateAlns[other] and mapqs[mate][0] == 0 and mapqs[other][0] == 60:
                    mapqs[mate][0] = 20

            prim = [(alns[0] if alns else None) for alns in mateAlns]
            isProper = False
            tlen = [0, 0]
            if prim[0] is not None and prim[1] is not None:
                left = min(prim[0].pos, prim[1].pos)
                right = max(prim[0].end, prim[1].end)
                span = right - left
                if prim[0].isReverse != prim[1].isReverse:
                    fwd = (prim[0] if not prim[0].isReverse else prim[1])
                    rev = (prim[1] if not prim[0].isReverse else prim[0])
                    if fwd.pos <= rev.pos and properLow <= span <= properHigh: isProper = True
                for mate in range(2):
                    other = 1 - mate
                    isLeft = (prim[mate].pos < prim[other].pos) or (prim[mate].pos == prim[other].pos and mate == 0)
                    tlen[mate] = (span if isLeft else -span)

            for mate in range(2):
                other = 1 - mate
                isFirst = ((mate == 0) != isRev)
                fwdSeq = readSeqs[mate][fragIndex].tobytes()
                fwdQual = readQuals[mate][fragIndex].tobytes()
                alns = mateAlns[mate]

                baseFlag = 0x1 | (0x40 if isFirst else 0x80)
                if isProper: baseFlag |= 0x2
                if prim[other] is None: baseFlag |= 0x8
                elif prim[other].isReverse: baseFlag |= 0x20

                if prim[other] is not None:
                    mateFields = ("=", prim[other].pos + 1)
                    mateCigarTag = "\tMC:Z:" + cigarString(prim[other], readLen, "S")
                else:
                    mateFields = ("*", 0)
                    mateCigarTag = ""

                if not alns:
                    stats["unmapped"] += 1
                    # an unmapped read is written in its sequenced orientation, at the position of its mate
                    isReadReverse = bool(isMateReverse[mate][fragIndex])
                    seq = (revcomp(fwdSeq) if isReadReverse else fwdSeq)
                    qual = (fwdQual[::-1] if isReadReverse else fwdQual)
                    flag = baseFlag | 0x4
                    if prim[other] is not None:
                        samFp.write("%s\t%i\t%s\t%i\t0\t*\t=\t%i\t0\t%s\t%s\tRG:Z:sim%s\n" %
                                    (qname, flag, simChrom, prim[other].pos + 1, prim[other].pos + 1,
                                     seq.decode(), qual.decode(), mateCigarTag))
                    else:
                        samFp.write("%s\t%i\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tRG:Z:sim\n" %
                                    (qname, flag, seq.decode(), qual.decode()))
                    continue

                if len(alns) > 1: stats["split"] += 1

                saEntries = []
                for (alnIndex, aln) in enumerate(alns):
                    saEntries.append("%s,%i,%s,%s,%i,0" % (simChrom, aln.pos + 1, ("-" if aln.isReverse else "+"),
                                                            cigarString(aln, readLen, "S"), mapqs[mate][alnIndex]))

                for (alnIndex, aln) in enumerate(alns):
                    isSupplement = (alnIndex > 0)
                    flag = baseFlag | (0x10 if aln.isReverse else 0) | (0x800 if isSupplement else 0)
                    if aln.isChainReverse:
                        seq = revcomp(fwdSeq)
                        qual = fwdQual[::-1]
                    else:
                        seq = fwdSeq
                        qual = fwdQual
                    if isSupplement:
                        # supplementary alignments are hard-clipped, as bwa-mem writes them by default
                        lead = aln.queryStart
                        trail = readLen - aln.queryEnd
                        if aln.isChainReverse: (lead, trail) = (trail, lead)
                        seq = seq[lead:readLen - trail]
                        qual = qual[lead:readLen - trail]
                    tags = "RG:Z:sim" + mateCigarTag
                    others = [s for (i, s) in enumerate(saEntries) if i != alnIndex]
                    if others: tags += "\tSA:Z:" + ";".join(others) + ";"
                    samFp.write("%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n" %
                                (qname, flag, simChrom, aln.pos + 1, mapqs[mate][alnIndex],
                                 cigarString(aln, readLen, ("H" if isSupplement else "S")),
                                 mateFields[0], mateFields[1], (0 if isSupplement and prim[other] is None else tlen[mate]),
                                 seq.decode(), qual.decode(), tags))

    samFp.close()
    log("wrote %(fragments)i fragments: %(complex)i reads projected through a variant or junk tail, "
        "%(split)i split reads, %(unmapped)i unmapped reads, %(chimeras)i chimeric fragments, "
        "%(junk)i reads with junk tails" % stats)


if __name__ == "__main__":
    main()
