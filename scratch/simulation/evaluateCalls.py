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
Compare a manta diploid VCF to the truth table written by simulateSVReads.py

Matching follows the usual conventions for SV benchmarking: a call matches a truth
variant of the same type when both breakends are within --max-distance bases and the
sizes agree to within --min-size-similarity. Inversions are expected in the single
record form written by convertInversion.py. A tandem duplication may be matched by
an insertion call of the same size, which is how small duplications are reported.
"""

import argparse
import gzip
import json
import sys


def openMaybeGz(path):
    if path.endswith(".gz"): return gzip.open(path, "rt")
    return open(path)


def getInfo(infoStr):
    ret = {}
    for item in infoStr.split(";"):
        if "=" in item:
            (k, v) = item.split("=", 1)
            ret[k] = v
        else:
            ret[item] = True
    return ret


class Call(object):
    def __init__(self, line):
        w = line.rstrip("\n").split("\t")
        self.chrom = w[0]
        self.pos = int(w[1])
        self.vid = w[2]
        self.ref = w[3]
        self.alt = w[4]
        self.qual = w[5]
        self.filter = w[6]
        self.info = getInfo(w[7])
        self.svType = self.info.get("SVTYPE", "NA")
        self.end = int(self.info.get("END", self.pos))
        self.isImprecise = ("IMPRECISE" in self.info)
        self.svLen = None
        if "SVLEN" in self.info:
            try: self.svLen = abs(int(self.info["SVLEN"].split(",")[0]))
            except ValueError: pass
        if self.svLen is None and self.svType in ("DEL", "DUP", "INV"):
            self.svLen = self.end - self.pos
        self.gt = "."
        if len(w) > 9:
            fmt = w[8].split(":")
            vals = w[9].split(":")
            if "GT" in fmt: self.gt = vals[fmt.index("GT")]
        self.isPass = (self.filter == "PASS")
        self.matched = None

    def altAlleleCount(self):
        return sum(1 for a in self.gt.replace("|", "/").split("/") if a == "1")


class Truth(object):
    def __init__(self, line):
        w = line.rstrip("\n").split("\t")
        self.chrom = w[0]
        self.pos = int(w[1])
        self.end = int(w[2])
        self.svType = w[3]
        self.svLen = abs(int(w[4]))
        self.genotype = w[5]
        self.hapMask = int(w[6])
        self.group = w[7]
        self.svClass = w[8]
        self.matched = None

    def altAlleleCount(self):
        return bin(self.hapMask).count("1")


def isTypeMatch(truth, call):
    if truth.svType == call.svType: return True
    if truth.svType == "DUP" and call.svType == "INS": return True
    return False


def matchScore(truth, call, opt):
    """
    return None if no match, otherwise a distance (smaller is better)
    """
    if truth.chrom != call.chrom: return None
    if not isTypeMatch(truth, call): return None

    isInsertionLike = (call.svType == "INS")
    if isInsertionLike:
        # insertion calls have a single reference position, a duplication can be reported at either end
        candidates = [truth.pos]
        if truth.svType == "DUP": candidates.append(truth.end)
        dist = min(abs(call.pos - c) for c in candidates)
        if dist > opt.max_distance: return None
    else:
        d1 = abs(call.pos - truth.pos)
        d2 = abs(call.end - truth.end)
        if max(d1, d2) > opt.max_distance: return None
        dist = d1 + d2

    if call.svLen is not None and call.svLen > 0:
        sim = min(call.svLen, truth.svLen) / float(max(call.svLen, truth.svLen))
        if sim < opt.min_size_similarity: return None
        return dist + (1. - sim) * 100.
    # insertion of unknown size (incompletely assembled)
    return dist + 50.


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--min-size", type=int, default=50)
    parser.add_argument("--max-distance", type=int, default=200)
    parser.add_argument("--min-size-similarity", type=float, default=0.7)
    parser.add_argument("--inv-collapse-distance", type=int, default=50,
                        help="INV records with both ends within this distance are treated as one inversion")
    parser.add_argument("--include-filtered", action="store_true", help="also use calls which are not PASS")
    parser.add_argument("--artifacts", default=None,
                        help="table of simulated low-fraction artifact variants, used to label false positives")
    parser.add_argument("--json", default=None, help="write summary to this file")
    parser.add_argument("--label", default="")
    parser.add_argument("--list-errors", action="store_true", help="print false negatives and false positives")
    opt = parser.parse_args()

    truths = []
    for line in open(opt.truth):
        if line.startswith("#"): continue
        t = Truth(line)
        if t.svLen < opt.min_size: continue
        truths.append(t)

    calls = []
    skippedTypes = {}
    for line in openMaybeGz(opt.vcf):
        if line.startswith("#"): continue
        c = Call(line)
        if not (c.isPass or opt.include_filtered): continue
        if c.svType not in ("DEL", "INS", "DUP", "INV"):
            skippedTypes[c.svType] = skippedTypes.get(c.svType, 0) + 1
            continue
        if c.svLen is not None and c.svLen < opt.min_size: continue
        calls.append(c)

    # an inversion is called as two junctions, which convertInversion.py writes as two INV records with (nearly)
    # the same coordinates. These are collapsed so that one inversion is counted once.
    collapsed = []
    for c in calls:
        isDuplicate = False
        if c.svType == "INV":
            for o in collapsed:
                if o.svType == "INV" and o.chrom == c.chrom and \
                        abs(o.pos - c.pos) <= opt.inv_collapse_distance and abs(o.end - c.end) <= opt.inv_collapse_distance:
                    isDuplicate = True
                    break
        if not isDuplicate: collapsed.append(c)
    calls = collapsed

    # greedy one-to-one matching, closest pairs first
    pairs = []
    for (ti, t) in enumerate(truths):
        for (ci, c) in enumerate(calls):
            if abs(c.pos - t.pos) > (opt.max_distance + t.svLen + 10): continue
            s = matchScore(t, c, opt)
            if s is not None: pairs.append((s, ti, ci))
    pairs.sort()
    for (s, ti, ci) in pairs:
        if truths[ti].matched is None and calls[ci].matched is None:
            truths[ti].matched = ci
            calls[ci].matched = ti

    def classKey(t):
        return t.svClass.replace("_fwd", "").replace("_rev", "")

    classes = {}
    for t in truths:
        d = classes.setdefault(classKey(t), {"total": 0, "found": 0, "gtMatch": 0})
        d["total"] += 1
        if t.matched is not None:
            d["found"] += 1
            if calls[t.matched].altAlleleCount() == t.altAlleleCount(): d["gtMatch"] += 1

    typeStats = {}
    for t in truths:
        d = typeStats.setdefault(t.svType, {"truth": 0, "tp": 0, "calls": 0, "fp": 0})
        d["truth"] += 1
        if t.matched is not None: d["tp"] += 1
    for c in calls:
        d = typeStats.setdefault(c.svType, {"truth": 0, "tp": 0, "calls": 0, "fp": 0})
        d["calls"] += 1
        if c.matched is None: d["fp"] += 1

    # label false positives which are calls of a simulated low-fraction artifact
    artifacts = []
    if opt.artifacts:
        for line in open(opt.artifacts):
            if line.startswith("#"): continue
            artifacts.append(Truth(line))
    artifactFP = 0
    artifactsCalled = set()
    for c in calls:
        c.isArtifact = False
        if c.matched is not None: continue
        for (ai, a) in enumerate(artifacts):
            if a.chrom == c.chrom and (a.pos - opt.max_distance) <= c.pos <= (a.end + opt.max_distance):
                c.isArtifact = True
                artifactsCalled.add(ai)
                break
        if c.isArtifact: artifactFP += 1

    totalTruth = len(truths)
    totalFound = sum(1 for t in truths if t.matched is not None)
    totalCalls = len(calls)
    totalFP = sum(1 for c in calls if c.matched is None)
    gtMatch = sum(d["gtMatch"] for d in classes.values())
    recall = totalFound / float(max(1, totalTruth))
    precision = (totalCalls - totalFP) / float(max(1, totalCalls))
    f1 = (2 * recall * precision / (recall + precision)) if (recall + precision) > 0 else 0.

    out = sys.stdout
    if opt.label: out.write("== %s ==\n" % (opt.label))
    out.write("%-22s %6s %6s %7s %8s\n" % ("truth class", "total", "found", "recall", "GT match"))
    for k in sorted(classes):
        d = classes[k]
        out.write("%-22s %6i %6i %7.3f %8s\n" %
                  (k, d["total"], d["found"], d["found"] / float(d["total"]),
                   ("%.3f" % (d["gtMatch"] / float(d["found"])) if d["found"] else "NA")))
    out.write("\n%-8s %6s %6s %7s %6s %6s %9s\n" % ("type", "truth", "TP", "recall", "calls", "FP", "precision"))
    for k in sorted(typeStats):
        d = typeStats[k]
        out.write("%-8s %6i %6i %7s %6i %6i %9s\n" %
                  (k, d["truth"], d["tp"], ("%.3f" % (d["tp"] / float(d["truth"])) if d["truth"] else "NA"),
                   d["calls"], d["fp"],
                   ("%.3f" % ((d["calls"] - d["fp"]) / float(d["calls"])) if d["calls"] else "NA")))
    out.write("\nALL: truth %i found %i recall %.3f | calls %i FP %i precision %.3f | F1 %.3f | GT concordance %s\n" %
              (totalTruth, totalFound, recall, totalCalls, totalFP, precision, f1,
               ("%.3f" % (gtMatch / float(totalFound)) if totalFound else "NA")))
    if artifacts:
        out.write("FP breakdown: %i at simulated artifacts (%i of %i artifacts called), %i elsewhere\n" %
                  (artifactFP, len(artifactsCalled), len(artifacts), totalFP - artifactFP))
    if skippedTypes:
        out.write("call types not evaluated: %s\n" % (", ".join("%s=%i" % kv for kv in sorted(skippedTypes.items()))))

    if opt.list_errors:
        out.write("\nfalse negatives:\n")
        for t in truths:
            if t.matched is None:
                out.write("  FN %s:%i-%i %s len=%i gt=%s %s %s\n" %
                          (t.chrom, t.pos, t.end, t.svType, t.svLen, t.genotype, t.group, t.svClass))
        out.write("false positives:\n")
        for c in calls:
            if c.matched is None:
                out.write("  FP %s:%i-%i %s len=%s gt=%s qual=%s filter=%s id=%s%s\n" %
                          (c.chrom, c.pos, c.end, c.svType, c.svLen, c.gt, c.qual, c.filter, c.vid,
                           (" [artifact]" if c.isArtifact else "")))

    if opt.json:
        with open(opt.json, "w") as fp:
            json.dump({"label": opt.label, "classes": classes, "types": typeStats,
                       "truth": totalTruth, "found": totalFound, "calls": totalCalls, "fp": totalFP,
                       "fpAtArtifacts": artifactFP, "artifacts": len(artifacts), "artifactsCalled": len(artifactsCalled),
                       "recall": recall, "precision": precision, "f1": f1,
                       "gtConcordance": (gtMatch / float(totalFound) if totalFound else None)}, fp, indent=1, sort_keys=True)


if __name__ == "__main__":
    main()
