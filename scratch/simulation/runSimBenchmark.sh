#!/usr/bin/env bash
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

#
# Simulate a diploid sample with known SVs, call it with an installed manta build, and compare to truth.
#
# Usage:
#   runSimBenchmark.sh MANTA_INSTALL_DIR REFERENCE_FASTA REGION WORK_DIR [LABEL] [-- extra configManta.py arguments]
#
# Environment:
#   SIM_SEED       simulation seed (default 7)
#   SIM_COVERAGE   simulated depth (default 30)
#   SIM_TE_GTF     optional repeatmasker style GTF used to pick mobile element insertion sequences
#   SIM_JOBS       cores given to the manta workflow (default 2)
#   SIM_REUSE      if set to 1, reuse the simulated alignments already in WORK_DIR
#

set -o nounset
set -o pipefail
set -o errexit

if [ $# -lt 4 ]; then
    sed -n '/^# Simulate/,/^$/p' "$0" 1>&2
    exit 2
fi

mantaDir=$1
refFasta=$2
region=$3
workDir=$4
shift 4
label=sim
if [ $# -gt 0 ] && [ "$1" != "--" ]; then label=$1; shift; fi
if [ $# -gt 0 ] && [ "$1" == "--" ]; then shift; fi

scriptDir=$(cd "$(dirname "$0")" && pwd)
samtools=${mantaDir}/libexec/samtools
seed=${SIM_SEED:-7}
coverage=${SIM_COVERAGE:-30}
jobs=${SIM_JOBS:-2}

mkdir -p "${workDir}"
prefix=${workDir}/sim

if [ "${SIM_REUSE:-0}" != "1" ] || [ ! -f "${prefix}.bam.bai" ]; then
    teArg=""
    if [ -n "${SIM_TE_GTF:-}" ]; then teArg="--te-gtf ${SIM_TE_GTF}"; fi
    python3 "${scriptDir}/simulateSVReads.py" --ref "${refFasta}" --region "${region}" --out-prefix "${prefix}" \
        --seed "${seed}" --coverage "${coverage}" ${teArg}
    "${samtools}" faidx "${prefix}.ref.fa"
    "${samtools}" sort -@ "${jobs}" -o "${prefix}.bam" "${prefix}.sam" 2>/dev/null
    "${samtools}" index "${prefix}.bam"
    rm -f "${prefix}.sam"
fi

runDir=${workDir}/manta.${label}
rm -rf "${runDir}"
python3 "${mantaDir}/bin/configManta.py" --bam "${prefix}.bam" --referenceFasta "${prefix}.ref.fa" \
    --runDir "${runDir}" "$@" >/dev/null
python3 "${runDir}/runWorkflow.py" -j "${jobs}" >"${runDir}.log" 2>&1 || {
    echo "manta workflow failed, see ${runDir}.log" 1>&2
    tail -n 30 "${runDir}.log" 1>&2
    exit 1
}

# single-record inversions make the comparison to truth simpler
python3 "${mantaDir}/libexec/convertInversion.py" "${samtools}" "${prefix}.ref.fa" \
    "${runDir}/results/variants/diploidSV.vcf.gz" >"${runDir}/diploidSV.inv.vcf"

python3 "${scriptDir}/evaluateCalls.py" --truth "${prefix}.truth.tsv" --vcf "${runDir}/diploidSV.inv.vcf" \
    --label "${label}" --json "${runDir}/eval.json" --list-errors >"${runDir}/eval.txt"
sed -n '1,/^false negatives:/p' "${runDir}/eval.txt" | sed '$d'
