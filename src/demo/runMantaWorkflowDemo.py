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
Execute small manta demonstration/verification runs
"""

def main() :
    import os,sys

    import difflib
    import gzip
    import re
    import subprocess

    #
    # initialize paths:
    #
    scriptDir=os.path.abspath(os.path.dirname(__file__))
    demoDir=os.path.abspath(os.path.join(scriptDir,os.pardir,"share","demo","manta"))
    dataDir=os.path.join(demoDir,"data")
    expectedBaseDir=os.path.join(demoDir,"expectedResults")

    analysisBaseDir="MantaDemoAnalysis"

    configScript=os.path.join(scriptDir,"configManta.py")

    # error logging
    #
    logfp = sys.stderr

    def dataFile(name) :
        return "\"%s\"" % (os.path.join(dataDir,name))

    # Each demo analysis has a label and a list of configuration arguments. The label is the name of the
    # analysis subdirectory and of the expected results subdirectory.
    #
    analyses = [
        # Simulated germline sample with known variants. See README.md.
        ("SimulatedGermline",
         ["--bam=" + dataFile("SimulatedGermline.chr21_region.bam"),
          "--referenceFasta=" + dataFile("SimulatedGermline.chr21_region.fa")]),

        # Real sequencing reads of a tumor cell line and its matched normal. The two alignment files are
        # analyzed as two samples in a joint germline analysis:
        ("HCC1954",
         ["--bam=" + dataFile("HCC1954.NORMAL.30x.compare.COST16011_region.bam"),
          "--bam=" + dataFile("G15512.HCC1954.1.COST16011_region.bam"),
          "--referenceFasta=" + dataFile("Homo_sapiens_assembly19.COST16011_region.fa"),
          "--region=8:107652000-107655000",
          "--region=11:94974000-94989000",
          "--exome"]),
    ]

    if not os.path.exists(configScript) :
        logfp.write("\n")
        logfp.write("ERROR: Manta workflow must be installed prior to running demo.\n")
        logfp.write("\n")
        sys.exit(2)

    if os.path.exists(analysisBaseDir) :
        logfp.write("\n")
        logfp.write("ERROR: Demo analysis directory already exists: '" + analysisBaseDir + "'\n")
        logfp.write("       Please remove/move this to rerun demo.\n")
        logfp.write("\n")
        sys.exit(2)

    rexclude = re.compile("##(fileDate|source|startTime|reference|cmdline)")

    def rstream(f) :
        """
        stream filtered lines from gzipped manta vcf
        """

        def rstreamFilter(line) :
            """
            predicate defining manta output lines to ignore
            """
            return (rexclude.match(line) is not None)

        rfp = gzip.open(f, 'rt')
        for line in rfp :
            if rstreamFilter(line) : continue
            yield line

    for (label, configArgs) in analyses :
        analysisDir=os.path.join(analysisBaseDir,label)
        expectedDir=os.path.join(expectedBaseDir,label)

        #
        # Step 1: configure demo
        #
        cmd  = "\"%s\"" % (configScript)
        for configArg in configArgs :
            cmd += " " + configArg
        cmd += " --runDir=\"%s\"" % (analysisDir)

        logfp.write("\n")
        logfp.write("**** Starting demo configuration and run: %s\n" % (label))
        logfp.write("**** Configuration cmd: %s\n" % (cmd))
        logfp.write("\n")

        config_retval=subprocess.call(cmd,shell=True)

        if config_retval != 0 :
            logfp.write("\n")
            logfp.write("ERROR: Demo configuration step failed\n")
            logfp.write("\n")
            sys.exit(1)
        else :
            logfp.write("\n")
            logfp.write("**** Completed demo configuration.\n")
            logfp.write("\n")

        #
        # Step 2: run demo (on single local core):
        #
        cmd=[sys.executable,os.path.join(analysisDir,"runWorkflow.py"),"-j","1","-g","4"]

        logfp.write("\n")
        logfp.write("**** Starting demo workflow execution.\n")
        logfp.write("**** Workflow cmd: '%s'\n" % (" ".join(cmd)))
        logfp.write("\n")

        run_retval=subprocess.call(cmd)

        if run_retval != 0 :
            logfp.write("\n")
            logfp.write("ERROR: Workflow execution step failed\n")
            logfp.write("\n")
            sys.exit(1)
        else :
            logfp.write("\n")
            logfp.write("**** Completed demo workflow execution.\n")
            logfp.write("\n")

        #
        # Step 3: Compare results to expected calls
        #
        resultsDir=os.path.join(analysisDir,"results","variants")
        logfp.write("\n")
        logfp.write("**** Starting comparison to expected results.\n")
        logfp.write("**** Expected results dir: '%s'\n" % (expectedDir))
        logfp.write("**** Demo results dir: '%s'\n" % (resultsDir))
        logfp.write("\n")

        expectedFiles = sorted([f for f in os.listdir(expectedDir) if os.path.isfile(os.path.join(expectedDir,f))])
        if len(expectedFiles) == 0 :
            logfp.write("\n")
            logfp.write("ERROR: No expected results found in '%s'\n" % (expectedDir))
            logfp.write("\n")
            sys.exit(1)

        for f in expectedFiles :
            efile=os.path.join(expectedDir,f)
            rfile=os.path.join(resultsDir,f)

            if not os.path.isfile(rfile) :
                logfp.write("\n")
                logfp.write("ERROR: Demo results file is missing: '%s'\n" % (rfile))
                logfp.write("\n")
                sys.exit(1)

            udiff_output = difflib.unified_diff(list(rstream(efile)),list(rstream(rfile)),fromfile=efile,tofile=rfile,n=0)
            udiff_str = "".join(list(udiff_output))

            if len(udiff_str) :
                logfp.write(udiff_str)

                logfp.write("\n")
                logfp.write("\n")
                logfp.write("ERROR: Found difference between demo and expected results in file '%s'.\n" % (f))
                logfp.write("       Expected file: '%s'\n" % (efile))
                logfp.write("       Demo results file: '%s'\n" % (rfile))
                logfp.write("\n")
                logfp.write("\n")

                sys.exit(1)

        logfp.write("\n")
        logfp.write("**** No differences between expected and computed results: %s\n" % (label))
        logfp.write("\n")

    logfp.write("\n")
    logfp.write("**** Demo/verification successfully completed\n")
    logfp.write("\n")



main()
