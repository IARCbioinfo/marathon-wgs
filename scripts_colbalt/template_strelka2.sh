#!/usr/bin/env bash
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2017 Illumina, Inc.
#
# Execute small somatic variant calling

# Template to be launched with launch_strelka2.sh (to replace TUMORBAM and NORMALBAM)

module load extenv/fg
module unload samtools
module load strelka/2.8.4
set -o nounset
set -o pipefail

TUMORBAM=/ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_post_al/tumor/TUMORBAMID.pa.bam
NORMALBAM=/ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_post_al/normal/symlinks/TUMORBAMID.normal.pa.bam
REFFILES=/ccc/work/cont007/fg0094/soudadel/MESO/calling_somatic/ref_files

#
# Step 1: configure
#

cmd="configureStrelkaSomaticWorkflow.py \
--tumorBam='$TUMORBAM' \
--normalBam='$NORMALBAM' \
--referenceFasta='$FG_BIOBANK/by-name/Homo_sapiens/hs38dh/hs38dh_all_chr.fasta' \
--callRegions='$REFFILES/hs38dh.bed.gz' \
--callMemMb=1024 \
--runDir=STRELKA/TUMORBAMID"

echo 1>&2
echo "**** Starting configuration and run." 1>&2
echo "**** Configuration cmd: '$cmd'" 1>&2
echo 1>&2
eval $cmd

if [ $? -ne 0 ]; then
    echo 1>&2
    echo "ERROR: configuration step failed" 1>&2
    echo 1>&2
    exit 1
else
    echo 1>&2
    echo "**** Completed configuration." 1>&2
    echo 1>&2
fi


#
# Step 2: run demo (on single local core):
#

cmd="STRELKA/TUMORBAMID/runWorkflow.py -m local -j 28"
echo 1>&2
echo "**** Starting workflow execution." 1>&2
echo "**** Workflow cmd: '$cmd'" 1>&2
echo 1>&2
$cmd
