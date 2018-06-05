#!/bin/bash
APPDIR=/ccc/cont007/dsku/lautrec/home/app/fg/fg/products

BAM=/ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_raw/MYBAMID.reliable.realign.recal.bam
$APPDIR/samtools-1.6/el7-x86_64-haswell/bin/samtools view -h ${BAM} | $APPDIR/bwakit-0.7.15/bin/k8 $APPDIR/bwakit-0.7.15/bin/bwa-postalt.js $APPDIR/bwakit-0.7.15/resource-GRCh38/hs38DH.fa.alt | \
$APPDIR/sambamba-0.6.5/el7-x86_64-generic/bin/sambamba view -S -f bam -l 0 /dev/stdin | \
$APPDIR/sambamba-0.6.5/el7-x86_64-generic/bin/sambamba sort -t 50 -m 6G --tmpdir=/ccc/scratch/cont007/fg0094/soudadel/MESO_tmp -o /ccc/store/cont007/fg0094/soudadel/MESO/post_align/MYBAMID_pa.bam /dev/stdin
