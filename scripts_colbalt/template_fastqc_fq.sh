#!/bin/bash

module load extenv/fg

APPDIR=/ccc/cont007/dsku/lautrec/home/app/fg/fg/products

$APPDIR/fastqc-0.11.3/fastqc -t 12 -o /ccc/store/cont007/fg0094/soudadel/MESO/fastqc -j java $FQFILE
