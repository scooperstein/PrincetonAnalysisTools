#!/bin/bash

## Resubmit failed analysis jobs to condor. To run on FNAL LPC, do:
## bash resubmit_failed_jobs.sh [jobdir]
##
## Author: Stephane Cooperstein
##

echo $1
for sampledir in $1/*
    do
        echo $sampledir
        for i in `seq 0 500`
            do
                if [ -f $sampledir/$i.submit ]; then
                    if [ ! -f $sampledir/output_*_$i.root ]; then
                        echo "output file not found, resubmitting $sampledir/$i.submit"
                        echo "condor_submit $sampledir/$i.submit"
                        condor_submit $sampledir/$i.submit
                    fi
                fi      
            done
    done
