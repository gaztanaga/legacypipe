#!/bin/bash
export ver="2"
export outdir="/scratch2/scratchdirs/kaylanb/dr3/runs/mzls_v${ver}"
export dr="/scratch2/scratchdirs/kaylanb/dr3/legacypipe/py/in_progress"
mkdir -p $dr $outdir
for brick in `cat brick_list.txt`;do
    export brick="$brick"
    bri=$(echo $brick | head -c 3)
    tractor_fits=$outdir/tractor/$bri/$brick/tractor-$brick.fits
    if [ -e "$tractor_fits" ]; 
    then
        echo $brick >> $dr/finished.txt
        echo skipping $brick, its done
    elif [ -e "$dr/inq_$brick.txt" ] && [ ! -e "$dr/ran_$brick.txt" ]
    then
        echo skipping $brick, its queued and has not ran
    elif [ ! -e "$dr/inq_$brick.txt" ] 
    then        
        echo submitting $brick
        sbatch submit_dr4.sh --export brick,ver,outdir,dr
        rm $dr/ran_$brick.txt
        touch $dr/inq_$brick.txt
    else
        echo AHHHHHHH should not reach this point!!!!!!!!!!!!!!!
    fi
done
