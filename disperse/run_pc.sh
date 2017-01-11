#!/bin/bash   

outdir="../data/disperse/out/"

for file in ../data/disperse/*.ascii
do
    ~/src/disperse/bin/delaunay_3D $file -outDir $outdir
    ~/src/disperse/bin/mse file $outdir$(basename $file)'.NDnet' -upSkl -nsig 5 -robustness -outDir $outdir

    ~/src/disperse/bin/skelconv $outdir$(basename $file)'.NDnet_s5.up.NDskl' -smooth 5 -breakdown -trimBelow robustness 0.4 -assemble 70 -to segs_ascii -outDir $outdir
    ~/src/disperse/bin/skelconv $outdir$(basename $file)'.NDnet_s5.up.NDskl' -smooth 5 -breakdown -trimBelow robustness 0.4 -assemble 70 -to crits_ascii -outDir $outdir
    ~/src/disperse/bin/skelconv $outdir$(basename $file)'.NDnet_s5.up.NDskl' -smooth 5 -breakdown -trimBelow robustness 0.4 -assemble 70 -to vtp -outDir $outdir
done

