#!/bin/bash

snapnum=(
11
12
13
15
17
19
22
25
30
)

z=(
9p72
8p93
8p22
6p97
5p92
5p03
3p95
3p10
2p07
)

echo "starting download..."


for index in ${!snapnum[*]}; do

  echo "${snapnum[$index]}, ${z[$index]}"

  wget --http-user=sussex --http-passwd=G787739L "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=
  select top 10000
          m_crit200, haloId
          from
            MPAHaloTrees..MRscPlanck1
          where
            snapnum = ${snapnum[$index]}
            and haloID = firstHaloInFOFgroupId
        order by m_crit200 desc

  " -O data/halos/z${z[$index]}_halos.csv -q

done
