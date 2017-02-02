#!/bin/bash

selection=(
stellarMass
sfr
)

selectionValue=(
0.1
1
)

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


for j in ${!snapnum[*]}; do
  for index in ${!snapnum[*]}; do

    echo "${snapnum[$index]}, ${z[$index]} | ${selection[$j]} > ${selectionValue[$j]}"

    wget --http-user=sussex --http-passwd=G787739L "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=
    select
      zn.galaxyId as zn_galaxyId,
      zn.haloId as zn_haloId,
      zn.fofCentralId as zn_fofCentralId,
      zn.x as zn_x, zn.y as zn_y, zn.z as zn_z,
      zn.velX as zn_velX, zn.velY as zn_velY, zn.velZ as zn_velZ,
      zn.stellarMass as zn_stellarMass,
      zn.sfr as zn_sfr,
      zn.blackHoleMass as zn_blackHoleMass,
      zn.quasarAccretionRate as zn_quasarAccretionRate,
      zn.radioAccretionRate as zn_radioAccretionRate,
      zn.np as zn_np, zn.centralMvir as zn_centralMvir,
      zn.mvir as zn_mvir, zn.rvir as zn_rvir,
      z0.z0_haloId as z0_haloId,
      z0.z0_mcrit200 as z0_mcrit200,
      z0.z0_centralId as z0_centralId,
      z0.z0_central_mcrit200 as z0_central_mcrit200,
      z0.z0_central_rcrit200 as z0_central_rcrit200,
      z0.z0_x as z0_x, z0.z0_y as z0_y, z0.z0_z as z0_z,
      z0.z0_central_x as z0_central_x,
      z0.z0_central_y as z0_central_y,
      z0.z0_central_z as z0_central_z
    from
      /* Select all galaxies above a given observational threshold for a particular redshift */
      (select
        galaxyId, haloId, fofCentralId, x, y, z,
        velX, velY, velZ, stellarMass, sfr,
        blackHoleMass, descendantId, treeRootID,
        quasarAccretionRate, radioAccretionRate,
        np, centralMvir, mvir, rvir
      from
        Henriques2015a..MRscPlanck1
      where
        snapnum = ${snapnum[$index]}
        and ${selection[$j]} > ${selectionValue[$j]}) zn

      left join

      (select
       prog.haloId as zn_haloId,
       z0_join.halos_haloId as z0_haloId,
       z0_join.halos_x as z0_x, z0_join.halos_y as z0_y, z0_join.halos_z as z0_z,
       z0_join.halos_m_crit200 as z0_mcrit200,
       z0_join.cen_x as z0_central_x, z0_join.cen_y as z0_central_y, z0_join.cen_z as z0_central_z,
       z0_join.cen_haloId as z0_centralId,
       z0_join.cen_m_crit200 as z0_central_mcrit200,
       z0_join.cen_r_crit200 as z0_central_rcrit200
       from
          MPAHaloTrees..MRscPlanck1 prog,
          (select
          cen.haloId as cen_haloId,
          cen.x as cen_x, cen.y as cen_y, cen.z as cen_z,
          cen.m_crit200 as cen_m_crit200, cen.r_crit200 as cen_r_crit200,
          halos.x as halos_x, halos.y as halos_y, halos.z as halos_z,
          halos.haloId as halos_haloId,
          halos.m_crit200 as halos_m_crit200,
          halos.lastProgenitorId as halos_lastProgenitorId

          from

            (select
            m_crit200, haloId, x, y, z,
            POWER((m_crit200 * 3.)/(200 * 27.7619760963 * 4 * Pi()), 1./3) * 0.67 as r_crit200
            from mpahalotrees..mrscplanck1
            where snapnum = 58
            and firstHaloInFOFgroupId = haloId) cen,

            (select
            x, y, z, haloId, m_crit200, lastProgenitorId, firstHaloInFOFgroupId
            from mpahalotrees..mrscplanck1
            where snapnum = 58) halos

          where
          cen.haloId = halos.firstHaloInFOFgroupId
          and halos.x between cen.x - cen.r_crit200 and cen.x %2B cen.r_crit200
          and halos.y between cen.y - cen.r_crit200 and cen.y %2B cen.r_crit200
          and halos.z between cen.z - cen.r_crit200 and cen.z %2B cen.r_crit200
          and POWER(POWER(halos.x - cen.x, 2) %2B POWER(halos.y - cen.y, 2) %2B POWER(halos.z - cen.z, 2), 0.5) <= cen.r_crit200) z0_join

        where
        prog.haloId between z0_join.halos_haloId and z0_join.halos_lastProgenitorId
        and prog.snapnum = ${snapnum[$index]}) z0

        on zn.haloId = z0.zn_haloId
    " -O data/r200/henriques2015a_z${z[$index]}_sfr_r200.csv -q

  done
done
