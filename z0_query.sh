#!/bin/bash

out_dir=/lustre/scratch/astro/cl478/protoclusters_data/
username=$1
password=$2

wget --http-user=${username} --http-passwd=${password} "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=
select
hen.x as hen_x, hen.y as hen_y, hen.z as hen_z,
hen.subHaloId as hen_subhaloId,
hen.fofSubHaloId as hen_fofSubHaloId,
hen.stellarMass as hen_stellarMass, hen.sfr as hen_sfr, hen.snapnum as hen_snapnum,
hen.centralMvir as hen_centralMvir, hen.centralRvir as hen_centralRvir,
cen.haloId as cen_haloId, cen.firstHaloInFOFgroupId as cen_firstHaloInFOFgroupId,
cen.x as cen_x, cen.y as cen_y, cen.z as cen_z,
cen.m_crit200 as cen_m_crit200, cen.r_crit200 as cen_r_crit200
from

  (select
      x, y, z,  galaxyId, subHaloId, fofSubHaloId, stellarMass,
      sfr, snapnum, centralMvir, centralRvir
      from henriques2015a..MRscPlanck1
      where
      snapnum = 58 and
      stellarMass >= 0.1) hen,

  (select
      m_crit200,
      POWER((m_crit200 * 3.)/(200 * 27.7619760963 * 4 * Pi()), 1./3) as r_crit200,
      haloId, firstHaloInFOFgroupId,
      x, y, z
      from
      mpahalotrees..mrscplanck1
      where snapnum = 58
      and firstHaloInFOFgroupId = haloId
      and m_crit200 > 6.73e3) cen

where
hen.x between cen.x - cen.r_crit200 and cen.x %2B cen.r_crit200
and hen.y between cen.y - cen.r_crit200 and cen.y %2B cen.r_crit200
and hen.z between cen.z - cen.r_crit200 and cen.z %2B cen.r_crit200
and POWER(POWER(hen.x - cen.x,2) %2B POWER(hen.y - cen.y,2) %2B POWER(hen.z - cen.z,2),0.5) <= cen.r_crit200
" -O ${out_dir}z0_gals.csv -q
