
## wget query

wget --http-user=sussex --http-passwd=G787739L "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=" -O temp.csv

-- z | snapnum | m_crit200
-- 9.72 | 11 | 5
-- 8.93 | 12 |
-- 8.22 | 13 |
-- 6.97 | 15 | 28
-- 5.91 | 17 | 57
-- 5.03 | 19 | 105
-- 3.95 | 22 | 250
-- 3.10 | 25 | 500
-- 2.07 | 30 | 1300


select
  zn.haloId as zn_haloId,
  zn.m_crit200 as zn_m_crit200,
  z0.haloId as z0_haloId,
  z0.lastProgenitorId as z0_lastProgenitorId,
  z0.m_crit200 as z0_m_crit200
from

  (select
    m_crit200, haloId
  from
    MPAHaloTrees..MRscPlanck1
  where
    snapnum = 11
    and haloID = firstHaloInFOFgroupId
    and m_crit200 > 5) zn,

  (select
    m_crit200, haloId, lastProgenitorId
   from
    MPAHaloTrees..MRscPlanck1
   where
    snapnum = 58) z0

where z0.haloId < zn.haloId
and z0.lastProgenitorId >= zn.haloId




 ## select largest halos at z = n, find z = 0 descendants
 # selects largest central halos at z = n, then finds all central z = n progenitors of z = 0 halos, does a left join between each of the tables


select
  zn.m_crit200 as zn_m_crit200,
  zn.haloId as zn_haloId,
  z0.z0_m_crit200 as z0_m_crit200,
  z0.z0_haloId as z0_haloId,
  z0.z0_lastProgenitorId as z0_lastProgenitorId
from

  (select
    m_crit200, haloId
  from
    MPAHaloTrees..MRscPlanck1
  where
    snapnum = 15
    and haloID = firstHaloInFOFgroupId
    and m_crit200 > 28) zn

left join
  (select
    zn_temp.haloId as zn_haloId,
    z0_temp.haloId as z0_haloId,
    z0_temp.lastProgenitorId as z0_lastProgenitorId,
    z0_temp.m_crit200 as z0_m_crit200
  from
    (select
      m_crit200, haloId
    from
      MPAHaloTrees..MRscPlanck1
    where
      snapnum = 15
      and haloID = firstHaloInFOFgroupId) zn_temp,

    (select
      m_crit200, haloId, lastProgenitorId
     from
      MPAHaloTrees..MRscPlanck1
     where
      snapnum = 58) z0_temp

   where zn_temp.haloId > z0_temp.haloId
   and zn_temp.haloId <= z0_temp.lastProgenitorId) z0

on zn.haloId = z0.zn_haloId



## select largest halos at z = 0, find z = n descendants
# can be used to subset by largest z = n descendant in post-processing
# can be used to select z = n main branch descendants (via mainLeafId)

select
  zn.haloId as zn_haloId,
  zn.m_crit200 as zn_m_crit200,
  z0.haloId as z0_haloId,
  z0.m_crit200 as z0_m_crit200,
  z0.mainLeafId as z0_mainLeafId
from
  (select
    m_crit200,
    haloId
  from
    MPAHaloTrees..MRscPlanck1
  where
    snapnum = 30
    and haloID = firstHaloInFOFgroupId) zn,

  (select
    m_crit200,
    haloId,
    lastProgenitorId,
    mainLeafId
   from
    MPAHaloTrees..MRscPlanck1
   where
    snapnum = 58
    and haloID = firstHaloInFOFgroupId
    and m_crit200 > 1e4) z0

 where zn.haloId > z0.haloId
 and zn.haloId <= z0.lastProgenitorId


## find all protocluster halos at a given redshift

wget --http-user=sussex --http-passwd=G787739L "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=
select
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
        POWER((m_crit200 * 3.)/(200 * 12.7 * 4 * Pi()), 1./3) * 0.67 as r_crit200
        from mpahalotrees..mrscplanck1
        where snapnum = 58
        and firstHaloInFOFgroupId = haloId
        and m_crit200 > 1e4) cen,

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
    and prog.snapnum = 30
" -O protocluster_halos_r200_z2p07.csv

"
