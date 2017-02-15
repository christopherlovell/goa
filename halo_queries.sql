
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
    snapnum = 11
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
