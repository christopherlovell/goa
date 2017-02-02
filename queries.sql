
##
## Queries
##

# All galaxy properties use Henriques2015a, all halo properties use MRscPlanck1

# ---- Galaxy properties
# snapnum | z
# 9 | 11.51
# 10 | 10.57
# 11 | 9.72
# 12 | 8.93
# 13 | 8.22
# 15 | 6.97
# 17 | 5.92
# 19 | 5.03
# 22 | 3.95
# 25 | 3.10
# 30 | 2.07
#
# ---- Selection criteria
# stellar mass: M* > 10e9 Msol | stellarmass > 0.1
# star formation rate: SFR > 1 Msol yr^-1 | sfr > 1
# * use 1-e9 stellar mass selection to subset for 10e10 stellar mass galaxies


### Clusters defined within r200

-- zn.coldGas as zn_coldGas, zn.hotGas as zn_hotGas,
-- zn.metalsHotGas as zn_metalsHotGas,
-- zn.vvir as zn_vvir, zn.vmax as zn_vmax,
-- zn.xrayLum as zn_xrayLum,
-- zn.phkey as zn_phKey,



wget --http-user=***** --http-passwd=******* "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=
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
    snapnum = 13
    and stellarMass > 0.1) zn

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
        POWER((m_crit200 * 3.)/(200 * 12.7 * 4 * Pi()), 1./3) * 0.67 as r_crit200
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
    and prog.snapnum = 13) z0

    on zn.haloId = z0.zn_haloId
" -O henriques2015a_z8p22_mstar_r200.csv

"




### Clusters defined as the entire FOF group at z = 0

-- wget --http-user=***** --http-passwd=****** "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=
-- select
--   zn.galaxyId as zn_galaxyId,
--   zn.haloId as zn_haloId,
--   zn.fofCentralId as zn_fofCentralId,
--   zn.phkey as zn_phKey,
--   zn.x as zn_x, zn.y as zn_y, zn.z as zn_z,
--   zn.stellarMass as zn_stellarMass,
--   zn.sfr as zn_sfr, zn.vvir as zn_vvir, zn.vmax as zn_vmax,
--   zn.blackHoleMass as zn_blackHoleMass, zn.xrayLum as zn_xrayLum,
--   zn.quasarAccretionRate as zn_quasarAccretionRate,
--   zn.radioAccretionRate as zn_radioAccretionRate,
--   zn.coldGas as zn_coldGas, zn.hotGas as zn_hotGas,
--   zn.metalsHotGas as zn_metalsHotGas,
--   z0.z0_haloId as z0_haloId,
--   z0.z0_mcrit200 as z0_mcrit200,
--   z0.z0_x as z0_x, z0.z0_y as z0_y, z0.z0_z as z0_z,
--   z0.z0_centralId as z0_centralId,
--   z0.z0_central_mcrit200 as z0_central_mcrit200
-- from
--   /* Select all galaxies above a given observational threshold for a particular redshift */
--   (select
--     galaxyId, haloId, fofCentralId, phkey, x, y, z,
--     stellarMass, sfr, vvir, vmax,
--     blackHoleMass, xrayLum, descendantId, treeRootID,
--     quasarAccretionRate, radioAccretionRate, coldGas, hotGas, metalsHotGas
--   from
--     Henriques2015a..MRscPlanck1
--   where
--     snapnum = 30
--     and sfr > 1) zn
--
--   left join
--
--   (select
--    prog.haloId as zn_haloId,
--    z0_join.z0_haloId as z0_haloId,
--    z0_join.z0_mcrit200 as z0_mcrit200,
--    z0_join.z0_x as z0_x, z0_join.z0_y as z0_y, z0_join.z0_z as z0_z,
--    z0_join.z0_central_haloId as z0_centralId,
--    z0_join.z0_central_mcrit200 as z0_central_mcrit200
--    from
--       MPAHaloTrees..MRscPlanck1 prog,
--       (select
--          z0_all.m_crit200 as z0_mcrit200,
--          z0_all.lastProgenitorId as z0_lastProgenitorId,
--          z0_all.haloId as z0_haloId,
--          z0_all.x as z0_x, z0_all.y as z0_y, z0_all.z as z0_z,
--          z0_central.haloId as z0_central_haloId,
--          z0_central.m_crit200 as z0_central_mcrit200
--         from
--             /* Select all halos at z = 0 */
--             (select
--               m_crit200, haloId, firstHaloInFOFgroupId,
--               lastProgenitorId, mainLeafId, x, y, z
--              from
--               MPAHaloTrees..MRscPlanck1
--              where
--               snapnum = 58) z0_all,
--
--             /* select all central halos at z = 0 */
--             (select
--               m_crit200, haloId,
--               lastProgenitorId, mainLeafId
--              from
--               MPAHaloTrees..MRscPlanck1
--              where
--               snapnum = 58
--               and haloID = firstHaloInFOFgroupId) z0_central
--          where
--            /* Match centrals  */
--            z0_central.haloId = z0_all.firstHaloInFOFgroupId) z0_join
--
--        where
--         prog.haloId between z0_join.z0_haloId and z0_join.z0_lastprogenitorId
--         and prog.snapnum = 30
--      ) z0
--
--       on zn.haloId = z0.zn_haloId
-- " -O henriques2015a_z2p07_sfr.csv
--
-- "
