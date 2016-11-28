
##
## Queries
##

# All galaxy properties use Henriques2015a, all halo properties use MRscPlanck1

# ---- Galaxy properties
# snapnum | z
# 16 | 6.42
# 19 | 5.03
# 22 | 3.95
# 25 | 3.10
# 30 | 2.07
#
# ---- Selection criteria
# stellar mass: M* > 10e9 Msol | stellarmass > 0.1
# star formation rate: SFR > 1 Msol yr^-1 | sfr > 1
# * use 1-e9 stellar mass selection to subset for 10e10 stellar mass galaxies


wget --http-user=***** --http-passwd=****** "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=
select
  zn.galaxyId as zn_galaxyId,
  zn.haloId as zn_haloId,
  zn.fofCentralId as zn_fofCentralId,
  zn.phkey as zn_phKey,
  zn.x as zn_x, zn.y as zn_y, zn.z as zn_z,
  zn.stellarMass as zn_stellarMass,
  zn.sfr as zn_sfr, zn.vvir as zn_vvir, zn.vmax as zn_vmax,
  zn.blackHoleMass as zn_blackHoleMass, zn.xrayLum as zn_xrayLum,
  zn.quasarAccretionRate as zn_quasarAccretionRate,
  zn.radioAccretionRate as zn_radioAccretionRate,
  zn.coldGas as zn_coldGas, zn.hotGas as zn_hotGas,
  zn.metalsHotGas as zn_metalsHotGas,
  z0.z0_haloId as z0_haloId,
  z0.z0_mcrit200 as z0_mcrit200,
  z0.z0_centralId as z0_centralId,
  z0.z0_central_mcrit200 as z0_central_mcrit200
from
  /* Select all galaxies above a given observational threshold for a particular redshift */
  (select
    galaxyId, haloId, fofCentralId, phkey, x, y, z,
    stellarMass, sfr, vvir, vmax,
    blackHoleMass, xrayLum, descendantId, treeRootID,
    quasarAccretionRate, radioAccretionRate, coldGas, hotGas, metalsHotGas
  from
    Henriques2015a..MRscPlanck1
  where
    snapnum = 30
    and sfr > 1) zn

  left join

  (select
   prog.haloId as zn_haloId,
   z0_join.z0_haloId as z0_haloId,
   z0_join.z0_mcrit200 as z0_mcrit200,
   z0_join.z0_central_haloId as z0_centralId,
   z0_join.z0_central_mcrit200 as z0_central_mcrit200
   from
      MPAHaloTrees..MRscPlanck1 prog,
      (select
         z0_all.m_crit200 as z0_mcrit200,
         z0_all.lastProgenitorId as z0_lastProgenitorId,
         z0_all.haloId as z0_haloId,
         z0_central.haloId as z0_central_haloId,
         z0_central.m_crit200 as z0_central_mcrit200
        from
            /* Select all halos at z = 0 */
            (select
              m_crit200, haloId, firstHaloInFOFgroupId,
              lastProgenitorId, mainLeafId
             from
              MPAHaloTrees..MRscPlanck1
             where
              snapnum = 58) z0_all,

            /* select all central halos at z = 0 */
            (select
              m_crit200, haloId,
              lastProgenitorId, mainLeafId
             from
              MPAHaloTrees..MRscPlanck1
             where
              snapnum = 58
              and haloID = firstHaloInFOFgroupId) z0_central
         where
           /* Match centrals  */
           z0_central.haloId = z0_all.firstHaloInFOFgroupId) z0_join

       where
        prog.haloId between z0_join.z0_haloId and z0_join.z0_lastprogenitorId
        and prog.snapnum = 30
     ) z0

      on zn.haloId = z0.zn_haloId
" -O henriques2015a_z2p07_sfr.csv

"
