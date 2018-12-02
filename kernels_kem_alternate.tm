-----
MU69 Flyby Kernel Set, for Hazard Team use at MU69 ORT1 (January 2018).
-----

This metakernel file is a copy of NHGV's default kernel set for the NH MU69
Prime trajectory flyby. It includes all SPICE files used by NHGV, including
frames kernels, trajectories for s/c and many planets and satellites, MU69,
frames kernels, instrument kernels, etc. This is a larger set than necessary
for KEM (e.g., it includes kernels for the Pluto system and other planets), but
it is guaranteed to be complete.

All of these files exist on ixion, in GV's kernel respository at
/home/html/nh/gv-dev/kernels. For code run on ixion, this metakernel can be
incorporated using FURNSH() with nothing else necessary.

Henry Throop 12-Jan-2018

-----

\begindata

PATH_SYMBOLS = ('DIR', 'DIR2')

PATH_VALUES  = ('/Users/throop/gv/dev/kernels', '/Users/throop/git/NH_rings')

KERNELS_TO_LOAD = (
'$DIR/pck00009.tpc',
'$DIR/2014MU69_SBP_170528_may25a.bsp',
'$DIR/nh_v220.tf',
'$DIR/new-horizons_1646.tsc',
'$DIR/gv_pck.tpc',
'$DIR/gv_naif_ids.txt',
'$DIR/pluto_solar_heliographic.tf',
'$DIR/gv_pluto_smallsats.tf',
'$DIR/gv_smallbodies.tf',
'$DIR/plu013.bsp',
'$DIR/Nix_ephem_v02.bsp',
'$DIR/Hydra_ephem_v02.bsp'
'$DIR/nh_nom_20060119_20150727_v03.bsp',
'$DIR/nh_pred_20141201_20190301_od122.bsp',
'$DIR/nh_ref_20150801_20190901_od128_tcm22.bsp'
'$DIR/merged_nhpc_2007_v001.bc',
'$DIR/NavSE_plu047_od122.bsp',
'$DIR/NavPE_de433_od122.bsp',
'$DIR/nh_scispi_2015_pred.bc',
'$DIR/nh_scispi_2015_recon.bc',
'$DIR/nh_lorri_wcs.bc'
'$DIR/merged_18359_v1_cmd.bc'
'$DIR/NavPE_de433_od128.bsp',
'$DIR/2014MU69_SBP_180112_rd2b.bsp',
'$DIR/nh_nom_20170201_20210714_v13_alternate.bsp'
'$DIR/NavSBE_2014MU69_od135.bsp',
'$DIR/NavPE_de433_od135.bsp'
'$DIR/gv_pluto_smallsats_lhr.tf',
'$DIR/gv_pck.tpc',
'$DIR2/mu69_ring_frames.tf',
)
