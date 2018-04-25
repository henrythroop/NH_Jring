\begindata

KERNELS_TO_LOAD = (
     'kernels/pck00010.tpc',
     'kernels/de418.bsp',

     'kernels/naif0012.tls',
     'kernels/new-horizons_1121.tsc',

     'kernels/nh_nom_20060119_20150727_v03.bsp',
     'kernels/nh_recon_e2j_v1.bsp'
     'kernels/jup310.bsp',

      )
\begintext


Kernels listed later have higher priority.

-----

This is a set of kernels just for the Jupiter flyby. It is minimal, so as to prevent any errors.
I need: 
  - Solar system kernel
  - Jupiter system kernel including Adrastea
  - NH .bsp kernel
  
jup310: Jupiter, Metis, Adrastea -- updated 2013 R Jacobsen. Replaces old
jup204 + jup260. [Changed this 25-Apr-2018]

de413:  Barycenters incl Jupiter bary. ** Brian C is using de418. Doesn't make a difference?
pck09:  pck10 is newer. Use it now.
naif0012: Latest
new-horizons_0664: Very old, from 2010. But still new enough not to matter. Updating to 1121 doesn't change things.

In Carcich kernel:


93230080	kernels/jup204.bsp
210829312	kernels/jup260.bsp
16456704	/home/html/nh/science_spice/kernels/spk/de418.bsp
210833408	/home/html/nh/science_spice/kernels/spk/jup260.bsp

2531328	/home/html/nh/science_spice/kernels/spk/nh_recon_e2j_v1.bsp *** This is the key one I was missing
359424	/home/html/nh/science_spice/kernels/spk/nh_recon_j2sep07_prelimv1.bsp
14873600	/home/html/nh/science_spice/kernels/spk/nh_recon_od077_v01.bsp
2447360	/home/html/nh/science_spice/kernels/spk/nh_recon_od117_v01.bsp
1724416	/home/html/nh/science_spice/kernels/spk/nh_recon_pluto_od122_v01.bsp
304297984	/home/html/nh/science_spice/kernels/ck/merged_nhpc_2007_v006.bc
3804160	/home/html/nh/science_spice/kernels/spk/nh_pred_20141201_20190301_od122.bsp
2156544	/home/html/nh/science_spice/kernels/spk/NavSE_plu047_od122.bsp
3385344	/home/html/nh/science_spice/kernels/spk/NavPE_de433_od122.bsp
162402304	/home/html/nh/science_spice/kernels/ck/nh_scispi_2015_pred.bc
104508416	/home/html/nh/science_spice/kernels/ck/nh_scispi_2015_recon.bc
465920	/home/html/nh/science_spice/kernels/ck/nh_lorri_wcs.bc
96345	/home/html/nh/science_spice/kernels/sclk/new-horizons_1121.tsc
5257	/home/html/nh/science_spice/kernels/lsk/naif0012.tls
12555	/home/html/nh/science_spice/kernels/pck/nh_targets_v001.tpc
126143	/home/html/nh/science_spice/kernels/pck/pck00010.tpc
4926	/home/html/nh/science_spice/kernels/pck/nh_pcnh_006.tpc
128248	/home/html/nh/science_spice/kernels/fk/nh_v220.tf
1698	/home/html/nh/science_spice/kernels/ik/nh_allinstruments_v002.ti
45993	/home/html/nh/science_spice/kernels/ik/nh_alice_v200.ti
54270	/home/html/nh/science_spice/kernels/ik/nh_lorri_v201.ti
35778	/home/html/nh/science_spice/kernels/ik/nh_pepssi_v110.ti
30607	/home/html/nh/science_spice/kernels/ik/nh_ralph_v100.ti
9036	/home/html/nh/science_spice/kernels/ik/nh_rex_v100.ti
7681	/home/html/nh/science_spice/kernels/ik/nh_sdc_v101.ti
17367	/home/html/nh/science_spice/kernels/ik/nh_swap_v100.ti
7413	/home/html/nh/science_spice/kernels/ik/nh_astr_v000.ti
6291	/home/html/nh/science_spice/kernels/ik/nh_fss_v000.ti
5865	/home/html/nh/science_spice/kernels/fk/nh_soc_misc_v001.tf
19456	/home/html/nh/science_spice/kernels/spk/nh_stars.bsp
4419	kernels/gv_pluto_smallsats_lhr.tf
6751	kernels/gv_pck.tpc
