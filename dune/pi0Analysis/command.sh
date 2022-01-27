#samweb def: prodgenie_nue_dune10kt_1x2x6_mcc11_lbl_reco
#lar -c run_pi0Analysis.fcl /pnfs/dune/tape_backed/dunepro/mcc11/protodune/mc/full-reconstructed/08/54/87/06/nue_dune10kt_1x2x6_13428997_0_20181124T053021_gen_g4_detsim_reco.root -n 1

#lar -c run_pi0Analysis.fcl prodgenie_nue_dune10kt_1x2x6_gen_g4_detsim_reco.root
#samweb definitions: https://dune-data.fnal.gov/mc/mcc11/index.html, FD Beamsim Requests
#e.g. prodgenie_anu_dune10kt_1x2x6_mcc11_lbl_reco



#lar -c run_pi0Analysis.fcl prodgenie_nue_dune10kt_1x2x6_gen_g4_detsim_reco.root -n 10
#lar -c run_pi0AnalysisPlots.fcl 


#lar -c run_pi0Analysis.fcl  -S outputList.txt -n 100		#these are mcc11 files
lar -c run_pi0Analysis.fcl  -S good_nue_files.txt -n 1000	#these are dunetpc v09_36_00_01 files I've made


#debugging
#lar -c run_pi0Analysis.fcl /pnfs/dune/scratch/users/mbrunett/FDsamples/files/v09_36_00_01/nue_dunefd_1x2x6_20k/geng4detsimreco/50424407_1281/RootOutput-292d-bc6b-51a5-cbf9_c2786763-c6ef-4f0a-923d-b45f25b1c816.root

#Samples
#../../../../../pi0Samples/nuesample_100K_standard/pi0AnalysisOutput_nue_100K_standard.root
