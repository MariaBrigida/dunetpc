#samweb def: prodgenie_nue_dune10kt_1x2x6_mcc11_lbl_reco
#lar -c run_pi0Analysis.fcl /pnfs/dune/tape_backed/dunepro/mcc11/protodune/mc/full-reconstructed/08/54/87/06/nue_dune10kt_1x2x6_13428997_0_20181124T053021_gen_g4_detsim_reco.root -n 1

#lar -c run_pi0Analysis.fcl prodgenie_nue_dune10kt_1x2x6_gen_g4_detsim_reco.root


lar -c run_pi0Analysis.fcl  -S outputList.txt -n 10
