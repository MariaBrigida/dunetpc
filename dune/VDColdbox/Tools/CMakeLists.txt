cet_build_plugin(VDColdboxChannelRanges  art::tool
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )
cet_build_plugin(VDColdboxChannelGroups  art::tool
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )


install_fhicl()
