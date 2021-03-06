# dunetpc/dune/Protodune/singlephase/Tool/CMakeLists.txt

# David Adams
# June 2018

# Add plugin for each tool.

simple_plugin(PdspOnlineChannel   "tool"
                dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                art_Utilities canvas
                cetlib cetlib_except
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(IcebergOnlineChannel   "tool"
                dune-raw-data_Services_ChannelMap_IcebergChannelMapService_service
                art_Utilities canvas
                cetlib cetlib_except
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ProtoduneOnlineChannel   "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ProtoDuneChannelRanges  "tool"
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ApaChannelRanges  "tool"
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(HalfApaChannelRanges  "tool"
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ProtoDuneChannelGroups  "tool"
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ApaChannelGroups  "tool"
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
