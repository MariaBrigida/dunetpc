# dunetpc/dune/DAQTriggerSim/Service/CMakeLists.txt
#
# Michael Baird
# March 2016
#
# CMake instructions for DAQTriggerSim/Services

simple_plugin(TemplateTriggerService   "service"
                ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                nurandom_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                art_Utilities canvas
                ${MF_MESSAGELOGGER}
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(HardwareMapperService   "service"
                larcorealg_Geometry
                dune_Geometry
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                ${ART_PERSISTENCY_PROVENANCE}
                ${MF_MESSAGELOGGER}
                ${ROOT_BASIC_LIB_LIST}
                art_Utilities canvas
             )

install_fhicl()
install_headers()
install_source()
