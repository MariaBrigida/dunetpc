art_make(
          MODULE_LIBRARIES
                        larcoreobj_SummaryData
                        larcore_Geometry_Geometry_service
			larcorealg_Geometry
                        nurandom_RandomUtils_NuRandomService_service
                        larsim_Simulation lardataobj_Simulation
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Utilities
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${SQLITE3}
                        ${ROOT_EG}
                        ${IFDH}
        )

install_headers()
install_fhicl()
install_source()
