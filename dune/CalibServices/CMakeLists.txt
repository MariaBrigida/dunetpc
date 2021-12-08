cet_build_plugin(XYZCalibServiceProtoDUNE art::service
              dune_Calib
	      ART_FRAMEWORK_CORE
	      ART_FRAMEWORK_PRINCIPAL
	      ART_FRAMEWORK_SERVICES_REGISTRY
	      art::Persistency_Common canvas
	      art::Persistency_Provenance canvas
	      ART_ROOT_IO_TFILE_SUPPORT
	      ROOT::Core
	      ART_ROOT_IO_TFILESERVICE_SERVICE
	      ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE
	      art::Utilities canvas
	      messagefacility::MF_MessageLogger
	      
  	      cetlib::cetlib cetlib_except
   	      ROOT_BASIC_LIB_LIST
              
)

cet_build_plugin(LifetimeCalibServiceProtoDUNE art::service
      dune_Calib
      ART_FRAMEWORK_CORE
      ART_FRAMEWORK_PRINCIPAL
      ART_FRAMEWORK_SERVICES_REGISTRY
      art::Persistency_Common canvas
      art::Persistency_Provenance canvas
      ART_ROOT_IO_TFILE_SUPPORT
      ROOT::Core
      ART_ROOT_IO_TFILESERVICE_SERVICE
      ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE
      art::Utilities canvas
      messagefacility::MF_MessageLogger

      cetlib::cetlib cetlib_except
      ROOT_BASIC_LIB_LIST

) 


install_headers()
install_fhicl()
install_source()
