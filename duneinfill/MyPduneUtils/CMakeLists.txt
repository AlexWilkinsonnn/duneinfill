
art_make( LIBRARY_NAME MyProtoDUNEUtilities
	  LIB_LIBRARIES
                        dunecore::ArtSupport
                        lardataobj_RawData
                        lardataobj_RecoBase
                        larreco_RecoAlg
                        larreco_Calorimetry
                        larreco_RecoAlg_PMAlg
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardataalg_DetectorInfo
                        lardata_ArtDataHelper
                        lardata_Utilities
                        larsim_Simulation 
                        lardataobj_Simulation
                        duneprototypes::ProtoDUNEDataUtils
                        geant4reweight::ReweightBaseLib
                        geant4reweight::PropBaseLib
                        dunecore::DuneObj
                        dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        larsim_MCCheater_PhotonBackTrackerService_service
                        larsim_MCCheater_BackTrackerService_service
                        larsim_MCCheater_ParticleInventoryService_service
                        nug4_ParticleNavigation 
                        nusimdata_SimulationBase
                        ${CLHEP}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ROOT_BASIC_LIB_LIST
                        ProtoDUNEUtilities
                        
	  MODULE_LIBRARIES
          MyProtoDUNEUtilities
          ProtoDUNEUtilities
)

install_headers()
install_fhicl()
install_source()
install_scripts()
