/**
 *  @file   PandoraPFANew/src/PandoraPFANewProcessor.cc
 * 
 *  @brief  Implementation of the pandora pfa new processor class.
 * 
 *  $Log: $
 */

#include "marlin/Exceptions.h"

#include "Api/PandoraApi.h"

#include "PandoraPFANewProcessor.h"

PandoraPFANewProcessor pandoraPFANewProcessor;

pandora::Pandora *PandoraPFANewProcessor::m_pPandora = NULL;

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraPFANewProcessor::PandoraPFANewProcessor() :
    Processor("PandoraPFANewProcessor"),
    m_nRun(0),
    m_nEvent(0)
{
    _description = "Pandora reconstructs clusters and particle flow objects";
    this->ProcessSteeringFile();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::init()
{
    try
    {
        streamlog_out(MESSAGE) << "PandoraPFANewProcessor - Init" << std::endl;
        this->FinaliseSteeringParameters();
        m_pPandora = new pandora::Pandora();

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_geometryCreator.CreateGeometry());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RegisterUserComponents());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPandora, m_settings.m_pandoraSettingsXmlFile));
    }
    catch (StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "Failed to initialize pandora pfa new processor: " << statusCodeException.ToString() << std::endl;
        throw;
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Failed to initialize pandora pfa new processor, unrecognized exception" << std::endl;
        throw;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::processRunHeader(LCRunHeader *pLCRunHeader)
{
    m_detectorName = pLCRunHeader->getDetectorName();
    streamlog_out(MESSAGE) << "Detector Name " << m_detectorName << ", Run " << ++m_nRun <<  std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::processEvent(LCEvent *pLCEvent)
{
    static int eventCounter = 0;

    if (eventCounter < m_settings.m_nEventsToSkip)
    {
        ++eventCounter;
        throw marlin::SkipEventException(this);
    }

    try
    {
        streamlog_out(MESSAGE) << "PandoraPFANewProcessor, Run " << m_nRun << ", Event " << ++m_nEvent << std::endl;

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_mcParticleCreator.CreateMCParticles(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_trackCreator.CreateTrackAssociations(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_trackCreator.CreateTracks(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_mcParticleCreator.CreateTrackToMCParticleRelationships(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_caloHitCreator.CreateCaloHits(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_mcParticleCreator.CreateCaloHitToMCParticleRelationships(pLCEvent));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPandora));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pfoCreator.CreateParticleFlowObjects(pLCEvent));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPandora));
        this->Reset();
    }
    catch (StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "StatusCodeException: " << statusCodeException.ToString() << std::endl;
        throw;
    }
    catch (lcio::EventException& eventException)
    {
        streamlog_out(ERROR) << "LCIO Event exception: " << eventException.what() << std::endl;
        throw;
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Failed to process event" << std::endl;

        if (STATUS_CODE_SUCCESS != PandoraApi::Reset(*m_pPandora))
        {
            streamlog_out(ERROR) << "Failed to reset Pandora, aborting" << std::endl;
            abort();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::check(LCEvent *pLCEvent)
{
    streamlog_out(MESSAGE) << "PandoraPFANewProcessor - Check" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::end()
{
    delete m_pPandora;
    streamlog_out(MESSAGE) << "PandoraPFANewProcessor - End" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::RegisterUserComponents() const
{
    // Insert user code here ...

    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterEnergyCorrectionFunction(*m_pPandora, "MyHadronicEnergyCorrection",
    //    pandora::HADRONIC, &PandoraPFANewProcessor::MyHadronicEnergyCorrection));

    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterParticleIdFunction(*m_pPandora, "MyParticleId",
    //    &PandoraPFANewProcessor::MyParticleId));

    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*m_pPandora, "MyAlgorithm", new MyAlgorithm::Factory));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::ProcessSteeringFile()
{
    // Insert user code here ...
    registerProcessorParameter("PandoraSettingsXmlFile",
                            "The pandora settings xml file",
                            m_settings.m_pandoraSettingsXmlFile,
                            std::string());

    // Input collections
    registerInputCollections(LCIO::TRACK,
                            "TrackCollections", 
                            "Names of the Track collections used for clustering",
                            m_trackCreator.m_settings.m_trackCollections,
                            StringVector());

    registerInputCollections(LCIO::VERTEX,
                            "V0VertexCollections", 
                            "Name of external V0 Vertex collections",
                            m_trackCreator.m_settings.m_v0VertexCollections,
                            StringVector());

    registerInputCollections(LCIO::VERTEX,
                            "KinkVertexCollections", 
                            "Name of external kink Vertex collections",
                            m_trackCreator.m_settings.m_kinkVertexCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "ECalCaloHitCollections", 
                            "Name of the ECAL calo hit collections",
                            m_caloHitCreator.m_settings.m_eCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "HCalCaloHitCollections", 
                            "Name of the HCAL calo hit collections",
                            m_caloHitCreator.m_settings.m_hCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "LCalCaloHitCollections", 
                            "Name of the LCAL calo hit collections",
                            m_caloHitCreator.m_settings.m_lCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "LHCalCaloHitCollections", 
                            "Name of the LHCAL calo hit collections",
                            m_caloHitCreator.m_settings.m_lHCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "MuonCaloHitCollections", 
                            "Name of the muon calo hit collections",
                            m_caloHitCreator.m_settings.m_muonCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::MCPARTICLE,
                            "MCParticleCollections", 
                            "Name of mc particle collections",
                            m_mcParticleCreator.m_settings.m_mcParticleCollections,
                            StringVector());

    registerInputCollections(LCIO::LCRELATION, 
                            "RelCaloHitCollections",
                            "SimCaloHit to CaloHit Relations Collection Name",
                            m_mcParticleCreator.m_settings.m_lcCaloHitRelationCollections,
                            StringVector());

    registerInputCollections(LCIO::LCRELATION, 
                            "RelTrackCollections",
                            "Track to MCParticle Relations Collection Name",
                            m_mcParticleCreator.m_settings.m_lcTrackRelationCollections,
                            StringVector());

    // Absorber properties
    registerProcessorParameter("AbsorberRadiationLength",
                            "The absorber radation length",
                            m_geometryCreator.m_settings.m_absorberRadiationLength,
                            float(1.));

    registerProcessorParameter("AbsorberInteractionLength",
                            "The absorber interaction length",
                            m_geometryCreator.m_settings.m_absorberInteractionLength,
                            float(1.));

    // Name of PFO collection written by MarlinPandora
    registerOutputCollection( LCIO::CLUSTER,
                              "ClusterCollectionName" , 
                              "Cluster Collection Name "  ,
                              m_pfoCreator.m_settings.m_clusterCollectionName,
                              std::string("PandoraPFANewClusters"));

    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                              "PFOCollectionName" , 
                              "PFO Collection Name "  ,
                              m_pfoCreator.m_settings.m_pfoCollectionName,
                              std::string("PandoraPFANewPFOs"));

    // Calibration constants
    registerProcessorParameter("ECalToMipCalibration",
                            "The calibration from deposited ECal energy to mip",
                            m_caloHitCreator.m_settings.m_eCalToMip,
                            float(1.));

    registerProcessorParameter("HCalToMipCalibration",
                            "The calibration from deposited HCal energy to mip",
                            m_caloHitCreator.m_settings.m_hCalToMip,
                            float(1.));

    registerProcessorParameter("ECalMipThreshold",
                            "Threshold for creating calo hits in the ECal, units mip",
                            m_caloHitCreator.m_settings.m_eCalMipThreshold,
                            float(0.));

    registerProcessorParameter("HCalMipThreshold",
                            "Threshold for creating calo hits in the HCal, units mip",
                            m_caloHitCreator.m_settings.m_hCalMipThreshold,
                            float(0.));

    registerProcessorParameter("ECalToEMGeVCalibration",
                            "The calibration from deposited ECal energy to EM energy",
                            m_caloHitCreator.m_settings.m_eCalToEMGeV,
                            float(1.));

    registerProcessorParameter("HCalToEMGeVCalibration",
                            "The calibration from deposited HCal energy to EM energy",
                            m_caloHitCreator.m_settings.m_hCalToEMGeV,
                            float(1.));

    registerProcessorParameter("ECalToHadGeVCalibration",
                            "The calibration from deposited ECal energy to hadronic energy",
                            m_caloHitCreator.m_settings.m_eCalToHadGeV,
                            float(1.));

    registerProcessorParameter("HCalToHadGeVCalibration",
                            "The calibration from deposited HCal energy to hadronic energy",
                            m_caloHitCreator.m_settings.m_hCalToHadGeV,
                            float(1.));

    registerProcessorParameter("MuonHitEnergy",
                            "The energy for a digital muon calorimeter hit, units GeV",
                            m_caloHitCreator.m_settings.m_muonHitEnergy,
                            float(0.5));

    registerProcessorParameter("MaxHCalHitHadronicEnergy",
                            "The maximum hadronic energy allowed for a single hcal hit",
                            m_caloHitCreator.m_settings.m_maxHCalHitHadronicEnergy,
                            float(1.));

    registerProcessorParameter("NOuterSamplingLayers",
                            "Number of layers from edge for hit to be flagged as an outer layer hit",
                            m_caloHitCreator.m_settings.m_nOuterSamplingLayers,
                            int(3));

    registerProcessorParameter("LayersFromEdgeMaxRearDistance",
                            "Maximum number of layers from candidate outer layer hit to rear of detector",
                            m_caloHitCreator.m_settings.m_layersFromEdgeMaxRearDistance,
                            float(250.f));

    // For calculating track properties
   registerProcessorParameter("MinTrackHits",
                            "Track quality cut: the minimum number of track hits",
                            m_trackCreator.m_settings.m_minTrackHits,
                            int(5));

   registerProcessorParameter("MaxTrackHits",
                            "Track quality cut: the maximum number of track hits",
                            m_trackCreator.m_settings.m_maxTrackHits,
                            int(5000));

    registerProcessorParameter("NumberOfHitsForTrackHelixFits",
                            "The number of hits to be used in helix fits at start/end of tracks",
                            m_trackCreator.m_settings.m_nHitsForHelixFits,
                            int(50));

    registerProcessorParameter("UseEndTrackHelixForECalProjection",
                            "==0 use full track, ==1 use last NumberOfHitsForTrackHelixFit hits",
                            m_trackCreator.m_settings.m_useEndTrackHelixForECalProjection,
                            int(1));

    // Track PFO usage parameters
    registerProcessorParameter("D0TrackCut",
                            "Track d0 cut used to determine whether track can be used to form pfo",
                            m_trackCreator.m_settings.m_d0TrackCut,
                            float(50.));

    registerProcessorParameter("Z0TrackCut",
                            "Track z0 cut used to determine whether track can be used to form pfo",
                            m_trackCreator.m_settings.m_z0TrackCut,
                            float(50.));

    registerProcessorParameter("MaxTpcInnerRDistance",
                            "Track cut on distance from tpc inner r to id whether track can form pfo",
                            m_trackCreator.m_settings.m_maxTpcInnerRDistance,
                            float(50.));

    registerProcessorParameter("UseNonVertexTracks",
                            "Whether can form pfos from tracks that don't start at vertex",
                            m_trackCreator.m_settings.m_usingNonVertexTracks,
                            int(1.));

    registerProcessorParameter("UseUnmatchedNonVertexTracks",
                            "Whether can form pfos from unmatched tracks that don't start at vertex",
                            m_trackCreator.m_settings.m_usingUnmatchedNonVertexTracks,
                            int(0.));

    registerProcessorParameter("UseUnmatchedVertexTracks",
                            "Whether can form pfos from unmatched tracks that start at vertex",
                            m_trackCreator.m_settings.m_usingUnmatchedVertexTracks,
                            int(1.));

    registerProcessorParameter("UnmatchedVertexTrackMaxEnergy",
                            "Maximum energy for unmatched vertex track",
                            m_trackCreator.m_settings.m_unmatchedVertexTrackMaxEnergy,
                            float(5.));

    registerProcessorParameter("D0UnmatchedVertexTrackCut",
                            "d0 cut used to determine whether unmatched vertex track can form pfo",
                            m_trackCreator.m_settings.m_d0UnmatchedVertexTrackCut,
                            float(5.));

    registerProcessorParameter("Z0UnmatchedVertexTrackCut",
                            "z0 cut used to determine whether unmatched vertex track can form pfo",
                            m_trackCreator.m_settings.m_z0UnmatchedVertexTrackCut,
                            float(5.));

    registerProcessorParameter("ZCutForNonVertexTracks",
                            "Non vtx track z cut to determine whether track can be used to form pfo",
                            m_trackCreator.m_settings.m_zCutForNonVertexTracks,
                            float(250.));

    // Track "reaches ecal" parameters
    registerProcessorParameter("ReachesECalNTpcHits",
                            "Minimum number of tpc hits to consider track as reaching ecal",
                            m_trackCreator.m_settings.m_reachesECalNTpcHits,
                            int(10));

    registerProcessorParameter("ReachesECalTpcOuterDistance",
                            "Max distance from track to tpc r max to id whether track reaches ecal",
                            m_trackCreator.m_settings.m_reachesECalTpcOuterDistance,
                            float(-100.));

    registerProcessorParameter("ReachesECalTpcZMaxDistance",
                            "Max distance from track to tpc z max to id whether track reaches ecal",
                            m_trackCreator.m_settings.m_reachesECalTpcZMaxDistance,
                            float(-50.));

    registerProcessorParameter("CurvatureToMomentumFactor",
                            "Constant relating track curvature in b field to momentum",
                            m_trackCreator.m_settings.m_curvatureToMomentumFactor,
                            float(0.3 / 2000.));

    // Track relationship parameters
    registerProcessorParameter("ShouldFormTrackRelationships",
                            "Whether to form pandora track relationships using v0 and kink info",
                            m_trackCreator.m_settings.m_shouldFormTrackRelationships,
                            int(1));

    // Additional geometry parameters
    registerProcessorParameter("ECalEndCapInnerSymmetryOrder",
                            "ECal end cap inner symmetry order (missing from ILD00 gear file)",
                            m_geometryCreator.m_settings.m_eCalEndCapInnerSymmetryOrder,
                            int(4));

    registerProcessorParameter("ECalEndCapInnerPhiCoordinate",
                            "ECal end cap inner phi coordinate (missing from ILD00 gear file)",
                            m_geometryCreator.m_settings.m_eCalEndCapInnerPhiCoordinate,
                            float(0.));

    registerProcessorParameter("HCalEndCapInnerSymmetryOrder",
                            "HCal end cap inner symmetry order (missing from ILD00 gear file)",
                            m_geometryCreator.m_settings.m_hCalEndCapInnerSymmetryOrder,
                            int(8));

    registerProcessorParameter("HCalEndCapInnerPhiCoordinate",
                            "HCal end cap inner phi coordinate (missing from ILD00 gear file)",
                            m_geometryCreator.m_settings.m_hCalEndCapInnerPhiCoordinate,
                            float(0.));

    // Number of events to skip
    registerProcessorParameter("NEventsToSkip",
                            "Number of events to skip at start of reconstruction job",
                            m_settings.m_nEventsToSkip,
                            int(0));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::FinaliseSteeringParameters()
{
    m_caloHitCreator.m_settings.m_absorberRadiationLength = m_geometryCreator.m_settings.m_absorberRadiationLength;
    m_caloHitCreator.m_settings.m_absorberInteractionLength = m_geometryCreator.m_settings.m_absorberInteractionLength;
    m_caloHitCreator.m_settings.m_hCalEndCapInnerSymmetryOrder = m_geometryCreator.m_settings.m_hCalEndCapInnerSymmetryOrder;
    m_caloHitCreator.m_settings.m_hCalEndCapInnerPhiCoordinate = m_geometryCreator.m_settings.m_hCalEndCapInnerPhiCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::Reset()
{
    m_caloHitCreator.Reset();
    m_trackCreator.Reset();
}
