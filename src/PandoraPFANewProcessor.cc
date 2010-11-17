/**
 *  @file   MarlinPandora/src/PandoraPFANewProcessor.cc
 * 
 *  @brief  Implementation of the pandora pfa new processor class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Exceptions.h"

#include "gear/BField.h"

#include "Api/PandoraApi.h"

#include "Utilities/HighGranularityPseudoLayerCalculator.h"

#include "ExternalClusteringAlgorithm.h"
#include "InteractionLengthCalculator.h"
#include "PandoraPFANewProcessor.h"
#include "SimpleBFieldCalculator.h"

#include <cstdlib>

PandoraPFANewProcessor pandoraPFANewProcessor;

pandora::Pandora *PandoraPFANewProcessor::m_pPandora = NULL;
EVENT::LCEvent *PandoraPFANewProcessor::m_pLcioEvent = NULL;

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraPFANewProcessor::PandoraPFANewProcessor() :
    Processor("PandoraPFANewProcessor"),
    m_pGeometryCreator(NULL),
    m_pCaloHitCreator(NULL),
    m_pTrackCreator(NULL),
    m_pMCParticleCreator(NULL),
    m_pPfoCreator(NULL),
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
        m_pGeometryCreator = new GeometryCreator(m_geometryCreatorSettings);
        m_pCaloHitCreator = new CaloHitCreator(m_caloHitCreatorSettings);
        m_pTrackCreator = new TrackCreator(m_trackCreatorSettings);
        m_pMCParticleCreator = new MCParticleCreator(m_mcParticleCreatorSettings);
        m_pPfoCreator = new PfoCreator(m_pfoCreatorSettings);

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->RegisterUserComponents());
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pGeometryCreator->CreateGeometry());
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPandora, m_settings.m_pandoraSettingsXmlFile));
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "Failed to initialize marlin pandora: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(ERROR) << "Failed to initialize marlin pandora: gear exception " << exception.what() << std::endl;
        throw exception;
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Failed to initialize marlin pandora: unrecognized exception" << std::endl;
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
    m_pLcioEvent = pLCEvent;

    if (eventCounter < m_settings.m_nEventsToSkip)
    {
        ++eventCounter;
        throw marlin::SkipEventException(this);
    }

    try
    {
        streamlog_out(MESSAGE) << "PandoraPFANewProcessor, Run " << m_nRun << ", Event " << ++m_nEvent << std::endl;

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateMCParticles(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTrackAssociations(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTracks(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateTrackToMCParticleRelationships(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pCaloHitCreator->CreateCaloHits(pLCEvent));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateCaloHitToMCParticleRelationships(pLCEvent));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPandora));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pPfoCreator->CreateParticleFlowObjects(pLCEvent));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPandora));
        this->Reset();
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "Marlin pandora failed to process event: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(ERROR) << "Marlin pandora failed to process event: gear exception " << exception.what() << std::endl;
        throw exception;
    }
    catch (EVENT::Exception &exception)
    {
        streamlog_out(ERROR) << "Marlin pandora failed to process event: lcio exception " << exception.what() << std::endl;
        throw exception;
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Marlin pandora failed to process event: unrecognized exception" << std::endl;
        throw;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::check(LCEvent */*pLCEvent*/)
{
    streamlog_out(MESSAGE) << "PandoraPFANewProcessor - Check" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::end()
{
    delete m_pPandora;
    delete m_pGeometryCreator;
    delete m_pCaloHitCreator;
    delete m_pTrackCreator;
    delete m_pMCParticleCreator;
    delete m_pPfoCreator;

    streamlog_out(MESSAGE) << "PandoraPFANewProcessor - End" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode PandoraPFANewProcessor::RegisterUserComponents() const
{
    //PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterEnergyCorrectionFunction(*m_pPandora,
    //    "MyHadronicEnergyCorrection", pandora::HADRONIC, &PandoraPFANewProcessor::MyHadronicEnergyCorrection));

    //PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterParticleIdFunction(*m_pPandora, "MyParticleId",
    //    &PandoraPFANewProcessor::MyParticleId));

    //PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*m_pPandora, "MyAlgorithm",
    //    new MyAlgorithm::Factory));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*m_pPandora, "ExternalClustering",
        new ExternalClusteringAlgorithm::Factory));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerCalculator(*m_pPandora,
        new pandora::HighGranularityPseudoLayerCalculator()));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetBFieldCalculator(*m_pPandora,
        new SimpleBFieldCalculator()));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::ProcessSteeringFile()
{
    registerProcessorParameter("PandoraSettingsXmlFile",
                            "The pandora settings xml file",
                            m_settings.m_pandoraSettingsXmlFile,
                            std::string());

    // Input collections
    registerInputCollections(LCIO::TRACK,
                            "TrackCollections", 
                            "Names of the Track collections used for clustering",
                            m_trackCreatorSettings.m_trackCollections,
                            StringVector());

    registerInputCollections(LCIO::VERTEX,
                            "KinkVertexCollections", 
                            "Name of external kink Vertex collections",
                            m_trackCreatorSettings.m_kinkVertexCollections,
                            StringVector());

    registerInputCollections(LCIO::VERTEX,
                            "ProngVertexCollections", 
                            "Name of external prong Vertex collections",
                            m_trackCreatorSettings.m_prongVertexCollections,
                            StringVector());

    registerInputCollections(LCIO::VERTEX,
                            "SplitVertexCollections", 
                            "Name of external split Vertex collections",
                            m_trackCreatorSettings.m_splitVertexCollections,
                            StringVector());

    registerInputCollections(LCIO::VERTEX,
                            "V0VertexCollections", 
                            "Name of external V0 Vertex collections",
                            m_trackCreatorSettings.m_v0VertexCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "ECalCaloHitCollections", 
                            "Name of the ECAL calo hit collections",
                            m_caloHitCreatorSettings.m_eCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "HCalCaloHitCollections", 
                            "Name of the HCAL calo hit collections",
                            m_caloHitCreatorSettings.m_hCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "LCalCaloHitCollections", 
                            "Name of the LCAL calo hit collections",
                            m_caloHitCreatorSettings.m_lCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "LHCalCaloHitCollections", 
                            "Name of the LHCAL calo hit collections",
                            m_caloHitCreatorSettings.m_lHCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "MuonCaloHitCollections", 
                            "Name of the muon calo hit collections",
                            m_caloHitCreatorSettings.m_muonCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::MCPARTICLE,
                            "MCParticleCollections", 
                            "Name of mc particle collections",
                            m_mcParticleCreatorSettings.m_mcParticleCollections,
                            StringVector());

    registerInputCollections(LCIO::LCRELATION, 
                            "RelCaloHitCollections",
                            "SimCaloHit to CaloHit Relations Collection Name",
                            m_mcParticleCreatorSettings.m_lcCaloHitRelationCollections,
                            StringVector());

    registerInputCollections(LCIO::LCRELATION, 
                            "RelTrackCollections",
                            "Track to MCParticle Relations Collection Name",
                            m_mcParticleCreatorSettings.m_lcTrackRelationCollections,
                            StringVector());

    // Absorber properties
    registerProcessorParameter("AbsorberRadiationLength",
                            "The absorber radation length",
                            m_geometryCreatorSettings.m_absorberRadiationLength,
                            float(1.));

    registerProcessorParameter("AbsorberInteractionLength",
                            "The absorber interaction length",
                            m_geometryCreatorSettings.m_absorberInteractionLength,
                            float(1.));

    // Name of PFO collection written by MarlinPandora
    registerOutputCollection( LCIO::CLUSTER,
                              "ClusterCollectionName" , 
                              "Cluster Collection Name "  ,
                              m_pfoCreatorSettings.m_clusterCollectionName,
                              std::string("PandoraPFANewClusters"));

    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                              "PFOCollectionName" , 
                              "PFO Collection Name "  ,
                              m_pfoCreatorSettings.m_pfoCollectionName,
                              std::string("PandoraPFANewPFOs"));

    // Calibration constants
    registerProcessorParameter("ECalToMipCalibration",
                            "The calibration from deposited ECal energy to mip",
                            m_caloHitCreatorSettings.m_eCalToMip,
                            float(1.));

    registerProcessorParameter("HCalToMipCalibration",
                            "The calibration from deposited HCal energy to mip",
                            m_caloHitCreatorSettings.m_hCalToMip,
                            float(1.));

    registerProcessorParameter("ECalMipThreshold",
                            "Threshold for creating calo hits in the ECal, units mip",
                            m_caloHitCreatorSettings.m_eCalMipThreshold,
                            float(0.));

    registerProcessorParameter("MuonToMipCalibration",
                            "The calibration from deposited Muon energy to mip",
                            m_caloHitCreatorSettings.m_muonToMip,
                            float(10.));

    registerProcessorParameter("HCalMipThreshold",
                            "Threshold for creating calo hits in the HCal, units mip",
                            m_caloHitCreatorSettings.m_hCalMipThreshold,
                            float(0.));

    registerProcessorParameter("ECalToEMGeVCalibration",
                            "The calibration from deposited ECal energy to EM energy",
                            m_caloHitCreatorSettings.m_eCalToEMGeV,
                            float(1.));

    registerProcessorParameter("HCalToEMGeVCalibration",
                            "The calibration from deposited HCal energy to EM energy",
                            m_caloHitCreatorSettings.m_hCalToEMGeV,
                            float(1.));

    registerProcessorParameter("ECalToHadGeVCalibration",
                            "The calibration from deposited ECal energy to hadronic energy",
                            m_caloHitCreatorSettings.m_eCalToHadGeV,
                            float(1.));

    registerProcessorParameter("HCalToHadGeVCalibration",
                            "The calibration from deposited HCal energy to hadronic energy",
                            m_caloHitCreatorSettings.m_hCalToHadGeV,
                            float(1.));

    registerProcessorParameter("DigitalMuonHits",
                            "Treat muon hits as digital",
                            m_caloHitCreatorSettings.m_muonDigitalHits,
                            int(1));

    registerProcessorParameter("MuonHitEnergy",
                            "The energy for a digital muon calorimeter hit, units GeV",
                            m_caloHitCreatorSettings.m_muonHitEnergy,
                            float(0.5));

    registerProcessorParameter("MaxHCalHitHadronicEnergy",
                            "The maximum hadronic energy allowed for a single hcal hit",
                            m_caloHitCreatorSettings.m_maxHCalHitHadronicEnergy,
                            float(1.));

    registerProcessorParameter("NOuterSamplingLayers",
                            "Number of layers from edge for hit to be flagged as an outer layer hit",
                            m_caloHitCreatorSettings.m_nOuterSamplingLayers,
                            int(3));

    registerProcessorParameter("LayersFromEdgeMaxRearDistance",
                            "Maximum number of layers from candidate outer layer hit to rear of detector",
                            m_caloHitCreatorSettings.m_layersFromEdgeMaxRearDistance,
                            float(250.f));

    // B-field parameters
    registerProcessorParameter("MuonBarrelBField",
                            "The bfield in the muon barrel, units Tesla",
                            SimpleBFieldCalculator::m_muonBarrelBField,
                            float(-1.5f));

    registerProcessorParameter("MuonEndCapBField",
                            "The bfield in the muon endcap, units Tesla",
                            SimpleBFieldCalculator::m_muonEndCapBField,
                            float(0.01f));

    // Average interaction length parameters
    registerProcessorParameter("AverageInteractionLengthTracker",
                            "average number interaction length per mm in the tracker",
                            InteractionLengthCalculator::Settings::m_avgIntLengthTracker,
                            float(0.f));

    registerProcessorParameter("AverageInteractionLengthCoil",
                            "average number interaction length per mm in the coil",
                            InteractionLengthCalculator::Settings::m_avgIntLengthCoil,
                            float(0.0025189));

    registerProcessorParameter("AverageInteractionLengthECalBarrel",
                            "average number interaction length per mm in the ECal Barrel",
                            InteractionLengthCalculator::Settings::m_avgIntLengthECalBarrel,
                            float(0.0040269f));

    registerProcessorParameter("AverageInteractionLengthHCalBarrel",
                            "average number interaction length per mm in the HCal Barrel",
                            InteractionLengthCalculator::Settings::m_avgIntLengthHCalBarrel,
                            float(0.0033041916f));

    registerProcessorParameter("AverageInteractionLengthECalEndCap",
                            "average number interaction length per mm in the ECal EndCap",
                            InteractionLengthCalculator::Settings::m_avgIntLengthECalEndCap,
                            float(0.0040269f));

    registerProcessorParameter("AverageInteractionLengthHCalEndCap",
                            "average number interaction length per mm in the HCal EndCap",
                            InteractionLengthCalculator::Settings::m_avgIntLengthHCalEndCap,
                            float(0.0033041916));

    registerProcessorParameter("AverageInteractionLengthMuonBarrel",
                            "average number interaction length per mm in the Muon Barrel",
                            InteractionLengthCalculator::Settings::m_avgIntLengthMuonBarrel,
                            float(0.0040269f));

    registerProcessorParameter("AverageInteractionLengthMuonEndCap",
                            "average number interaction length per mm in the Muon EndCap",
                            InteractionLengthCalculator::Settings::m_avgIntLengthMuonEndCap,
                            float(0.0040269f));

    // Track relationship parameters
    registerProcessorParameter("ShouldFormTrackRelationships",
                            "Whether to form pandora track relationships using v0 and kink info",
                            m_trackCreatorSettings.m_shouldFormTrackRelationships,
                            int(1));

    // Initial track hit specifications
   registerProcessorParameter("MinTrackHits",
                            "Track quality cut: the minimum number of track hits",
                            m_trackCreatorSettings.m_minTrackHits,
                            int(5));

   registerProcessorParameter("MinFtdTrackHits",
                            "Track quality cut: the minimum number of ftd track hits for ftd only tracks",
                            m_trackCreatorSettings.m_minTrackHits,
                            int(4));

   registerProcessorParameter("MaxTrackHits",
                            "Track quality cut: the maximum number of track hits",
                            m_trackCreatorSettings.m_maxTrackHits,
                            int(5000));

    // Track PFO usage parameters
    registerProcessorParameter("D0TrackCut",
                            "Track d0 cut used to determine whether track can be used to form pfo",
                            m_trackCreatorSettings.m_d0TrackCut,
                            float(50.));

    registerProcessorParameter("Z0TrackCut",
                            "Track z0 cut used to determine whether track can be used to form pfo",
                            m_trackCreatorSettings.m_z0TrackCut,
                            float(50.));

    registerProcessorParameter("UseNonVertexTracks",
                            "Whether can form pfos from tracks that don't start at vertex",
                            m_trackCreatorSettings.m_usingNonVertexTracks,
                            int(1.));

    registerProcessorParameter("UseUnmatchedNonVertexTracks",
                            "Whether can form pfos from unmatched tracks that don't start at vertex",
                            m_trackCreatorSettings.m_usingUnmatchedNonVertexTracks,
                            int(0.));

    registerProcessorParameter("UseUnmatchedVertexTracks",
                            "Whether can form pfos from unmatched tracks that start at vertex",
                            m_trackCreatorSettings.m_usingUnmatchedVertexTracks,
                            int(1.));

    registerProcessorParameter("UnmatchedVertexTrackMaxEnergy",
                            "Maximum energy for unmatched vertex track",
                            m_trackCreatorSettings.m_unmatchedVertexTrackMaxEnergy,
                            float(5.));

    registerProcessorParameter("D0UnmatchedVertexTrackCut",
                            "d0 cut used to determine whether unmatched vertex track can form pfo",
                            m_trackCreatorSettings.m_d0UnmatchedVertexTrackCut,
                            float(5.));

    registerProcessorParameter("Z0UnmatchedVertexTrackCut",
                            "z0 cut used to determine whether unmatched vertex track can form pfo",
                            m_trackCreatorSettings.m_z0UnmatchedVertexTrackCut,
                            float(5.));

    registerProcessorParameter("ZCutForNonVertexTracks",
                            "Non vtx track z cut to determine whether track can be used to form pfo",
                            m_trackCreatorSettings.m_zCutForNonVertexTracks,
                            float(250.));

    // Track "reaches ecal" parameters
    registerProcessorParameter("ReachesECalNTpcHits",
                            "Minimum number of tpc hits to consider track as reaching ecal",
                            m_trackCreatorSettings.m_reachesECalNTpcHits,
                            int(11));

    registerProcessorParameter("ReachesECalNFtdHits",
                            "Minimum number of ftd hits to consider track as reaching ecal",
                            m_trackCreatorSettings.m_reachesECalNFtdHits,
                            int(4));

    registerProcessorParameter("ReachesECalTpcOuterDistance",
                            "Max distance from track to tpc r max to id whether track reaches ecal",
                            m_trackCreatorSettings.m_reachesECalTpcOuterDistance,
                            float(-100.));

    registerProcessorParameter("ReachesECalTpcZMaxDistance",
                            "Max distance from track to tpc z max to id whether track reaches ecal",
                            m_trackCreatorSettings.m_reachesECalTpcZMaxDistance,
                            float(-50.));

    registerProcessorParameter("ReachesECalFtdZMaxDistance",
                            "Max distance from track hit to ftd z position to identify ftd hits",
                            m_trackCreatorSettings.m_reachesECalFtdZMaxDistance,
                            float(1.));

    registerProcessorParameter("CurvatureToMomentumFactor",
                            "Constant relating track curvature in b field to momentum",
                            m_trackCreatorSettings.m_curvatureToMomentumFactor,
                            float(0.3 / 2000.));

    registerProcessorParameter("MinTrackECalDistanceFromIp",
                            "Sanity check on separation between ip and track projected ecal position",
                            m_trackCreatorSettings.m_minTrackECalDistanceFromIp,
                            float(100.));

    // Final track quality parameters
    registerProcessorParameter("MaxTrackSigmaPOverP",
                            "Cut on fractional track momentum error",
                            m_trackCreatorSettings.m_maxTrackSigmaPOverP,
                            float(0.15));

    registerProcessorParameter("MinMomentumForTrackHitChecks",
                            "Min track momentum required to perform final quality checks on number of hits",
                            m_trackCreatorSettings.m_minMomentumForTrackHitChecks,
                            float(1.));

    registerProcessorParameter("TpcMembraneMaxZ",
                            "Tpc membrane max z coordinate",
                            m_trackCreatorSettings.m_tpcMembraneMaxZ,
                            float(10.));

    registerProcessorParameter("MinTpcHitFractionOfExpected",
                            "Cut on fractional of expected number of TPC hits",
                            m_trackCreatorSettings.m_minTpcHitFractionOfExpected,
                            float(0.20));

    registerProcessorParameter("MinFtdHitsForTpcHitFraction",
                            "Cut on minimum number of FTD hits of TPC hit fraction to be applied",
                            m_trackCreatorSettings.m_minFtdHitsForTpcHitFraction,
                            int(2));

    registerProcessorParameter("MaxTpcInnerRDistance",
                            "Track cut on distance from tpc inner r to id whether track can form pfo",
                            m_trackCreatorSettings.m_maxTpcInnerRDistance,
                            float(50.));

    // Additional geometry parameters
    registerProcessorParameter("ECalEndCapInnerSymmetryOrder",
                            "ECal end cap inner symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_eCalEndCapInnerSymmetryOrder,
                            int(4));

    registerProcessorParameter("ECalEndCapInnerPhiCoordinate",
                            "ECal end cap inner phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_eCalEndCapInnerPhiCoordinate,
                            float(0.));

    registerProcessorParameter("ECalEndCapOuterSymmetryOrder",
                            "ECal end cap outer symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_eCalEndCapOuterSymmetryOrder,
                            int(8));

    registerProcessorParameter("ECalEndCapOuterPhiCoordinate",
                            "ECal end cap outer phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_eCalEndCapOuterPhiCoordinate,
                            float(0.));

    registerProcessorParameter("HCalEndCapInnerSymmetryOrder",
                            "HCal end cap inner symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalEndCapInnerSymmetryOrder,
                            int(4));

    registerProcessorParameter("HCalEndCapInnerPhiCoordinate",
                            "HCal end cap inner phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalEndCapInnerPhiCoordinate,
                            float(0.));

    registerProcessorParameter("HCalEndCapOuterSymmetryOrder",
                            "HCal end cap outer symmetry order (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalEndCapOuterSymmetryOrder,
                            int(16));

    registerProcessorParameter("HCalEndCapOuterPhiCoordinate",
                            "HCal end cap outer phi coordinate (missing from ILD gear files)",
                            m_geometryCreatorSettings.m_hCalEndCapOuterPhiCoordinate,
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
    // ATTN: This function seems to be necessary for operations that cannot easily be performed at construction of the processor,
    // when the steering file is parsed e.g. the call to GEAR to get the inner bfield
    m_caloHitCreatorSettings.m_absorberRadiationLength = m_geometryCreatorSettings.m_absorberRadiationLength;
    m_caloHitCreatorSettings.m_absorberInteractionLength = m_geometryCreatorSettings.m_absorberInteractionLength;
    m_caloHitCreatorSettings.m_hCalEndCapInnerSymmetryOrder = m_geometryCreatorSettings.m_hCalEndCapInnerSymmetryOrder;
    m_caloHitCreatorSettings.m_hCalEndCapInnerPhiCoordinate = m_geometryCreatorSettings.m_hCalEndCapInnerPhiCoordinate;

    m_trackCreatorSettings.m_prongSplitVertexCollections = m_trackCreatorSettings.m_prongVertexCollections;
    m_trackCreatorSettings.m_prongSplitVertexCollections.insert(m_trackCreatorSettings.m_prongSplitVertexCollections.end(),
        m_trackCreatorSettings.m_splitVertexCollections.begin(), m_trackCreatorSettings.m_splitVertexCollections.end());

    SimpleBFieldCalculator::m_innerBField = marlin::Global::GEAR->getBField().at(gear::Vector3D(0., 0., 0.)).z();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::Reset()
{
    m_pCaloHitCreator->Reset();
    m_pTrackCreator->Reset();
}