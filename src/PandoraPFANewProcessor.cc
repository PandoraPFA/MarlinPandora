/**
 *  @file   PandoraPFANew/src/PandoraPFANewProcessor.cc
 * 
 *  @brief  Implementation of the pandora pfa new processor class.
 * 
 *  $Log: $
 */

#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "EVENT/MCParticle.h"

#include "UTIL/CellIDDecoder.h"

#include <cmath>

// User algorithm includes here


#include "Api/PandoraApi.h"

#include "PandoraPFANewProcessor.h"

PandoraPFANewProcessor pandoraPFANewProcessor;

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
        std::cout << "PandoraPFANewProcessor - Init" << std::endl;

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateGeometry());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RegisterUserAlgorithmFactories());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(m_pandora, m_settings.m_pandoraSettingsXmlFile));
    }
    catch (StatusCodeException &statusCodeException)
    {
        std::cout << "Failed to initialize pandora pfa new processor: " << statusCodeException.ToString() << std::endl;
        throw;
    }
    catch (...)
    {
        std::cout << "Failed to initialize pandora pfa new processor, unrecognized exception" << std::endl;
        throw;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::processRunHeader(LCRunHeader *pLCRunHeader)
{
    m_detectorName = pLCRunHeader->getDetectorName();
    std::cout << "Detector Name " << m_detectorName << ", Run " << ++m_nRun <<  std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::processEvent(LCEvent *pLCEvent)
{
    try
    {
        std::cout << "Run " << m_nRun << ", Event " << ++m_nEvent << std::endl;

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateMCParticles(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateTracks(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateCaloHits(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(m_pandora));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ProcessParticleFlowObjects(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(m_pandora));
    }
    catch (StatusCodeException &statusCodeException)
    {
        std::cout << "Failed to process event: " << statusCodeException.ToString() << std::endl;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(m_pandora));
    }
    catch (...)
    {
        std::cout << "Failed to process event, unrecognized exception" << std::endl;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(m_pandora));        
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::check(LCEvent *pLCEvent)
{
    std::cout << "PandoraPFANewProcessor - Check" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::end()
{
    std::cout << "PandoraPFANewProcessor - End" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateGeometry()
{
    // Insert user code here ...
    PandoraApi::Geometry::Parameters geometryParameters;

    geometryParameters.m_mainTrackerInnerRadius = 1;
    geometryParameters.m_mainTrackerOuterRadius = 2;
    geometryParameters.m_mainTrackerZExtent     = 3;

    geometryParameters.m_nRadLengthsInZGap      = 4;
    geometryParameters.m_nIntLengthsInZGap      = 5;
    geometryParameters.m_nRadLengthsInRadialGap = 6;
    geometryParameters.m_nIntLengthsInRadialGap = 7;

    PandoraApi::Geometry::Parameters::LayerParameters layerParameters;
    layerParameters.m_distanceFromIp        = 1;
    layerParameters.m_nRadiationLengths     = 2;
    layerParameters.m_nInteractionLengths   = 3;

    geometryParameters.m_eCalBarrelParameters.m_innerDistanceFromIp = 1;
    geometryParameters.m_eCalBarrelParameters.m_innerSymmetry       = 2;
    geometryParameters.m_eCalBarrelParameters.m_innerAngle          = 3;
    geometryParameters.m_eCalBarrelParameters.m_outerDistanceFromIp = 4;
    geometryParameters.m_eCalBarrelParameters.m_outerSymmetry       = 5;
    geometryParameters.m_eCalBarrelParameters.m_outerAngle          = 6;
    geometryParameters.m_eCalBarrelParameters.m_nLayers             = 1;
    geometryParameters.m_eCalBarrelParameters.m_layerParametersList.push_back(layerParameters);

    geometryParameters.m_hCalBarrelParameters.m_innerDistanceFromIp = 1;
    geometryParameters.m_hCalBarrelParameters.m_innerSymmetry       = 2;
    geometryParameters.m_hCalBarrelParameters.m_innerAngle          = 3;
    geometryParameters.m_hCalBarrelParameters.m_outerDistanceFromIp = 4;
    geometryParameters.m_hCalBarrelParameters.m_outerSymmetry       = 5;
    geometryParameters.m_hCalBarrelParameters.m_outerAngle          = 6;
    geometryParameters.m_hCalBarrelParameters.m_nLayers             = 1;
    geometryParameters.m_hCalBarrelParameters.m_layerParametersList.push_back(layerParameters);

    geometryParameters.m_eCalEndCapParameters.m_innerDistanceFromIp = 1;
    geometryParameters.m_eCalEndCapParameters.m_innerSymmetry       = 2;
    geometryParameters.m_eCalEndCapParameters.m_innerAngle          = 3;
    geometryParameters.m_eCalEndCapParameters.m_outerDistanceFromIp = 4;
    geometryParameters.m_eCalEndCapParameters.m_outerSymmetry       = 5;
    geometryParameters.m_eCalEndCapParameters.m_outerAngle          = 6;
    geometryParameters.m_eCalEndCapParameters.m_nLayers             = 1;
    geometryParameters.m_eCalEndCapParameters.m_layerParametersList.push_back(layerParameters);

    geometryParameters.m_hCalEndCapParameters.m_innerDistanceFromIp = 1;
    geometryParameters.m_hCalEndCapParameters.m_innerSymmetry       = 2;
    geometryParameters.m_hCalEndCapParameters.m_innerAngle          = 3;
    geometryParameters.m_hCalEndCapParameters.m_outerDistanceFromIp = 4;
    geometryParameters.m_hCalEndCapParameters.m_outerSymmetry       = 5;
    geometryParameters.m_hCalEndCapParameters.m_outerAngle          = 6;
    geometryParameters.m_hCalEndCapParameters.m_nLayers             = 1;
    geometryParameters.m_hCalEndCapParameters.m_layerParametersList.push_back(layerParameters);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::Create(m_pandora, geometryParameters));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::RegisterUserAlgorithmFactories()
{
    // Insert user code here ...

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateMCParticles(const LCEvent *const pLCEvent)
{
    // Insert user code here ...
    for (StringVector::const_iterator iter = m_settings.m_mcParticleCollections.begin(), iterEnd = m_settings.m_mcParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pMCParticleCollection = pLCEvent->getCollection(*iter);

            const int numberMCParticles = pMCParticleCollection->getNumberOfElements();

            for(int i = 0; i < numberMCParticles; ++i)
            {
                MCParticle* pMcParticle = dynamic_cast<MCParticle*>(pMCParticleCollection->getElementAt(i));

                double innerRadius = 0.;
                double outerRadius = 0.;
                double momentum    = 0.;

                for(int i = 0; i < 3; ++i)
                {
                    innerRadius += pow(pMcParticle->getVertex()[i], 2);
                    outerRadius += pow(pMcParticle->getEndpoint()[i], 2);
                    momentum    += pow(pMcParticle->getMomentum()[i], 2);
                }

                innerRadius = std::sqrt(innerRadius);
                outerRadius = std::sqrt(outerRadius);
                momentum    = std::sqrt(momentum);
     
                PandoraApi::MCParticle::Parameters mcParticleParameters;
                mcParticleParameters.m_energy = pMcParticle->getEnergy();
                mcParticleParameters.m_momentum = momentum;
                mcParticleParameters.m_innerRadius = innerRadius;
                mcParticleParameters.m_outerRadius = outerRadius;
                mcParticleParameters.m_pParentAddress = pMcParticle;

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters));

                // Create parent-daughter relationships
                for(MCParticleVec::const_iterator itDaughter = pMcParticle->getDaughters().begin(),
                    itDaughterEnd = pMcParticle->getDaughters().end(); itDaughter != itDaughterEnd; ++itDaughter)
                {
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(m_pandora, pMcParticle,
                        *itDaughter));
                }
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            std::cout << "Failed to extract MCParticles: " << statusCodeException.ToString() << std::endl;
        }
        catch (...)
        {
            std::cout << "Failed to extract MCParticles, unrecognised exception" << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateTracks(const LCEvent *const pLCEvent)
{
    // Insert user code here ...
    for (StringVector::const_iterator iter = m_settings.m_trackCollections.begin(), 
        iterEnd = m_settings.m_trackCollections.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pTrackCollection = pLCEvent->getCollection(*iter);
            
            for (int i = 0; i < pTrackCollection->getNumberOfElements(); ++i)
            {
                Track* pTrack = dynamic_cast<Track*>(pTrackCollection->getElementAt(i));

                PandoraApi::Track::Parameters trackParameters;
                trackParameters.m_d0 = 1;
                trackParameters.m_z0 = 2;
                trackParameters.m_momentumAtDca = pandora::CartesianVector(1, 2, 3);

                trackParameters.m_trackStateAtStart = pandora::TrackState(1, 2, 3, 4, 5, 6);
                trackParameters.m_trackStateAtEnd = pandora::TrackState(1, 2, 3, 4, 5, 6);
                trackParameters.m_trackStateAtECal = pandora::TrackState(1, 2, 3, 4, 5, 6);

                trackParameters.m_calorimeterProjections.push_back(pandora::TrackState(1, 2, 3, 4, 5, 6));
                trackParameters.m_calorimeterProjections.push_back(pandora::TrackState(7, 8, 9, 10, 11, 12));

                trackParameters.m_reachesECal = true;
                trackParameters.m_pParentAddress = pTrack;

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(m_pandora, trackParameters));
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            std::cout << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
        }
        catch (...)
        {
            std::cout << "Failed to extract a track, unrecognised exception" << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateCaloHits(const LCEvent *const pLCEvent)
{
    // Insert user code here ...
    for (StringVector::const_iterator iter = m_settings.m_caloHitCollections.begin(), 
        iterEnd = m_settings.m_caloHitCollections.end(); iter != iterEnd; ++iter)
    {
        try
        {
            CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");

            const LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);

            for (int i = 0; i < pCaloHitCollection->getNumberOfElements(); ++i)
            {
                CalorimeterHit* pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                PandoraApi::CaloHit::Parameters caloHitParameters;
                caloHitParameters.m_positionVector = pandora::CartesianVector(1, 2, 3);
                caloHitParameters.m_normalVector = pandora::CartesianVector(4, 5, 6);

                caloHitParameters.m_cellSizeU = 1;
                caloHitParameters.m_cellSizeV = 2;
                caloHitParameters.m_cellSizeZ = 3;

                caloHitParameters.m_nRadiationLengths = 4;
                caloHitParameters.m_nInteractionLengths = 5;

                caloHitParameters.m_energy = pCaloHit->getEnergy();
                caloHitParameters.m_time = 7;

                caloHitParameters.m_isDigital = false;
                caloHitParameters.m_hitType = pandora::ECAL;
                caloHitParameters.m_detectorRegion = pandora::BARREL;

                caloHitParameters.m_layer = cellIdDecoder(pCaloHit)["K-1"];
                caloHitParameters.m_pParentAddress = pCaloHit;

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(m_pandora, caloHitParameters));
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            std::cout << "Failed to extract a calo hit: " << statusCodeException.ToString() << std::endl;
        }
        catch (...)
        {
            std::cout << "Failed to extract a calo hit, unrecognised exception" << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::ProcessParticleFlowObjects(const LCEvent *const pLCEvent)
{
    // Insert user code here ...
    PandoraApi::ParticleFlowObjectList particleFlowObjectList;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=,
        PandoraApi::GetParticleFlowObjects(m_pandora, particleFlowObjectList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode PandoraPFANewProcessor::CreateMCTrees(const LCEvent *const pLCEvent)
{
    // Insert user code here ...

    return STATUS_CODE_SUCCESS;
}

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
                            m_settings.m_trackCollections,
                            StringVector());

    registerInputCollections(LCIO::VERTEX,
                            "V0VertexCollections", 
                            "Name of external V0 Vertex collections",
                            m_settings.m_v0VertexCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "CaloHitcollections", 
                            "Name of the HCAL collection used to form clusters",
                            m_settings.m_caloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::MCPARTICLE,
                            "MCParticleCollections", 
                            "Name of mc particle collections",
                            m_settings.m_mcParticleCollections,
                            StringVector());
}
