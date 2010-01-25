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
#include "EVENT/Vertex.h"
#include "EVENT/MCParticle.h"
#include "EVENT/SimCalorimeterHit.h"

#include "IMPL/ClusterImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "UTIL/CellIDDecoder.h"
#include "UTIL/LCRelationNavigator.h"

#include "marlin/Global.h"

#include "gear/BField.h"
#include "gear/GEAR.h"
#include "gear/GearParameters.h"
#include "gear/CalorimeterParameters.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/LayerLayout.h"

#include "CalorimeterHitType.h"
#include "ClusterShapes.h"
#include "HelixClass.h"

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
        streamlog_out(MESSAGE) << "PandoraPFANewProcessor - Init" << std::endl;
        m_pPandora = new pandora::Pandora();

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateGeometry());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RegisterUserAlgorithmFactories());
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
    try
    {
        streamlog_out(MESSAGE) << "Run " << m_nRun << ", Event " << ++m_nEvent << std::endl;

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateMCParticles(pLCEvent));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateTrackAssociations(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateTracks(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateTrackToMCParticleRelationships(pLCEvent));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateCaloHits(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateCaloHitToMCParticleRelationships(pLCEvent));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPandora));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ProcessParticleFlowObjects(pLCEvent));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPandora));
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

StatusCode PandoraPFANewProcessor::CreateGeometry() const
{
    try
    {
        // Insert user code here ...
        PandoraApi::Geometry::Parameters geometryParameters;

        const gear::TPCParameters &tpcParameters    = marlin::Global::GEAR->getTPCParameters();
        const gear::PadRowLayout2D &tpcPadLayout    = tpcParameters.getPadLayout();
        geometryParameters.m_mainTrackerInnerRadius = tpcPadLayout.getPlaneExtent()[0];
        geometryParameters.m_mainTrackerOuterRadius = tpcPadLayout.getPlaneExtent()[1];
        geometryParameters.m_mainTrackerZExtent     = tpcParameters.getMaxDriftLength();

        const gear::GearParameters &coilParameters  = marlin::Global::GEAR->getGearParameters("CoilParameters");
        geometryParameters.m_coilInnerRadius        = coilParameters.getDoubleVal("Coil_cryostat_inner_radius");
        geometryParameters.m_coilOuterRadius        = coilParameters.getDoubleVal("Coil_cryostat_outer_radius");
        geometryParameters.m_coilZExtent            = coilParameters.getDoubleVal("Coil_cryostat_half_z");
        geometryParameters.m_bField                 = marlin::Global::GEAR->getBField().at(gear::Vector3D(0., 0., 0.)).z();

        geometryParameters.m_nRadLengthsInZGap      = 0;
        geometryParameters.m_nIntLengthsInZGap      = 0;
        geometryParameters.m_nRadLengthsInRadialGap = 0;
        geometryParameters.m_nIntLengthsInRadialGap = 0;

        const gear::CalorimeterParameters &eCalBarrelParameters = marlin::Global::GEAR->getEcalBarrelParameters();
        const gear::CalorimeterParameters &eCalEndCapParameters = marlin::Global::GEAR->getEcalEndcapParameters();
        const gear::CalorimeterParameters &hCalBarrelParameters = marlin::Global::GEAR->getHcalBarrelParameters();
        const gear::CalorimeterParameters &hCalEndCapParameters = marlin::Global::GEAR->getHcalEndcapParameters();

        // Initialize settings to gear defaults
        SetDefaultSubDetectorParameters(eCalBarrelParameters, geometryParameters.m_eCalBarrelParameters);
        SetDefaultSubDetectorParameters(eCalEndCapParameters, geometryParameters.m_eCalEndCapParameters);
        SetDefaultSubDetectorParameters(hCalBarrelParameters, geometryParameters.m_hCalBarrelParameters);
        SetDefaultSubDetectorParameters(hCalEndCapParameters, geometryParameters.m_hCalEndCapParameters);

        // Non-default values (and those missing from GEAR parameters file)...
        geometryParameters.m_eCalEndCapParameters.m_innerSymmetryOrder = 4;
        geometryParameters.m_hCalBarrelParameters.m_outerPhiCoordinate = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_phi0");
        geometryParameters.m_hCalBarrelParameters.m_outerSymmetryOrder = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_order");

        // Addition subdetectors
        this->SetAdditionalSubDetectorParameters(geometryParameters);

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::Create(*m_pPandora, geometryParameters));
    }
    catch (gear::UnknownParameterException &e)
    {
        streamlog_out(ERROR) << "Failed to extract geometry information from gear." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::SetDefaultSubDetectorParameters(const gear::CalorimeterParameters &inputParameters,
    PandoraApi::GeometryParameters::SubDetectorParameters &subDetectorParameters) const
{
    const gear::LayerLayout &layerLayout = inputParameters.getLayerLayout();

    subDetectorParameters.m_innerRCoordinate    = inputParameters.getExtent()[0];
    subDetectorParameters.m_innerZCoordinate    = inputParameters.getExtent()[2];
    subDetectorParameters.m_innerPhiCoordinate  = inputParameters.getPhi0();
    subDetectorParameters.m_innerSymmetryOrder  = inputParameters.getSymmetryOrder();
    subDetectorParameters.m_outerRCoordinate    = inputParameters.getExtent()[1];
    subDetectorParameters.m_outerZCoordinate    = inputParameters.getExtent()[3];
    subDetectorParameters.m_outerPhiCoordinate  = inputParameters.getPhi0();
    subDetectorParameters.m_outerSymmetryOrder  = inputParameters.getSymmetryOrder();
    subDetectorParameters.m_nLayers             = layerLayout.getNLayers();

    for(int i = 0; i < layerLayout.getNLayers(); ++i)
    {
        PandoraApi::Geometry::Parameters::LayerParameters layerParameters;
        layerParameters.m_closestDistanceToIp   = layerLayout.getDistance(i) + (0.5 * (layerLayout.getThickness(i) + layerLayout.getAbsorberThickness(i)));
        layerParameters.m_nRadiationLengths     = m_settings.m_absorberRadiationLength * layerLayout.getAbsorberThickness(i);
        layerParameters.m_nInteractionLengths   = m_settings.m_absorberInteractionLength * layerLayout.getAbsorberThickness(i);
        subDetectorParameters.m_layerParametersList.push_back(layerParameters);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::SetAdditionalSubDetectorParameters(PandoraApi::GeometryParameters &geometryParameters) const
{
    PandoraApi::Geometry::Parameters::SubDetectorParameters yokeBarrelParameters;
    const gear::CalorimeterParameters &yokeBarrelInputParameters = marlin::Global::GEAR->getYokeBarrelParameters();
    SetDefaultSubDetectorParameters(yokeBarrelInputParameters, yokeBarrelParameters);
    geometryParameters.m_additionalSubDetectors["YokeBarrel"] = yokeBarrelParameters;

    PandoraApi::Geometry::Parameters::SubDetectorParameters yokeEndcapParameters;
    const gear::CalorimeterParameters &yokeEndcapInputParameters = marlin::Global::GEAR->getYokeEndcapParameters();
    SetDefaultSubDetectorParameters(yokeEndcapInputParameters, yokeEndcapParameters);
    geometryParameters.m_additionalSubDetectors["YokeEndcap"] = yokeEndcapParameters;

    PandoraApi::Geometry::Parameters::SubDetectorParameters eCalPlugParameters;
    const gear::CalorimeterParameters &eCalPlugInputParameters = marlin::Global::GEAR->getEcalPlugParameters();
    SetDefaultSubDetectorParameters(eCalPlugInputParameters, eCalPlugParameters);
    geometryParameters.m_additionalSubDetectors["ECalPlug"] = eCalPlugParameters;

    PandoraApi::Geometry::Parameters::SubDetectorParameters hCalRingParameters;
    const gear::CalorimeterParameters &hCalRingInputParameters = marlin::Global::GEAR->getHcalRingParameters();
    SetDefaultSubDetectorParameters(hCalRingInputParameters, hCalRingParameters);
    geometryParameters.m_additionalSubDetectors["HCalRing"] = hCalRingParameters;

    PandoraApi::Geometry::Parameters::SubDetectorParameters lCalParameters;
    const gear::CalorimeterParameters &lCalInputParameters = marlin::Global::GEAR->getLcalParameters();
    SetDefaultSubDetectorParameters(lCalInputParameters, lCalParameters);
    geometryParameters.m_additionalSubDetectors["LCal"] = lCalParameters;

    PandoraApi::Geometry::Parameters::SubDetectorParameters lHCalParameters;
    const gear::CalorimeterParameters &lHCalInputParameters = marlin::Global::GEAR->getLHcalParameters();
    SetDefaultSubDetectorParameters(lHCalInputParameters, lHCalParameters);
    geometryParameters.m_additionalSubDetectors["LHCal"] = lHCalParameters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::RegisterUserAlgorithmFactories() const
{
    // Insert user code here ...

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateMCParticles(const LCEvent *const pLCEvent) const
{
    // Insert user code here ...
    for (StringVector::const_iterator iter = m_settings.m_mcParticleCollections.begin(), iterEnd = m_settings.m_mcParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pMCParticleCollection = pLCEvent->getCollection(*iter);

            for(int i = 0, iMax = pMCParticleCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    MCParticle *pMcParticle = dynamic_cast<MCParticle*>(pMCParticleCollection->getElementAt(i));

                    double innerRadius = 0.;
                    double outerRadius = 0.;
                    pandora::CartesianVector momentum( pMcParticle->getMomentum()[0], pMcParticle->getMomentum()[1], pMcParticle->getMomentum()[2] );

                    for(int i = 0; i < 3; ++i)
                    {
                        innerRadius += pow(pMcParticle->getVertex()[i], 2);
                        outerRadius += pow(pMcParticle->getEndpoint()[i], 2);
                    }

                    innerRadius = std::sqrt(innerRadius);
                    outerRadius = std::sqrt(outerRadius);
         
                    PandoraApi::MCParticle::Parameters mcParticleParameters;
                    mcParticleParameters.m_energy = pMcParticle->getEnergy();
                    mcParticleParameters.m_particleId = pMcParticle->getPDG();
                    mcParticleParameters.m_momentum = momentum;
                    mcParticleParameters.m_innerRadius = innerRadius;
                    mcParticleParameters.m_outerRadius = outerRadius;
                    mcParticleParameters.m_pParentAddress = pMcParticle;

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*m_pPandora, mcParticleParameters));

                    // Create parent-daughter relationships
                    for(MCParticleVec::const_iterator itDaughter = pMcParticle->getDaughters().begin(),
                        itDaughterEnd = pMcParticle->getDaughters().end(); itDaughter != itDaughterEnd; ++itDaughter)
                    {
                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*m_pPandora, pMcParticle,
                            *itDaughter));
                    }
                }
                catch (StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract MCParticle: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract MCParticle, unrecognised exception" << std::endl;
                }
            }
        }
        catch (...)
        {
            streamlog_out(WARNING) << "Failed to extract MCParticles collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateTrackAssociations(const LCEvent *const pLCEvent) const
{
    // Insert user code here ...
    for (StringVector::const_iterator iter = m_settings.m_v0VertexCollections.begin(), iterEnd = m_settings.m_v0VertexCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pV0Collection = pLCEvent->getCollection(*iter);

            for (int i = 0, iMax = pV0Collection->getNumberOfElements(); i < iMax; ++i)
            {
                Vertex *pVertex = dynamic_cast<Vertex*>(pV0Collection->getElementAt(i));

                ReconstructedParticle *pReconstructedParticle = pVertex->getAssociatedParticle();
                TrackVec trackVec = pReconstructedParticle->getTracks();
                for(unsigned int iTrack = 0; iTrack < trackVec.size(); ++iTrack)
                {
                    try
                    {
                        Track *pTrack = trackVec[iTrack];
                        TrackerHitVec trackerHitVec = pTrack->getTrackerHits();
                        const float nTrackHits(trackerHitVec.size());

                        streamlog_out(DEBUG) << "  V0Track " << iTrack
                                             << ", nTrackHits " << nTrackHits
                                             << ", ptrack " << pTrack << std::endl;
                    }
                    catch (...)
                    {
                        streamlog_out(WARNING) << "Failed to extract v0 vertex, unrecognised exception" << std::endl;
                    }
                }
            }
        }
        catch (...)
        {
            streamlog_out(WARNING) << "Failed to extract v0 vertex collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateTracks(const LCEvent *const pLCEvent)
{
    // Insert user code here ...
    m_trackVector.clear();

    for (StringVector::const_iterator iter = m_settings.m_trackCollections.begin(), 
        iterEnd = m_settings.m_trackCollections.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pTrackCollection = pLCEvent->getCollection(*iter);

            for (int i = 0, iMax = pTrackCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    Track *pTrack = dynamic_cast<Track*>(pTrackCollection->getElementAt(i));

                    const int nTrackHits(static_cast<int>(pTrack->getTrackerHits().size()));

                    if ((nTrackHits < m_settings.m_minTrackHits) || (nTrackHits > m_settings.m_maxTrackHits))
                        continue;

                    // Proceed to create the pandora track
                    PandoraApi::Track::Parameters trackParameters;
                    trackParameters.m_d0 = pTrack->getD0();
                    trackParameters.m_z0 = pTrack->getZ0();
                    trackParameters.m_pParentAddress = pTrack;

                    // For now, assume tracks are charged pions
                    trackParameters.m_mass = 0.13957018;

                    const float signedCurvature(pTrack->getOmega());

                    if (0. != signedCurvature)
                        trackParameters.m_chargeSign = static_cast<int>(signedCurvature / std::fabs(signedCurvature));

                    this->FitHelices(pTrack, trackParameters);
                    trackParameters.m_reachesECal = this->ReachesECAL(pTrack);

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(*m_pPandora, trackParameters));
                    m_trackVector.push_back(pTrack);
                }
                catch (StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract a track, unrecognised exception" << std::endl;
                }
            }
        }
        catch (...)
        {
            streamlog_out(WARNING) << "Failed to extract track collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::FitHelices(const Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    static const float bField(marlin::Global::GEAR->getBField().at(gear::Vector3D(0., 0., 0.)).z());

    // Fit from track parameters to determine momentum at dca
    HelixClass *pHelixFit = new HelixClass();
    pHelixFit->Initialize_Canonical(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(), pTrack->getTanLambda(), bField);

    trackParameters.m_momentumAtDca = pandora::CartesianVector(pHelixFit->getMomentum()[0], pHelixFit->getMomentum()[1], pHelixFit->getMomentum()[2]);

    // Fit start and end of tracks
    TrackerHitVec trackerHitvec(pTrack->getTrackerHits());
    const int nTrackHits = trackerHitvec.size();
    const int nTrackHitsForFit = std::min(m_settings.m_nHitsForHelixFits, nTrackHits);

    // Order hits by increasing z
    for (int iz = 0 ; iz < nTrackHits - 1; ++iz)
    {
        for (int jz = 0; jz < nTrackHits - iz - 1; ++jz)
        {
            if(trackerHitvec[jz]->getPosition()[2] > trackerHitvec[jz + 1]->getPosition()[2])
            {
                TrackerHit *pTempTrackerHit = trackerHitvec[jz];
                trackerHitvec[jz] = trackerHitvec[jz + 1];
                trackerHitvec[jz + 1] = pTempTrackerHit;
            }
        }
    }

    // Arrays for helix fits
    float xf[nTrackHitsForFit], yf[nTrackHitsForFit], zf[nTrackHitsForFit], rf[nTrackHitsForFit], af[nTrackHitsForFit];
    float xb[nTrackHitsForFit], yb[nTrackHitsForFit], zb[nTrackHitsForFit], rb[nTrackHitsForFit], ab[nTrackHitsForFit];

    for(int i = 0; i < nTrackHitsForFit; ++i)
    {
        xf[i] = trackerHitvec[i]->getPosition()[0];
        yf[i] = trackerHitvec[i]->getPosition()[1];
        zf[i] = trackerHitvec[i]->getPosition()[2];
        rf[i] = std::sqrt(xf[i] * xf[i] + yf[i] * yf[i]);
        af[i] = 0;

        int j = nTrackHits - 1 - i;
        xb[i] = trackerHitvec[j]->getPosition()[0];
        yb[i] = trackerHitvec[j]->getPosition()[1];
        zb[i] = trackerHitvec[j]->getPosition()[2];
        rb[i] = std::sqrt(xb[i] * xb[i] + yb[i] * yb[i]);
        ab[i] = 0;
    }

    // Find z extremes of track
    const float zMin(zf[0]);
    const float zMax(zb[0]);
    const float rMin(std::sqrt(xf[0] * xf[0] + yf[0] * yf[0]));
    const float rMax(std::sqrt(xb[0] * xb[0] + yb[0] * yb[0]));
    const int signPz(this->GetTrackSignPz(zMin, zMax, rMin, rMax));

    // Helix from first nTrackHitsForFit (i.e. lowest z)
    float par[5], dpar[5], chi2, distmax;
    ClusterShapes clusterShapesF(nTrackHitsForFit, af, xf, yf, zf);
    clusterShapesF.FitHelix(500, 0, 1, par, dpar, chi2, distmax);
    HelixClass *pHelix1 = new HelixClass();
    pHelix1->Initialize_BZ(par[0], par[1], par[2], par[3], par[4], bField, signPz, zf[0]);

    // Helix from last nTrackHitsForFit (i.e. highest z)
    ClusterShapes clusterShapesB(nTrackHitsForFit, ab, xb, yb, zb);
    clusterShapesB.FitHelix(500, 0, 1, par, dpar, chi2, distmax);
    HelixClass *pHelix2 = new HelixClass();
    pHelix2->Initialize_BZ(par[0], par[1], par[2], par[3], par[4], bField, signPz, zb[0]);

    // Label as start and end depending on assigned sign of Pz
    HelixClass *const pHelixStart = (signPz < 0) ? pHelix2 : pHelix1;
    HelixClass *const pHelixEnd   = (signPz < 0) ? pHelix1 : pHelix2;

    trackParameters.m_trackStateAtStart = pandora::TrackState(pHelixStart->getReferencePoint()[0], pHelixStart->getReferencePoint()[1],
        pHelixStart->getReferencePoint()[2], pHelixStart->getMomentum()[0], pHelixStart->getMomentum()[1], pHelixStart->getMomentum()[2]);

    trackParameters.m_trackStateAtEnd = pandora::TrackState(pHelixEnd->getReferencePoint()[0], pHelixEnd->getReferencePoint()[1],
        pHelixEnd->getReferencePoint()[2], pHelixEnd->getMomentum()[0], pHelixEnd->getMomentum()[1], pHelixEnd->getMomentum()[2]);

    // Get track state at ecal surface
    HelixClass* pHelixToProject = pHelixFit;

    if(0 != m_settings.m_useEndTrackHelixForECalProjection)
        pHelixToProject = pHelixEnd;

    float referencePoint[3] = {pHelixToProject->getReferencePoint()[0], pHelixToProject->getReferencePoint()[1],
        pHelixToProject->getReferencePoint()[2]};

    if(0 != m_settings.m_useDcaAsReferencePointForProjection)
    {
        const float trackPhi0(pTrack->getPhi());
        const float trackD0(pTrack->getD0());
        referencePoint[0] = ( trackD0 * sin(trackPhi0));
        referencePoint[1] = (-trackD0 * cos(trackPhi0));
        referencePoint[2] = pTrack->getZ0();
    }

    trackParameters.m_trackStateAtECal = this->GetECalProjection(pHelixToProject, referencePoint, signPz);

    streamlog_out(DEBUG) << "TrackStateAtStart: " << std::endl << trackParameters.m_trackStateAtStart.Get() << std::endl
                         << "TrackStateAtEnd: "   << std::endl << trackParameters.m_trackStateAtEnd.Get()   << std::endl
                         << "TrackStateAtECal: "  << std::endl << trackParameters.m_trackStateAtECal.Get()  << std::endl;

    delete pHelix1;
    delete pHelix2;
    delete pHelixFit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int PandoraPFANewProcessor::GetTrackSignPz(float zMin, float zMax, float rMin, float rMax) const
{
    // Need to decide track direction +ve/-ve in z. For now assume track originates from close to IP.
    // TODO: NEED TO ADD V0, KINK, PRONG, BACK-SCATTER INFORMATION HERE
    int signPz = 0;
    if((fabs(zMin) < fabs(zMax)) && (rMin < rMax))
    {
        signPz = +1;
    }
    else if((fabs(zMax) < fabs(zMin)) && (rMax < rMin))
    {
        signPz = -1;
    }

    // Tracks that cross TPC central plane
    if(zMin < 0 && zMax > 0)
    {
        if(rMin < rMax)
            signPz = +1;

        if(rMin > rMax)
            signPz = -1;

        // TODO: if this is a track from IP check momentumAtDca assignment is correct
        streamlog_out(WARNING) << "PandoraPFANewProcessor::GetTrackSignPz CROSSES TPC check code..." << std::endl;
    }

    // If above conditions not satisfied, default is to order in z
    if(0 == signPz)
    {
        if(fabs(zMin) < fabs(zMax))
            signPz = +1;

        if(fabs(zMin) > fabs(zMax))
            signPz = -1;

        // TODO: these tracks should be associated with a V0, etc... If not something could be wrong
        streamlog_out(WARNING) << "PandoraPFANewProcessor::GetTrackSignPz DEFAULT TRACK DIRECTION..." << std::endl;

        if(0 == signPz)
            throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    return signPz;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::TrackState PandoraPFANewProcessor::GetECalProjection(HelixClass *const pHelix, float referencePoint[3], int signPz) const
{
    static const gear::CalorimeterParameters &ecalBarrelParameters = marlin::Global::GEAR->getEcalBarrelParameters();
    static const gear::CalorimeterParameters &ecalEndCapParameters = marlin::Global::GEAR->getEcalEndcapParameters();

    static const float phi0(ecalBarrelParameters.getPhi0());
    static const int ecalSymmetryOrder(ecalBarrelParameters.getSymmetryOrder());
    static const float rOfBarrel(ecalBarrelParameters.getExtent()[0]);
    static const float zOfEndCap(ecalEndCapParameters.getExtent()[2]);

    float bestEcalProjection[3];

    // First project to endcap
    float minTime = pHelix->getPointInZ(static_cast<float>(signPz) * zOfEndCap, referencePoint, bestEcalProjection);

    // Then project to barrel surface(s)
    float barrelProjection[3];
    static const float pi(std::acos(-1.));

    if (ecalSymmetryOrder > 0)
    {
        // Polygon
        float twopi_n = 2. * pi / (static_cast<float>(ecalSymmetryOrder));

        for (int i = 0; i < ecalSymmetryOrder; ++i)
        {
            float phi = twopi_n * static_cast<float>(i) + phi0;
            float xx = rOfBarrel * cos(phi);
            float yy = rOfBarrel * sin(phi);
            float ax = cos(phi + 0.5*pi);
            float ay = sin(phi + 0.5*pi);
            float tt = pHelix->getPointInXY(xx, yy , ax, ay, referencePoint, barrelProjection);

            // If helix intersects this plane before current best use this point
            if (tt < minTime)
            {
                minTime = tt;
                bestEcalProjection[0] = barrelProjection[0];
                bestEcalProjection[1] = barrelProjection[1];
                bestEcalProjection[2] = barrelProjection[2];
            }
        }
    }
    else
    {
        // Cylinder
        float tt = pHelix->getPointOnCircle(rOfBarrel, referencePoint, barrelProjection);

        if (tt < minTime)
        {
            minTime = tt;
            bestEcalProjection[0] = barrelProjection[0];
            bestEcalProjection[1] = barrelProjection[1];
            bestEcalProjection[2] = barrelProjection[2];
        }
    }

    float extrapolatedMomentum[3];
    pHelix->getExtrapolatedMomentum(bestEcalProjection, extrapolatedMomentum);

    return pandora::TrackState(bestEcalProjection[0], bestEcalProjection[1], bestEcalProjection[2],
        extrapolatedMomentum[0], extrapolatedMomentum[1], extrapolatedMomentum[2]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PandoraPFANewProcessor::ReachesECAL(const Track *const pTrack)
{
    static const gear::TPCParameters &tpcParameters = marlin::Global::GEAR->getTPCParameters();
    static const gear::PadRowLayout2D &tpcPadLayout = tpcParameters.getPadLayout();
    static const float tpcInnerR(tpcPadLayout.getPlaneExtent()[0]);
    static const float tpcOuterR(tpcPadLayout.getPlaneExtent()[1]);
    static const float tpcZmax(tpcParameters.getMaxDriftLength());
    static const float tpcMaxRow(tpcPadLayout.getNRows());

    if (0 == tpcMaxRow)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    static const float tpcRowHeight((tpcOuterR - tpcInnerR) / tpcMaxRow);

    TrackerHitVec trackerHitVec(pTrack->getTrackerHits());
    const int nTrackHits(trackerHitVec.size());

    int nTpcOuter = 0, nTpcEnd = 0;

    for (int i = 0; i < nTrackHits; ++i)
    {
        const float x(trackerHitVec[i]->getPosition()[0]);
        const float y(trackerHitVec[i]->getPosition()[1]);
        const float z(trackerHitVec[i]->getPosition()[2]);
        const float r(std::sqrt(x * x + y * y));

        // HitTypes: 1 = vtx, 2 = etd/ftd, 4 = sit/set, 5 = tpc
        int hitType = trackerHitVec[i]->getType() / 100;

        if(hitType == 5)
        {
            if(r > (tpcOuterR - 20 * tpcRowHeight))
            {
                if (++nTpcOuter > 5)
                    return true;
            }

            if(fabs(z) > (tpcZmax - 20 * tpcRowHeight))
            {
                if (++nTpcEnd > 5)
                    return true;
            }
        }
        else if((hitType == 4 && r > tpcOuterR) || (hitType == 2 && fabs(z) > tpcZmax) || (hitType == 2 && trackerHitVec[i]->getType() > 204))
        {
            return true;
        }
    }

    streamlog_out(DEBUG) << " nTrackHits : " << nTrackHits
                         << " vtxHits : " << pTrack->getSubdetectorHitNumbers()[0]
                         << " ftdHits : " << pTrack->getSubdetectorHitNumbers()[1]
                         << " sitHits : " << pTrack->getSubdetectorHitNumbers()[2]
                         << " tpcHits : " << pTrack->getSubdetectorHitNumbers()[3]
                         << " nTpcOuter : " << nTpcOuter << " nTpcEnd : " << nTpcEnd << std::endl;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateTrackToMCParticleRelationships(const LCEvent *const pLCEvent) const
{
    typedef std::map<MCParticle *, float> MCParticleToWeightMap;
    MCParticleToWeightMap mcParticleToWeightMap;

    for (StringVector::const_iterator iter = m_settings.m_lcTrackRelationCollections.begin(), iterEnd = m_settings.m_lcTrackRelationCollections.end();
         iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pMCRelationCollection = pLCEvent->getCollection(*iter);
            LCRelationNavigator navigate(pMCRelationCollection);

            for (TrackVector::const_iterator trackIter = m_trackVector.begin(), trackIterEnd = m_trackVector.end();
                trackIter != trackIterEnd; ++trackIter)
            {
                try
                {
                    mcParticleToWeightMap.clear();

                    MCParticle* mcParticle = NULL;
                    
                    const LCObjectVec &objectVec = navigate.getRelatedToObjects(*trackIter);

                    if (objectVec.size() > 0) 
                    {
                        mcParticle = dynamic_cast<MCParticle*>(objectVec[0]);
                        mcParticleToWeightMap[mcParticle] += 1.0;
                    }

                    for (MCParticleToWeightMap::const_iterator mcParticleIter = mcParticleToWeightMap.begin(),
                        mcParticleIterEnd = mcParticleToWeightMap.end(); mcParticleIter != mcParticleIterEnd; ++mcParticleIter)
                    {
                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackToMCParticleRelationship(*m_pPandora,
                            *trackIter, mcParticleIter->first, mcParticleIter->second));
                    }
                }
                catch (StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract track to mc particle relationship: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract track to mc particle relationship, unrecognised exception" << std::endl;
                }
            }
        }
        catch(...)
        {
            streamlog_out(WARNING) << "Failed to extract track to mc particle relationships collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateCaloHits(const LCEvent *const pLCEvent)
{
    CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");
    m_calorimeterHitVector.clear();

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateHCalCaloHits(pLCEvent));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateECalCaloHits(const LCEvent *const pLCEvent)
{
    static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getEcalEndcapParameters().getLayerLayout());
    static const gear::LayerLayout &barrelLayerLayout(marlin::Global::GEAR->getEcalBarrelParameters().getLayerLayout());

    static const float endCapZCoordinate(marlin::Global::GEAR->getEcalEndcapParameters().getExtent()[2]);
    static const unsigned int barrelSymmetryOrder(marlin::Global::GEAR->getEcalBarrelParameters().getSymmetryOrder());
    static const float barrelPhi0(marlin::Global::GEAR->getEcalBarrelParameters().getPhi0());

    for (StringVector::const_iterator iter = m_settings.m_eCalCaloHitCollections.begin(),
        iterEnd = m_settings.m_eCalCaloHitCollections.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);

            for (int i = 0, iMax = pCaloHitCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::ECAL;
                    caloHitParameters.m_detectorRegion = (fabs(pCaloHit->getPosition()[2]) < endCapZCoordinate) ?
                        pandora::BARREL : pandora::ENDCAP;

                    this->GetCommonCaloHitProperties(pCaloHit, cellIdDecoder, caloHitParameters);

                    float absorberCorrection(1.);

                    if (pandora::ENDCAP == caloHitParameters.m_detectorRegion.Get())
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);
                    }
                    else
                    {
                        const unsigned int staveNumber(cellIdDecoder(pCaloHit)["S-1"]);
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, barrelSymmetryOrder, barrelPhi0, staveNumber,
                            caloHitParameters, absorberCorrection);
                    }

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_eCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
                        continue;

                    caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * pCaloHit->getEnergy();
                    caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract ecal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract ecal calo hit, unrecognised exception" << std::endl;
                }
            }
        }
        catch (...)
        {
            streamlog_out(WARNING) << "Failed to extract ecal calo hit collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateHCalCaloHits(const LCEvent *const pLCEvent)
{
    static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout());
    static const gear::LayerLayout &barrelLayerLayout(marlin::Global::GEAR->getHcalBarrelParameters().getLayerLayout());

    static const float endCapZCoordinate(marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[2]);
    static const unsigned int barrelSymmetryOrder(marlin::Global::GEAR->getHcalBarrelParameters().getSymmetryOrder());
    static const float barrelPhi0(marlin::Global::GEAR->getHcalBarrelParameters().getPhi0());

    for (StringVector::const_iterator iter = m_settings.m_hCalCaloHitCollections.begin(),
        iterEnd = m_settings.m_hCalCaloHitCollections.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);

            for (int i = 0, iMax = pCaloHitCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::HCAL;

                    this->GetCommonCaloHitProperties(pCaloHit, cellIdDecoder, caloHitParameters);
                    caloHitParameters.m_detectorRegion = (fabs(pCaloHit->getPosition()[2]) < endCapZCoordinate) ?
                        pandora::BARREL : pandora::ENDCAP;

                    float absorberCorrection(1.);

                    if (pandora::ENDCAP == caloHitParameters.m_detectorRegion.Get())
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);
                    }
                    else
                    {
                        const unsigned int staveNumber(cellIdDecoder(pCaloHit)["S-1"]);
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, barrelSymmetryOrder, barrelPhi0, staveNumber,
                            caloHitParameters, absorberCorrection);
                    }

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = m_settings.m_hCalToHadGeV * pCaloHit->getEnergy();
                    caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract hcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract hcal calo hit, unrecognised exception" << std::endl;
                }
            }
        }
        catch (...)
        {
            streamlog_out(WARNING) << "Failed to extract hcal calo hit collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::GetCommonCaloHitProperties(CalorimeterHit *const pCaloHit, CellIDDecoder<CalorimeterHit> &cellIdDecoder,
    PandoraApi::CaloHit::Parameters &caloHitParameters) const
{
    const float *pCaloHitPosition(pCaloHit->getPosition());
    caloHitParameters.m_positionVector = pandora::CartesianVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

    caloHitParameters.m_pParentAddress = pCaloHit;
    caloHitParameters.m_isDigital = false;

    caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
    caloHitParameters.m_time = pCaloHit->getTime();

    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)["K-1"];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::GetEndCapCaloHitProperties(CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
    PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const
{
    const unsigned int physicalLayer(caloHitParameters.m_layer.Get());

    caloHitParameters.m_cellSizeU = layerLayout.getCellSize0(physicalLayer);
    caloHitParameters.m_cellSizeV = layerLayout.getCellSize1(physicalLayer);
    caloHitParameters.m_cellThickness = layerLayout.getThickness(physicalLayer);

    const float layerAbsorberThickness(layerLayout.getAbsorberThickness(std::max(0, static_cast<int>(physicalLayer) - 1)));

    if (0 == layerAbsorberThickness)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    caloHitParameters.m_nRadiationLengths = m_settings.m_absorberRadiationLength * layerAbsorberThickness;
    caloHitParameters.m_nInteractionLengths = m_settings.m_absorberInteractionLength * layerAbsorberThickness;

    absorberCorrection = layerLayout.getAbsorberThickness(0) / layerAbsorberThickness;

    caloHitParameters.m_normalVector = (pCaloHit->getPosition()[2] > 0) ? pandora::CartesianVector(0, 0, 1) : pandora::CartesianVector(0, 0, -1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::GetBarrelCaloHitProperties(CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
    unsigned int barrelSymmetryOrder, float barrelPhi0, unsigned int staveNumber,
    PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const
{
    const unsigned int physicalLayer(caloHitParameters.m_layer.Get());

    caloHitParameters.m_cellSizeU = layerLayout.getCellSize0(physicalLayer);
    caloHitParameters.m_cellSizeV = layerLayout.getCellSize1(physicalLayer);
    caloHitParameters.m_cellThickness = layerLayout.getThickness(physicalLayer);

    const float layerAbsorberThickness(layerLayout.getAbsorberThickness(std::max(0, static_cast<int>(physicalLayer) - 1)));

    if (0 == layerAbsorberThickness)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    caloHitParameters.m_nRadiationLengths = m_settings.m_absorberRadiationLength * layerAbsorberThickness;
    caloHitParameters.m_nInteractionLengths = m_settings.m_absorberInteractionLength * layerAbsorberThickness;

    absorberCorrection = layerLayout.getAbsorberThickness(0) / layerAbsorberThickness;

    if (barrelSymmetryOrder > 0)
    {
        static const float pi(std::acos(-1.));
        const float phi = barrelPhi0 + (2. * pi * static_cast<float>(staveNumber) / static_cast<float>(barrelSymmetryOrder));
        caloHitParameters.m_normalVector = pandora::CartesianVector(-std::sin(phi), std::cos(phi), 0);
    }
    else
    {
        const float *pCaloHitPosition(pCaloHit->getPosition());

        if (pCaloHitPosition[1] != 0)
        {
            const float phi = barrelPhi0 + std::atan(pCaloHitPosition[0] / pCaloHitPosition[1]);
            caloHitParameters.m_normalVector = pandora::CartesianVector(std::sin(phi), std::cos(phi), 0);
        }
        else
        {
            caloHitParameters.m_normalVector = (pCaloHitPosition[0] > 0) ? pandora::CartesianVector(1, 0, 0) : pandora::CartesianVector(-1, 0, 0);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateCaloHitToMCParticleRelationships(const LCEvent *const pLCEvent) const
{
    typedef std::map<MCParticle *, float> MCParticleToEnergyWeightMap;
    MCParticleToEnergyWeightMap mcParticleToEnergyWeightMap;

    for (StringVector::const_iterator iter = m_settings.m_lcCaloHitRelationCollections.begin(), iterEnd = m_settings.m_lcCaloHitRelationCollections.end();
         iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pMCRelationCollection = pLCEvent->getCollection(*iter);
            LCRelationNavigator navigate(pMCRelationCollection);

            for (CalorimeterHitVector::const_iterator caloHitIter = m_calorimeterHitVector.begin(),
                caloHitIterEnd = m_calorimeterHitVector.end(); caloHitIter != caloHitIterEnd; ++caloHitIter)
            {
                try
                {
                    mcParticleToEnergyWeightMap.clear();
                    const LCObjectVec &objectVec = navigate.getRelatedToObjects(*caloHitIter);

                    for(LCObjectVec::const_iterator itRel = objectVec.begin(), itRelEnd = objectVec.end(); itRel != itRelEnd; ++itRel)
                    {
                        SimCalorimeterHit *pSimHit = dynamic_cast<SimCalorimeterHit *>(*itRel);

                        if(pSimHit == NULL)
                            return STATUS_CODE_FAILURE;

                        for(int iCont = 0, iEnd = pSimHit->getNMCContributions(); iCont < iEnd; ++iCont)
                            mcParticleToEnergyWeightMap[pSimHit->getParticleCont(iCont)] += pSimHit->getEnergyCont(iCont);
                    }

                    for (MCParticleToEnergyWeightMap::const_iterator mcParticleIter = mcParticleToEnergyWeightMap.begin(),
                        mcParticleIterEnd = mcParticleToEnergyWeightMap.end(); mcParticleIter != mcParticleIterEnd; ++mcParticleIter)
                    {
                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(*m_pPandora,
                            *caloHitIter, mcParticleIter->first, mcParticleIter->second));
                    }
                }
                catch (StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract calo hit to mc particle relationship: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract calo hit to mc particle relationship, unrecognised exception" << std::endl;
                }
            }
        }
        catch(...)
        {
            streamlog_out(WARNING) << "Failed to extract calo hit to mc particle relationships collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::ProcessParticleFlowObjects( LCEvent * pLCEvent)
{
    pandora::ParticleFlowObjectList particleFlowObjectList;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraApi::GetParticleFlowObjects(*m_pPandora,
        particleFlowObjectList));

    LCCollectionVec *pClusterCollection = new LCCollectionVec(LCIO::CLUSTER);
    LCCollectionVec *pReconstructedParticleCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

    LCFlagImpl lcFlagImpl(pClusterCollection->getFlag());
    lcFlagImpl.setBit(LCIO::CLBIT_HITS);
    pClusterCollection->setFlag(lcFlagImpl.getFlag());

    std::vector<std::string> subDetectorNames ;
    subDetectorNames.push_back("ecal") ; const unsigned int ecal_Index(0) ;
    subDetectorNames.push_back("hcal") ; const unsigned int hcal_Index(1) ;
    subDetectorNames.push_back("yoke") ; const unsigned int yoke_Index(2) ;
    subDetectorNames.push_back("lcal") ; const unsigned int lcal_Index(3) ;
    subDetectorNames.push_back("lhcal"); const unsigned int lhcal_Index(4);
    subDetectorNames.push_back("bcal") ; const unsigned int bcal_Index(5) ;

    pClusterCollection->parameters().setValues("ClusterSubdetectorNames" , subDetectorNames);

    // Create lcio "reconstructed particles" from the pandora "particle flow objects"
    for (pandora::ParticleFlowObjectList::iterator itPFO = particleFlowObjectList.begin(), itPFOEnd = particleFlowObjectList.end();
         itPFO != itPFOEnd; ++itPFO)
    {
        ReconstructedParticleImpl *pReconstructedParticle= new ReconstructedParticleImpl();

        pandora::ClusterAddressList clusterAddressList = (*itPFO)->GetClusterAddressList();
        pandora::TrackAddressList trackAddressList = (*itPFO)->GetTrackAddressList();

        // Create lcio clusters
        for (pandora::ClusterAddressList::iterator itCluster = clusterAddressList.begin(), itClusterEnd = clusterAddressList.end();
            itCluster != itClusterEnd; ++itCluster)
        {
            ClusterImpl *pCluster = new ClusterImpl();

            const unsigned int nHitsInCluster((*itCluster).size());

            float clusterEnergy(0.);
            float *pHitE = new float[nHitsInCluster];
            float *pHitX = new float[nHitsInCluster];
            float *pHitY = new float[nHitsInCluster];
            float *pHitZ = new float[nHitsInCluster];

            for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit)
            {
                CalorimeterHit *pCalorimeterHit = (CalorimeterHit*)((*itCluster)[iHit]);
                pCluster->addHit(pCalorimeterHit, 1.0);

                const float caloHitEnergy(pCalorimeterHit->getEnergy());
                clusterEnergy += caloHitEnergy;

                pHitE[iHit] = caloHitEnergy;
                pHitX[iHit] = pCalorimeterHit->getPosition()[0];
                pHitY[iHit] = pCalorimeterHit->getPosition()[1];
                pHitZ[iHit] = pCalorimeterHit->getPosition()[2];

                std::vector<float> &subDetectorEnergies = pCluster->subdetectorEnergies();
                subDetectorEnergies.resize(subDetectorNames.size());

                switch(CHT(pCalorimeterHit->getType()).caloID())
                {
                    case CHT::ecal:  subDetectorEnergies[ecal_Index ] += caloHitEnergy; break;
                    case CHT::hcal:  subDetectorEnergies[hcal_Index ] += caloHitEnergy; break;
                    case CHT::yoke:  subDetectorEnergies[yoke_Index ] += caloHitEnergy; break;
                    case CHT::lcal:  subDetectorEnergies[lcal_Index ] += caloHitEnergy; break;
                    case CHT::lhcal: subDetectorEnergies[lhcal_Index] += caloHitEnergy; break;
                    case CHT::bcal:  subDetectorEnergies[bcal_Index ] += caloHitEnergy; break;
                    default: streamlog_out(DEBUG) << " no subdetector found for hit with type: " << pCalorimeterHit->getType() << std::endl;
                }
            }

            pCluster->setEnergy(clusterEnergy);

            ClusterShapes *pClusterShapes = new ClusterShapes(nHitsInCluster, pHitE, pHitX, pHitY, pHitZ);
            pCluster->setPosition(pClusterShapes->getCentreOfGravity());
            pCluster->setIPhi(std::atan2(pClusterShapes->getEigenVecInertia()[1], pClusterShapes->getEigenVecInertia()[0]));
            pCluster->setITheta(std::acos(pClusterShapes->getEigenVecInertia()[2]));

            pClusterCollection->addElement(pCluster);
            pReconstructedParticle->addCluster(pCluster);

            delete pClusterShapes;
            delete[] pHitE; delete[] pHitX; delete[] pHitY; delete[] pHitZ;
        }

        // Add tracks to the lcio reconstructed particles
        for (pandora::TrackAddressList::iterator itTrack = trackAddressList.begin(), itTrackEnd = trackAddressList.end(); itTrack != itTrackEnd;
            ++itTrack)
        {
            pReconstructedParticle->addTrack((Track*)(*itTrack));
        }

        float momentum[3] = { (*itPFO)->GetMomentum().GetX(), (*itPFO)->GetMomentum().GetY(), (*itPFO)->GetMomentum().GetZ() };
        pReconstructedParticle->setMomentum(momentum);
        pReconstructedParticle->setEnergy((*itPFO)->GetEnergy());
        pReconstructedParticle->setMass((*itPFO)->GetMass());
        pReconstructedParticle->setCharge((*itPFO)->GetChargeSign());
        pReconstructedParticle->setType((*itPFO)->GetParticleId());

        pReconstructedParticleCollection->addElement(pReconstructedParticle);
    }

    pLCEvent->addCollection(pClusterCollection, m_settings.m_clusterCollectionName.c_str());
    pLCEvent->addCollection(pReconstructedParticleCollection, m_settings.m_pfoCollectionName.c_str());

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
                            "ECalCaloHitCollections", 
                            "Name of the ECAL calo hit collections",
                            m_settings.m_eCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::CALORIMETERHIT,
                            "HCalCaloHitCollections", 
                            "Name of the HCAL calo hit collections",
                            m_settings.m_hCalCaloHitCollections,
                            StringVector());

    registerInputCollections(LCIO::MCPARTICLE,
                            "MCParticleCollections", 
                            "Name of mc particle collections",
                            m_settings.m_mcParticleCollections,
                            StringVector());

    registerInputCollections(LCIO::LCRELATION, 
                            "RelCaloHitCollections",
                            "SimCaloHit to CaloHit Relations Collection Name",
                            m_settings.m_lcCaloHitRelationCollections,
                            StringVector());

    registerInputCollections(LCIO::LCRELATION, 
                            "RelTrackCollections",
                            "Track to MCParticle Relations Collection Name",
                            m_settings.m_lcTrackRelationCollections,
                            StringVector());

    // Absorber properties
    registerProcessorParameter("AbsorberRadiationLength",
                            "The absorber radation length",
                            m_settings.m_absorberRadiationLength,
                            float(1.));

    registerProcessorParameter("AbsorberInteractionLength",
                            "The absorber interaction length",
                            m_settings.m_absorberInteractionLength,
                            float(1.));

    // Name of PFO collection written by MarlinPandora
    registerOutputCollection( LCIO::CLUSTER,
                              "ClusterCollectionName" , 
                              "Cluster Collection Name "  ,
                              m_settings.m_clusterCollectionName,
                              std::string("PandoraPFANewClusters"));

    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                              "PFOCollectionName" , 
                              "PFO Collection Name "  ,
                              m_settings.m_pfoCollectionName,
                              std::string("PandoraPFANewPFOs"));

    // Calibration constants
    registerProcessorParameter("ECalToMipCalibration",
                            "The calibration from deposited ECal energy to mip",
                            m_settings.m_eCalToMip,
                            float(1.));

    registerProcessorParameter("HCalToMipCalibration",
                            "The calibration from deposited HCal energy to mip",
                            m_settings.m_hCalToMip,
                            float(1.));

    registerProcessorParameter("ECalMipThreshold",
                            "Threshold for creating calo hits in the ECal, units mip",
                            m_settings.m_eCalMipThreshold,
                            float(0.));

    registerProcessorParameter("HCalMipThreshold",
                            "Threshold for creating calo hits in the HCal, units mip",
                            m_settings.m_hCalMipThreshold,
                            float(0.));

    registerProcessorParameter("ECalToEMGeVCalibration",
                            "The calibration from deposited ECal energy to EM energy",
                            m_settings.m_eCalToEMGeV,
                            float(1.));

    registerProcessorParameter("HCalToEMGeVCalibration",
                            "The calibration from deposited HCal energy to EM energy",
                            m_settings.m_hCalToEMGeV,
                            float(1.));

    registerProcessorParameter("ECalToHadGeVCalibration",
                            "The calibration from deposited ECal energy to hadronic energy",
                            m_settings.m_eCalToHadGeV,
                            float(1.));

    registerProcessorParameter("HCalToHadGeVCalibration",
                            "The calibration from deposited HCal energy to hadronic energy",
                            m_settings.m_hCalToHadGeV,
                            float(1.));

    // For calculating track properties
   registerProcessorParameter("MinTrackHits",
                            "Track quality cut: the minimum number of track hits",
                            m_settings.m_minTrackHits,
                            int(5));

   registerProcessorParameter("MaxTrackHits",
                            "Track quality cut: the maximum number of track hits",
                            m_settings.m_maxTrackHits,
                            int(5000));

    registerProcessorParameter("NumberOfHitsForTrackHelixFits",
                            "The number of hits to be used in helix fits at start/end of tracks",
                            m_settings.m_nHitsForHelixFits,
                            int(50));

    registerProcessorParameter("UseEndTrackHelixForECalProjection",
                            "==0 use full track, ==1 use last NumberOfHitsForTrackHelixFit hits",
                            m_settings.m_useEndTrackHelixForECalProjection,
                            int(1));

    registerProcessorParameter("UseDcaForReferenceInECalProjection",
                            "==0 use helix reference point, ==1 use DCA as reference point",
                            m_settings.m_useDcaAsReferencePointForProjection,
                            int(1));
}
