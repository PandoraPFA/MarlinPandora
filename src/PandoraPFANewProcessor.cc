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
#include "EVENT/SimCalorimeterHit.h"

#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/ClusterImpl.h"

#include "UTIL/CellIDDecoder.h"
#include "UTIL/LCRelationNavigator.h"

#include "IMPL/LCCollectionVec.h"

#include "marlin/Global.h"

#include "gear/BField.h"
#include "gear/GEAR.h"
#include "gear/GearParameters.h"
#include "gear/CalorimeterParameters.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/LayerLayout.h"

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

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateGeometry());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RegisterUserAlgorithmFactories());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(m_pandora, m_settings.m_pandoraSettingsXmlFile));
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

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateTracks(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateTrackToMCParticleRelationships(pLCEvent));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateCaloHits(pLCEvent));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateCaloHitToMCParticleRelationships(pLCEvent));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(m_pandora));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ProcessParticleFlowObjects(pLCEvent));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(m_pandora));
    }
    catch (StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "Failed to process event: " << statusCodeException.ToString() << std::endl;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(m_pandora));
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Failed to process event, unrecognized exception" << std::endl;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(m_pandora));        
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

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::Create(m_pandora, geometryParameters));
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
    geometryParameters.m_additionalSubDetectors.push_back(yokeBarrelParameters);

    PandoraApi::Geometry::Parameters::SubDetectorParameters yokeEndcapParameters;
    const gear::CalorimeterParameters &yokeEndcapInputParameters = marlin::Global::GEAR->getYokeEndcapParameters();
    SetDefaultSubDetectorParameters(yokeEndcapInputParameters, yokeEndcapParameters);
    geometryParameters.m_additionalSubDetectors.push_back(yokeEndcapParameters);

    PandoraApi::Geometry::Parameters::SubDetectorParameters eCalPlugParameters;
    const gear::CalorimeterParameters &eCalPlugInputParameters = marlin::Global::GEAR->getEcalPlugParameters();
    SetDefaultSubDetectorParameters(eCalPlugInputParameters, eCalPlugParameters);
    geometryParameters.m_additionalSubDetectors.push_back(eCalPlugParameters);

    PandoraApi::Geometry::Parameters::SubDetectorParameters hCalRingParameters;
    const gear::CalorimeterParameters &hCalRingInputParameters = marlin::Global::GEAR->getHcalRingParameters();
    SetDefaultSubDetectorParameters(hCalRingInputParameters, hCalRingParameters);
    geometryParameters.m_additionalSubDetectors.push_back(hCalRingParameters);

    PandoraApi::Geometry::Parameters::SubDetectorParameters lCalParameters;
    const gear::CalorimeterParameters &lCalInputParameters = marlin::Global::GEAR->getLcalParameters();
    SetDefaultSubDetectorParameters(lCalInputParameters, lCalParameters);
    geometryParameters.m_additionalSubDetectors.push_back(lCalParameters);

    PandoraApi::Geometry::Parameters::SubDetectorParameters lHCalParameters;
    const gear::CalorimeterParameters &lHCalInputParameters = marlin::Global::GEAR->getLHcalParameters();
    SetDefaultSubDetectorParameters(lHCalInputParameters, lHCalParameters);
    geometryParameters.m_additionalSubDetectors.push_back(lHCalParameters);
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
                mcParticleParameters.m_particleId = pMcParticle->getPDG();
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
            streamlog_out(ERROR) << "Failed to extract MCParticles: " << statusCodeException.ToString() << std::endl;
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Failed to extract MCParticles, unrecognised exception" << std::endl;
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
            
            for (int i = 0; i < pTrackCollection->getNumberOfElements(); ++i)
            {
                Track* pTrack = dynamic_cast<Track*>(pTrackCollection->getElementAt(i));
                m_trackVector.push_back( pTrack );

                PandoraApi::Track::Parameters trackParameters;
                trackParameters.m_d0 = pTrack->getD0();
                trackParameters.m_z0 = pTrack->getZ0();
                trackParameters.m_pParentAddress = pTrack;

                this->FitHelices(pTrack, trackParameters);
                trackParameters.m_reachesECal = this->ReachedECAL(pTrack);

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(m_pandora, trackParameters));
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Failed to extract a track, unrecognised exception" << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::FitHelices(const Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    static const float bField(marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z());

    // Fit from track parameters to determine momentum at dca
    HelixClass *pHelixFit = new HelixClass();
    pHelixFit->Initialize_Canonical(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(), pTrack->getTanLambda(), bField);

    trackParameters.m_momentumAtDca = pandora::CartesianVector(pHelixFit->getMomentum()[0], pHelixFit->getMomentum()[1], pHelixFit->getMomentum()[2]);
    delete pHelixFit;

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

    streamlog_out(DEBUG) << "TrackStateAtStart: " << std::endl << trackParameters.m_trackStateAtStart.Get() << std::endl
                         << "TrackStateAtEnd: "   << std::endl << trackParameters.m_trackStateAtEnd.Get()   << std::endl;

    // Get track state at ecal surface
    this->ProjectTrackToECal(pHelixEnd, signPz, trackParameters);

    delete pHelix1;
    delete pHelix2;
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
        streamlog_out(DEBUG) << "PandoraPFANewProcessor::FitHelices CROSSES TPC check code..." << std::endl;
    }

    // If above conditions not satisfied, default is to order in z
    if(0 == signPz)
    {
        if(fabs(zMin) < fabs(zMax))
            signPz = +1;

        if(fabs(zMin) > fabs(zMax))
            signPz = -1;

        // TODO: these tracks should be associated with a V0, etc... If not something could be wrong
        streamlog_out(DEBUG) << "PandoraPFANewProcessor::FitHelices DEFAULT TRACK DIRECTION..." << std::endl;

        if(0 == signPz)
            throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    return signPz;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFANewProcessor::ProjectTrackToECal(HelixClass *const pHelixEnd, int signPz, PandoraApi::Track::Parameters &trackParameters) const
{
    static const gear::CalorimeterParameters &ecalBarrelParameters = marlin::Global::GEAR->getEcalBarrelParameters();
    static const gear::CalorimeterParameters &ecalEndCapParameters = marlin::Global::GEAR->getEcalEndcapParameters();

    static const float phi0(ecalBarrelParameters.getPhi0());
    static const int ecalSymmetryOrder(ecalBarrelParameters.getSymmetryOrder());
    static const float rOfBarrel(ecalBarrelParameters.getExtent()[0]);
    static const float zOfEndCap(static_cast<float>(signPz) * ecalEndCapParameters.getExtent()[2]);

    float bestEcalProjection[3];

    // First project to endcap
    float referencePoint[3] = {pHelixEnd->getReferencePoint()[0], pHelixEnd->getReferencePoint()[1], pHelixEnd->getReferencePoint()[2]};
    float minTime = pHelixEnd->getPointInZ(zOfEndCap, referencePoint, bestEcalProjection);

    // Then project to barrel surface(s)
    float barrelProjection[3];
    static const float pi = std::acos(-1.);

    if (ecalSymmetryOrder > 0)
    {
        // Polygon
        float twopi_n = 2. * pi/(static_cast<float>(ecalSymmetryOrder));

        for (int i = 0; i < ecalSymmetryOrder; ++i)
        {
            float phi = twopi_n * ((float)i) + phi0;
            float xx = rOfBarrel * cos(phi);
            float yy = rOfBarrel * sin(phi);
            float ax = cos(phi + 0.5*pi);
            float ay = sin(phi + 0.5*pi);
            float tt = pHelixEnd->getPointInXY(xx, yy , ax, ay, referencePoint, barrelProjection);

            // if helix intersects this plane before current best use this point
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
        float tt = pHelixEnd->getPointOnCircle(rOfBarrel, referencePoint, barrelProjection);

        if (tt < minTime)
        {
            minTime = tt;
            bestEcalProjection[0] = barrelProjection[0];
            bestEcalProjection[1] = barrelProjection[1];
            bestEcalProjection[2] = barrelProjection[2];
        }
    }

    float extrapolatedMomentum[3];
    pHelixEnd->getExtrapolatedMomentum(bestEcalProjection, extrapolatedMomentum);

    trackParameters.m_trackStateAtECal = pandora::TrackState(bestEcalProjection[0], bestEcalProjection[1], bestEcalProjection[2],
        extrapolatedMomentum[0], extrapolatedMomentum[1], extrapolatedMomentum[2]);

    streamlog_out(DEBUG) << "TrackStateAtECal: " << std::endl << trackParameters.m_trackStateAtECal.Get() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PandoraPFANewProcessor::ReachedECAL(const Track *const pTrack)
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
            if((r > tpcOuterR) - (20 * tpcRowHeight))
            {
                if (++nTpcOuter > 5)
                    return true;
            }

            if((fabs(z) > tpcZmax) - (20 * tpcRowHeight))
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

            for (TrackVector::const_iterator trackIter = m_trackVector.begin(),
                     trackIterEnd = m_trackVector.end(); trackIter != trackIterEnd; ++trackIter)
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
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackToMCParticleRelationship(m_pandora,
                        *trackIter, mcParticleIter->first, mcParticleIter->second));
                }
            }
        }
        catch(...)
        {
            streamlog_out(ERROR) << "Failed to extract track to mc particle relationships from collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::CreateCaloHits(const LCEvent *const pLCEvent)
{
    // Insert user code here ...
    m_calorimeterHitVector.clear();

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
                m_calorimeterHitVector.push_back(pCaloHit);

                PandoraApi::CaloHit::Parameters caloHitParameters;

                const float *pCaloHitPosition(pCaloHit->getPosition());
                caloHitParameters.m_positionVector = pandora::CartesianVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

                // TODO complete these parameters


                caloHitParameters.m_normalVector = pandora::CartesianVector(4, 5, 6);

                caloHitParameters.m_cellSizeU = 1;
                caloHitParameters.m_cellSizeV = 2;
                caloHitParameters.m_cellSizeZ = 3;

                caloHitParameters.m_nRadiationLengths = 4;
                caloHitParameters.m_nInteractionLengths = 5;

                caloHitParameters.m_isDigital = false;
                caloHitParameters.m_hitType = pandora::ECAL;
                caloHitParameters.m_detectorRegion = pandora::BARREL;

                caloHitParameters.m_time = pCaloHit->getTime();
                caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();

                caloHitParameters.m_layer = cellIdDecoder(pCaloHit)["K-1"];
                caloHitParameters.m_pParentAddress = pCaloHit;

                if (pandora::ECAL == caloHitParameters.m_hitType.Get())
                {
                    caloHitParameters.m_mipEquivalentEnergy = m_settings.m_eCalToMip * pCaloHit->getEnergy();

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
                        continue;

                    caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * pCaloHit->getEnergy();
                    caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeV * pCaloHit->getEnergy();
                }
                else if (pandora::HCAL == caloHitParameters.m_hitType.Get())
                {
                    caloHitParameters.m_mipEquivalentEnergy = m_settings.m_hCalToMip * pCaloHit->getEnergy();

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * pCaloHit->getEnergy();
                    caloHitParameters.m_hadronicEnergy = m_settings.m_hCalToHadGeV * pCaloHit->getEnergy();
                }

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(m_pandora, caloHitParameters));
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            streamlog_out(ERROR) << "Failed to extract a calo hit: " << statusCodeException.ToString() << std::endl;
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Failed to extract a calo hit, unrecognised exception" << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
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
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(m_pandora,
                        *caloHitIter, mcParticleIter->first, mcParticleIter->second));
                }
            }
        }
        catch(...)
        {
            streamlog_out(ERROR) << "Failed to extract calo hit to mc particle relationships from collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PandoraPFANewProcessor::ProcessParticleFlowObjects( LCEvent * pLCEvent)
{
    // get the particle flow objects
    pandora::ParticleFlowObjectList particleFlowObjectList;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraApi::GetParticleFlowObjects(m_pandora,
        particleFlowObjectList));

    LCCollectionVec *pReconstructedParticleCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

    // get particle flow objects and create "reconstructed particles"
    for (pandora::ParticleFlowObjectList::iterator itPFO = particleFlowObjectList.begin(), itPFOEnd = particleFlowObjectList.end();
         itPFO != itPFOEnd; ++itPFO)
    {
        ReconstructedParticleImpl *pReconstructedParticle= new ReconstructedParticleImpl();

        pandora::ClusterAddressList clusterAddressList = (*itPFO)->GetClusterAddressList();
        pandora::TrackAddressList trackAddressList = (*itPFO)->GetTrackAddressList();

        // make LCIO clusters
        for (pandora::ClusterAddressList::iterator itCluster = clusterAddressList.begin(), itClusterEnd = clusterAddressList.end();
            itCluster != itClusterEnd; ++itCluster)
        {
            ClusterImpl *pCluster = new ClusterImpl();
            for (pandora::CaloHitAddressList::iterator itHit = (*itCluster).begin(), itHitEnd = (*itCluster).end(); itHit != itHitEnd; ++itHit)
            {
                pCluster->addHit((CalorimeterHit*)(*itHit), 1.0); // transform from Uid (=void*) to a CalorimeterHit*
            }

            pReconstructedParticle->addCluster(pCluster);
        }

        // add tracks
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
                            "CaloHitCollections", 
                            "Name of the HCAL collection used to form clusters",
                            m_settings.m_caloHitCollections,
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
    registerProcessorParameter("NumberOfHitsForTrackHelixFits",
                            "The number of hits to be used in helix fits at start/end of tracks",
                            m_settings.m_nHitsForHelixFits,
                            int(50));
}
