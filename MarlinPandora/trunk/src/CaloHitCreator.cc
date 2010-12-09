/**
 *  @file   MarlinPandora/src/CaloHitCreator.cc
 * 
 *  @brief  Implementation of the calo hit creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "gear/GearParameters.h"
#include "gear/CalorimeterParameters.h"
#include "gear/GearDistanceProperties.h"
#include "gear/GearPointProperties.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/LayerLayout.h"

#include "UTIL/CellIDDecoder.h"

#include "CaloHitCreator.h"
#include "PathLengthCalculator.h"
#include "PandoraPFANewProcessor.h"

#include <cmath>
#include <limits>

CalorimeterHitVector CaloHitCreator::m_calorimeterHitVector;

CaloHitCreator::CaloHitCreator(const Settings &settings) :
    m_settings(settings),
    m_pPandora(PandoraPFANewProcessor::GetPandora()),
    m_pPathLengthCalculator(PathLengthCalculator::GetInstance()),
    m_eCalBarrelOuterZ(marlin::Global::GEAR->getEcalBarrelParameters().getExtent()[3]),
    m_hCalBarrelOuterZ(marlin::Global::GEAR->getHcalBarrelParameters().getExtent()[3]),
    m_muonBarrelOuterZ(marlin::Global::GEAR->getYokeBarrelParameters().getExtent()[3]),
    m_coilOuterR(marlin::Global::GEAR->getGearParameters("CoilParameters").getDoubleVal("Coil_cryostat_outer_radius")),
    m_eCalBarrelInnerPhi0(marlin::Global::GEAR->getEcalBarrelParameters().getPhi0()),
    m_eCalBarrelInnerSymmetry(marlin::Global::GEAR->getEcalBarrelParameters().getSymmetryOrder()),
    m_hCalBarrelInnerPhi0(marlin::Global::GEAR->getHcalBarrelParameters().getPhi0()),
    m_hCalBarrelInnerSymmetry(marlin::Global::GEAR->getHcalBarrelParameters().getSymmetryOrder()),
    m_muonBarrelInnerPhi0(marlin::Global::GEAR->getYokeBarrelParameters().getPhi0()),
    m_muonBarrelInnerSymmetry(marlin::Global::GEAR->getYokeBarrelParameters().getSymmetryOrder()),
    m_hCalEndCapOuterR(marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[1]),
    m_hCalEndCapOuterZ(marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[3]),
    m_hCalBarrelOuterR(marlin::Global::GEAR->getHcalBarrelParameters().getExtent()[1]),
    m_hCalBarrelOuterPhi0(marlin::Global::GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_phi0")),
    m_hCalBarrelOuterSymmetry(marlin::Global::GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_order"))

{
    const gear::LayerLayout &hCalEndCapLayerLayout(marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout());
    const gear::LayerLayout &hCalBarrelLayerLayout(marlin::Global::GEAR->getHcalBarrelParameters().getLayerLayout()); 

    m_hCalEndCapLayerThickness = hCalEndCapLayerLayout.getThickness(hCalEndCapLayerLayout.getNLayers() - 1);
    m_hCalBarrelLayerThickness = hCalBarrelLayerLayout.getThickness(hCalBarrelLayerLayout.getNLayers() - 1);

    if ((0.f == m_hCalEndCapLayerThickness) || (0.f == m_hCalBarrelLayerThickness))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitCreator::~CaloHitCreator()
{
    delete m_pPathLengthCalculator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    UTIL::CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateMuonCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateLCalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateLHCalCaloHits(pLCEvent));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateECalCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_eCalCaloHitCollections.begin(), iterEnd = m_settings.m_eCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

            static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getEcalEndcapParameters().getLayerLayout());
            static const gear::LayerLayout &barrelLayerLayout(marlin::Global::GEAR->getEcalBarrelParameters().getLayerLayout()); 

            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding((layerCodingString.find("K-1") == std::string::npos) ? "K" : "K-1");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::ECAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = false;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    if (std::fabs(pCaloHit->getPosition()[2]) < m_eCalBarrelOuterZ)
                    {
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, m_eCalBarrelInnerSymmetry, m_eCalBarrelInnerPhi0,
                            cellIdDecoder(pCaloHit)["S-1"], caloHitParameters, absorberCorrection);
                        caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVBarrel * pCaloHit->getEnergy();
                    }
                    else
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);
                        caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVEndCap * pCaloHit->getEnergy();
                    }

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_eCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
                        continue;

                    caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract ecal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract ecal calo hit: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract ecal calo hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateHCalCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_hCalCaloHitCollections.begin(), iterEnd = m_settings.m_hCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

            static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout());
            static const gear::LayerLayout &barrelLayerLayout(marlin::Global::GEAR->getHcalBarrelParameters().getLayerLayout());

            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding((layerCodingString.find("K-1") == std::string::npos) ? "K" : "K-1");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::HCAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = (this->GetNLayersFromEdge(pCaloHit) <= m_settings.m_nOuterSamplingLayers);
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    if (std::fabs(pCaloHit->getPosition()[2]) < m_hCalBarrelOuterZ)
                    {
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, m_hCalBarrelInnerSymmetry, m_hCalBarrelInnerPhi0,
                            m_hCalBarrelInnerSymmetry - int(cellIdDecoder(pCaloHit)["S-1"] / 2), caloHitParameters, absorberCorrection);
                    }
                    else
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);
                    }

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = std::min(m_settings.m_hCalToHadGeV * pCaloHit->getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
                    caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract hcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract hcal calo hit: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract hcal calo hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateMuonCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_muonCaloHitCollections.begin(), iterEnd = m_settings.m_muonCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

            static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getYokeEndcapParameters().getLayerLayout());
            static const gear::LayerLayout &barrelLayerLayout(marlin::Global::GEAR->getYokeBarrelParameters().getLayerLayout()); 
            static const gear::LayerLayout &plugLayerLayout(marlin::Global::GEAR->getYokePlugParameters().getLayerLayout());

            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding((layerCodingString.find("K-1") == std::string::npos) ? "K" : "K-1");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::MUON;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = true;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    const float radius(std::sqrt(pCaloHit->getPosition()[0] * pCaloHit->getPosition()[0] +
                        pCaloHit->getPosition()[1] * pCaloHit->getPosition()[1]));

                    const bool isWithinCoil(radius < m_coilOuterR);
                    const bool isInBarrelRegion(std::fabs(pCaloHit->getPosition()[2]) < m_muonBarrelOuterZ);

                    float absorberCorrection(1.);

                    if (isInBarrelRegion && isWithinCoil)
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, plugLayerLayout, caloHitParameters, absorberCorrection);
                    }
                    else if (isInBarrelRegion)
                    {
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, m_muonBarrelInnerSymmetry, m_muonBarrelInnerPhi0,
                            cellIdDecoder(pCaloHit)["S-1"], caloHitParameters, absorberCorrection);
                    }
                    else
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);
                    }

                    if (m_settings.m_muonDigitalHits > 0)
                    {
                        caloHitParameters.m_isDigital = true;
                        caloHitParameters.m_inputEnergy = m_settings.m_muonHitEnergy;
                        caloHitParameters.m_hadronicEnergy = m_settings.m_muonHitEnergy;
                        caloHitParameters.m_electromagneticEnergy = m_settings.m_muonHitEnergy;
                        caloHitParameters.m_mipEquivalentEnergy = 1.f;
                    }
                    else
                    {
                        caloHitParameters.m_isDigital = false;
                        caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
                        caloHitParameters.m_hadronicEnergy = pCaloHit->getEnergy();
                        caloHitParameters.m_electromagneticEnergy = pCaloHit->getEnergy();
                        caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_muonToMip;
                    }

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract muon hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract muon hit: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract muon hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateLCalCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_lCalCaloHitCollections.begin(), iterEnd = m_settings.m_lCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

            static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getLcalParameters().getLayerLayout()); 

            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding((layerCodingString.find("K-1") == std::string::npos) ? "K" : "K-1");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::ECAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = false;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);
                    this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_eCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
                        continue;

                    caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * pCaloHit->getEnergy();
                    caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVEndCap * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract lcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract lcal calo hit: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract lcal calo hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateLHCalCaloHits(const EVENT::LCEvent *const pLCEvent)
{
    for (StringVector::const_iterator iter = m_settings.m_lHCalCaloHitCollections.begin(), iterEnd = m_settings.m_lHCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

            static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getLHcalParameters().getLayerLayout());

            UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding((layerCodingString.find("K-1") == std::string::npos) ? "K" : "K-1");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    EVENT::CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::HCAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = (this->GetNLayersFromEdge(pCaloHit) <= m_settings.m_nOuterSamplingLayers);
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);
                    this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = std::min(m_settings.m_hCalToHadGeV * pCaloHit->getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
                    caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract lhcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract lhcal calo hit: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(MESSAGE) << "Failed to extract lhcal calo hit collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetCommonCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const
{
    const float *pCaloHitPosition(pCaloHit->getPosition());
    const pandora::CartesianVector positionVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

    caloHitParameters.m_positionVector = positionVector;
    caloHitParameters.m_expectedDirection = positionVector.GetUnitVector();
    caloHitParameters.m_pParentAddress = (void*)pCaloHit;
    caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
    caloHitParameters.m_time = pCaloHit->getTime();

    // TODO When available, use gear to calculate path length properties.
    // const gear::Vector3D positionIP3D(0.f, 0.f, 0.f);
    // const gear::Vector3D positionVector3D(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);
    // radiationLengthsFromIp = marlin::Global::GEAR->getDistanceProperties().getNRadlen(positionIP3D, positionVector3D);
    // interactionLengthsFromIp = marlin::Global::GEAR->getDistanceProperties().getNIntlen(positionIP3D, positionVector3D);

    float radiationLengthsFromIp(0.f), interactionLengthsFromIp(0.f);

    try
    {
        m_pPathLengthCalculator->GetPathLengths(pCaloHit, radiationLengthsFromIp, interactionLengthsFromIp);
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "Failed to calculate calo hit path lengths from ip: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }

    caloHitParameters.m_nRadiationLengthsFromIp = radiationLengthsFromIp;
    caloHitParameters.m_nInteractionLengthsFromIp = interactionLengthsFromIp;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetEndCapCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
    PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const
{
    caloHitParameters.m_detectorRegion = pandora::ENDCAP;

    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), layerLayout.getNLayers() - 1));
    caloHitParameters.m_cellSizeU = layerLayout.getCellSize0(physicalLayer);
    caloHitParameters.m_cellSizeV = layerLayout.getCellSize1(physicalLayer);
    caloHitParameters.m_cellThickness = layerLayout.getThickness(physicalLayer);

    const float layerAbsorberThickness(layerLayout.getAbsorberThickness(physicalLayer));
    caloHitParameters.m_nCellRadiationLengths = m_settings.m_absorberRadiationLength * layerAbsorberThickness;
    caloHitParameters.m_nCellInteractionLengths = m_settings.m_absorberInteractionLength * layerAbsorberThickness;

    absorberCorrection = 1.;
    for (unsigned int i = 0, iMax = layerLayout.getNLayers(); i < iMax; ++i)
    {
        const float absorberThickness(layerLayout.getAbsorberThickness(i));

        if (absorberThickness <= 0.)
            continue;

        if (layerAbsorberThickness > 0.)
            absorberCorrection = absorberThickness / layerAbsorberThickness;

        break;
    }

    caloHitParameters.m_cellNormalVector = (pCaloHit->getPosition()[2] > 0) ? pandora::CartesianVector(0, 0, 1) :
        pandora::CartesianVector(0, 0, -1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetBarrelCaloHitProperties(const EVENT::CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
    unsigned int barrelSymmetryOrder, float barrelPhi0, unsigned int staveNumber, PandoraApi::CaloHit::Parameters &caloHitParameters,
    float &absorberCorrection) const
{
    caloHitParameters.m_detectorRegion = pandora::BARREL;

    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), layerLayout.getNLayers() - 1));
    caloHitParameters.m_cellSizeU = layerLayout.getCellSize0(physicalLayer);
    caloHitParameters.m_cellSizeV = layerLayout.getCellSize1(physicalLayer);
    caloHitParameters.m_cellThickness = layerLayout.getThickness(physicalLayer);

    const float layerAbsorberThickness(layerLayout.getAbsorberThickness(physicalLayer));
    caloHitParameters.m_nCellRadiationLengths = m_settings.m_absorberRadiationLength * layerAbsorberThickness;
    caloHitParameters.m_nCellInteractionLengths = m_settings.m_absorberInteractionLength * layerAbsorberThickness;

    absorberCorrection = 1.;
    for (unsigned int i = 0, iMax = layerLayout.getNLayers(); i < iMax; ++i)
    {
        const float absorberThickness(layerLayout.getAbsorberThickness(i));

        if (absorberThickness <= 0.)
            continue;

        if (layerAbsorberThickness > 0.)
            absorberCorrection = absorberThickness / layerAbsorberThickness;

        break;
    }

    if (barrelSymmetryOrder > 2)
    {
        static const float pi(std::acos(-1.));
        const float phi = barrelPhi0 + (2. * pi * static_cast<float>(staveNumber) / static_cast<float>(barrelSymmetryOrder));
        caloHitParameters.m_cellNormalVector = pandora::CartesianVector(-std::sin(phi), std::cos(phi), 0);
    }
    else
    {
        const float *pCaloHitPosition(pCaloHit->getPosition());

        if (pCaloHitPosition[1] != 0)
        {
            const float phi = barrelPhi0 + std::atan(pCaloHitPosition[0] / pCaloHitPosition[1]);
            caloHitParameters.m_cellNormalVector = pandora::CartesianVector(std::sin(phi), std::cos(phi), 0);
        }
        else
        {
            caloHitParameters.m_cellNormalVector = (pCaloHitPosition[0] > 0) ? pandora::CartesianVector(1, 0, 0) :
                pandora::CartesianVector(-1, 0, 0);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int CaloHitCreator::GetNLayersFromEdge(const EVENT::CalorimeterHit *const pCaloHit) const
{
    // Calo hit coordinate calculations
    const float barrelMaximumRadius(this->GetMaximumRadius(pCaloHit, m_hCalBarrelOuterSymmetry, m_hCalBarrelOuterPhi0));
    const float endCapMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_hCalEndCapInnerSymmetryOrder, m_settings.m_hCalEndCapInnerPhiCoordinate));
    const float caloHitAbsZ(std::fabs(pCaloHit->getPosition()[2]));

    // Distance from radial outer
    float radialDistanceToEdge(std::numeric_limits<float>::max());

    if (caloHitAbsZ < m_eCalBarrelOuterZ)
    {
        radialDistanceToEdge = (m_hCalBarrelOuterR - barrelMaximumRadius) / m_hCalBarrelLayerThickness;
    }
    else
    {
        radialDistanceToEdge = (m_hCalEndCapOuterR - endCapMaximumRadius) / m_hCalEndCapLayerThickness;
    }

    // Distance from rear of endcap outer
    float rearDistanceToEdge(std::numeric_limits<float>::max());

    if (caloHitAbsZ >= m_eCalBarrelOuterZ)
    {
        rearDistanceToEdge = (m_hCalEndCapOuterZ - caloHitAbsZ) / m_hCalEndCapLayerThickness;
    }
    else
    {
        const float rearDistance((m_eCalBarrelOuterZ - caloHitAbsZ) / m_hCalBarrelLayerThickness);

        if (rearDistance < m_settings.m_layersFromEdgeMaxRearDistance)
        {
            const float overlapDistance((m_hCalEndCapOuterR - endCapMaximumRadius) / m_hCalEndCapLayerThickness);
            rearDistanceToEdge = std::max(rearDistance, overlapDistance);
        }
    }

    return static_cast<int>(std::min(radialDistanceToEdge, rearDistanceToEdge));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CaloHitCreator::GetMaximumRadius(const EVENT::CalorimeterHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const
{
    const float *pCaloHitPosition(pCaloHit->getPosition());

    if (symmetryOrder <= 2)
        return std::sqrt((pCaloHitPosition[0] * pCaloHitPosition[0]) + (pCaloHitPosition[1] * pCaloHitPosition[1]));

    float maximumRadius(0.f);
    static const float twoPi = static_cast<float>(2. * std::acos(-1.));

    for (unsigned int i = 0; i < symmetryOrder; ++i)
    {
        const float phi = phi0 + i * twoPi / static_cast<float>(symmetryOrder);
        float radius = pCaloHitPosition[0] * std::cos(phi) + pCaloHitPosition[1] * std::sin(phi);

        if (radius > maximumRadius)
            maximumRadius = radius;
    }

    return maximumRadius;
}
