/**
 *  @file   PandoraPFANew/src/CaloHitCreator.cc
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
#include "InteractionLengthCalculator.h"
#include "PandoraPFANewProcessor.h"

#include <cmath>
#include <limits>

CalorimeterHitVector CaloHitCreator::m_calorimeterHitVector;

CaloHitCreator::CaloHitCreator(const Settings &settings) :
    m_pInteractionLengthCalculator(InteractionLengthCalculator::GetInstance()),
    m_settings(settings)
{
    try
    {
        m_pPandora = PandoraPFANewProcessor::GetPandora();

        // ECal parameters
        m_eCalEndCapInnerZ = marlin::Global::GEAR->getEcalEndcapParameters().getExtent()[2];
        m_eCalBarrelOuterZ = marlin::Global::GEAR->getEcalBarrelParameters().getExtent()[3];
        m_eCalBarrelInnerPhi0 = marlin::Global::GEAR->getEcalBarrelParameters().getPhi0();
        m_eCalBarrelInnerSymmetry = marlin::Global::GEAR->getEcalBarrelParameters().getSymmetryOrder();

        // HCal parameters
        m_hCalEndCapInnerZ = marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[2];
        m_hCalEndCapOuterR = marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[1];
        m_hCalEndCapOuterZ = marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[3];

        m_hCalBarrelOuterR = marlin::Global::GEAR->getHcalBarrelParameters().getExtent()[1];
        m_hCalBarrelInnerPhi0 = marlin::Global::GEAR->getHcalBarrelParameters().getPhi0();
        m_hCalBarrelInnerSymmetry = marlin::Global::GEAR->getHcalBarrelParameters().getSymmetryOrder();
        m_hCalBarrelOuterPhi0 = marlin::Global::GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_phi0");
        m_hCalBarrelOuterSymmetry = marlin::Global::GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_order");

        const gear::LayerLayout &hCalEndCapLayerLayout(marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout());
        const gear::LayerLayout &hCalBarrelLayerLayout(marlin::Global::GEAR->getHcalBarrelParameters().getLayerLayout()); 

        m_hCalEndCapLayerThickness = hCalEndCapLayerLayout.getThickness(hCalEndCapLayerLayout.getNLayers() - 1);
        m_hCalBarrelLayerThickness = hCalBarrelLayerLayout.getThickness(hCalBarrelLayerLayout.getNLayers() - 1);

        if ((0.f == m_hCalEndCapLayerThickness) || (0.f == m_hCalBarrelLayerThickness))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        // Muon yoke parameters
        m_muonEndCapInnerZ = marlin::Global::GEAR->getYokeEndcapParameters().getExtent()[2];
        m_muonBarrelInnerPhi0 = marlin::Global::GEAR->getYokeBarrelParameters().getPhi0();
        m_muonBarrelInnerSymmetry = marlin::Global::GEAR->getYokeBarrelParameters().getSymmetryOrder();
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Failed to initialize calo hit creator." << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitCreator::~CaloHitCreator()
{
    delete m_pInteractionLengthCalculator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitCreator::CreateCaloHits(const LCEvent *const pLCEvent)
{
    CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateHCalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateMuonCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateLCalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateLHCalCaloHits(pLCEvent));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitCreator::CreateECalCaloHits(const LCEvent *const pLCEvent)
{
    static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getEcalEndcapParameters().getLayerLayout());
    static const gear::LayerLayout &barrelLayerLayout(marlin::Global::GEAR->getEcalBarrelParameters().getLayerLayout()); 

    for (StringVector::const_iterator iter = m_settings.m_eCalCaloHitCollections.begin(), iterEnd = m_settings.m_eCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

            CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding((layerCodingString.find("K-1") == std::string::npos) ? "K" : "K-1");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::ECAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_detectorRegion = (std::fabs(pCaloHit->getPosition()[2]) < m_eCalEndCapInnerZ) ? pandora::BARREL : pandora::ENDCAP;
                    caloHitParameters.m_isInOuterSamplingLayer = false;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    if (pandora::ENDCAP == caloHitParameters.m_detectorRegion.Get())
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);
                    }
                    else
                    {
                        const unsigned int staveNumber(cellIdDecoder(pCaloHit)["S-1"]);
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, m_eCalBarrelInnerSymmetry, m_eCalBarrelInnerPhi0,
                            staveNumber, caloHitParameters, absorberCorrection);
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
            streamlog_out(MESSAGE) << "Failed to extract ecal calo hit collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitCreator::CreateHCalCaloHits(const LCEvent *const pLCEvent)
{
    static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout());
    static const gear::LayerLayout &barrelLayerLayout(marlin::Global::GEAR->getHcalBarrelParameters().getLayerLayout());

    for (StringVector::const_iterator iter = m_settings.m_hCalCaloHitCollections.begin(), iterEnd = m_settings.m_hCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

            CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding((layerCodingString.find("K-1") == std::string::npos) ? "K" : "K-1");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::HCAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_detectorRegion = (std::fabs(pCaloHit->getPosition()[2]) < m_hCalEndCapInnerZ) ? pandora::BARREL : pandora::ENDCAP;
                    caloHitParameters.m_isInOuterSamplingLayer = (this->GetNLayersFromEdge(pCaloHit) <= m_settings.m_nOuterSamplingLayers);
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    if (pandora::ENDCAP == caloHitParameters.m_detectorRegion.Get())
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);
                    }
                    else
                    {
                        // TODO Proper fix for barrel normal vectors
                        const unsigned int staveNumber(cellIdDecoder(pCaloHit)["S-1"]);
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, m_hCalBarrelInnerSymmetry, m_hCalBarrelInnerPhi0,
                            m_hCalBarrelInnerSymmetry - int(staveNumber / 2), caloHitParameters, absorberCorrection);
                    }

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = std::min(m_settings.m_hCalToHadGeV * pCaloHit->getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
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
            streamlog_out(MESSAGE) << "Failed to extract hcal calo hit collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitCreator::CreateMuonCaloHits(const LCEvent *const pLCEvent)
{
    static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getYokeEndcapParameters().getLayerLayout());
    static const gear::LayerLayout &barrelLayerLayout(marlin::Global::GEAR->getYokeBarrelParameters().getLayerLayout()); 

    for (StringVector::const_iterator iter = m_settings.m_muonCaloHitCollections.begin(), iterEnd = m_settings.m_muonCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

            CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding((layerCodingString.find("K-1") == std::string::npos) ? "K" : "K-1");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::MUON;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_detectorRegion = (std::fabs(pCaloHit->getPosition()[2]) < m_muonEndCapInnerZ) ? pandora::BARREL : pandora::ENDCAP;
                    caloHitParameters.m_isInOuterSamplingLayer = true;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    if (pandora::ENDCAP == caloHitParameters.m_detectorRegion.Get())
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);
                    }
                    else
                    {
                        const unsigned int staveNumber(cellIdDecoder(pCaloHit)["S-1"]);
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, m_muonBarrelInnerSymmetry, m_muonBarrelInnerPhi0,
                            staveNumber, caloHitParameters, absorberCorrection);
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

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract muon hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract muon hit, unrecognised exception" << std::endl;
                }
            }
        }
        catch (...)
        {
            streamlog_out(MESSAGE) << "Failed to extract muon hit collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitCreator::CreateLCalCaloHits(const LCEvent *const pLCEvent)
{
    static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getLcalParameters().getLayerLayout()); 

    for (StringVector::const_iterator iter = m_settings.m_lCalCaloHitCollections.begin(), iterEnd = m_settings.m_lCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

            CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding((layerCodingString.find("K-1") == std::string::npos) ? "K" : "K-1");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::ECAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_detectorRegion = pandora::ENDCAP;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = false;

                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);
                    this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);

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
                    streamlog_out(ERROR) << "Failed to extract lcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract lcal calo hit, unrecognised exception" << std::endl;
                }
            }
        }
        catch (...)
        {
            streamlog_out(MESSAGE) << "Failed to extract lcal calo hit collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitCreator::CreateLHCalCaloHits(const LCEvent *const pLCEvent)
{
    static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getLHcalParameters().getLayerLayout());

    for (StringVector::const_iterator iter = m_settings.m_lHCalCaloHitCollections.begin(), iterEnd = m_settings.m_lHCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pCaloHitCollection = pLCEvent->getCollection(*iter);
            const int nElements(pCaloHitCollection->getNumberOfElements());

            if (0 == nElements)
                continue;

            CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);
            const std::string layerCodingString(pCaloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
            const std::string layerCoding((layerCodingString.find("K-1") == std::string::npos) ? "K" : "K-1");

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::HCAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_detectorRegion = pandora::ENDCAP;
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

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(pCaloHit);
                }
                catch (StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract lhcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract lhcal calo hit, unrecognised exception" << std::endl;
                }
            }
        }
        catch (...)
        {
            streamlog_out(MESSAGE) << "Failed to extract lhcal calo hit collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetCommonCaloHitProperties(CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const
{
    const float *pCaloHitPosition(pCaloHit->getPosition());
    caloHitParameters.m_positionVector = pandora::CartesianVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

    caloHitParameters.m_pParentAddress = pCaloHit;

    caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
    caloHitParameters.m_time = pCaloHit->getTime();

    float interactionLengthsFromIp = 0.f;

    try
    {
        const gear::Vector3D positionIP(0,0,0);
        const float* pos = pCaloHit->getPosition();
        const gear::Vector3D positionHit(pos[0], pos[1], pos[2] );
        interactionLengthsFromIp = marlin::Global::GEAR->getDistanceProperties().getNIntlen(positionIP, positionHit);
    }
    catch (gear::Exception &excpt)
    {
        m_pInteractionLengthCalculator->ComputeInteractionLengthsFromIP(pCaloHit, interactionLengthsFromIp);
    }

    caloHitParameters.m_nInteractionLengthsFromIp = interactionLengthsFromIp;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetEndCapCaloHitProperties(CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
    PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const
{
    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), layerLayout.getNLayers() - 1));

    caloHitParameters.m_cellSizeU = layerLayout.getCellSize0(physicalLayer);
    caloHitParameters.m_cellSizeV = layerLayout.getCellSize1(physicalLayer);
    caloHitParameters.m_cellThickness = layerLayout.getThickness(physicalLayer);

    const float layerAbsorberThickness(layerLayout.getAbsorberThickness(std::max(0, physicalLayer - 1)));

    if (0 == layerAbsorberThickness)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    caloHitParameters.m_nRadiationLengths = m_settings.m_absorberRadiationLength * layerAbsorberThickness;
    caloHitParameters.m_nInteractionLengths = m_settings.m_absorberInteractionLength * layerAbsorberThickness;

    absorberCorrection = layerLayout.getAbsorberThickness(0) / layerAbsorberThickness;

    caloHitParameters.m_normalVector = (pCaloHit->getPosition()[2] > 0) ? pandora::CartesianVector(0, 0, 1) : pandora::CartesianVector(0, 0, -1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetBarrelCaloHitProperties(CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
    unsigned int barrelSymmetryOrder, float barrelPhi0, unsigned int staveNumber, PandoraApi::CaloHit::Parameters &caloHitParameters,
    float &absorberCorrection) const
{
    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), layerLayout.getNLayers() - 1));

    caloHitParameters.m_cellSizeU = layerLayout.getCellSize0(physicalLayer);
    caloHitParameters.m_cellSizeV = layerLayout.getCellSize1(physicalLayer);
    caloHitParameters.m_cellThickness = layerLayout.getThickness(physicalLayer);

    const float layerAbsorberThickness(layerLayout.getAbsorberThickness(std::max(0, physicalLayer - 1)));

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

int CaloHitCreator::GetNLayersFromEdge(CalorimeterHit *const pCaloHit) const
{
    // Calo hit coordinate calculations
    const float barrelMaximumRadius(this->GetMaximumRadius(pCaloHit, m_hCalBarrelOuterSymmetry, m_hCalBarrelOuterPhi0));
    const float endCapMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_hCalEndCapInnerSymmetryOrder, m_settings.m_hCalEndCapInnerPhiCoordinate));
    const float caloHitAbsZ(std::fabs(pCaloHit->getPosition()[2]));

    // Distance from radial outer
    float radialDistanceToEdge(std::numeric_limits<float>::max());

    if (caloHitAbsZ < m_eCalEndCapInnerZ)
    {
        radialDistanceToEdge = (m_hCalBarrelOuterR - barrelMaximumRadius) / m_hCalBarrelLayerThickness;
    }
    else
    {
        radialDistanceToEdge = (m_hCalEndCapOuterR - endCapMaximumRadius) / m_hCalEndCapLayerThickness;
    }

    // Distance from rear of endcap outer
    float rearDistanceToEdge(std::numeric_limits<float>::max());

    if (caloHitAbsZ >= m_eCalEndCapInnerZ)
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

float CaloHitCreator::GetMaximumRadius(CalorimeterHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const
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
