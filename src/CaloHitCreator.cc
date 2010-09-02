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
#include "PandoraPFANewProcessor.h"

#include <cmath>
#include <limits>

CalorimeterHitVector CaloHitCreator::m_calorimeterHitVector;

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitCreator::CreateCaloHits(const LCEvent *const pLCEvent)
{
    CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateHCalCaloHits(pLCEvent));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateLCalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateLHCalCaloHits(pLCEvent));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateMuonCaloHits(pLCEvent));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitCreator::CreateECalCaloHits(const LCEvent *const pLCEvent)
{
    static pandora::Pandora *pPandora = PandoraPFANewProcessor::GetPandora();

    static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getEcalEndcapParameters().getLayerLayout());
    static const gear::LayerLayout &barrelLayerLayout(marlin::Global::GEAR->getEcalBarrelParameters().getLayerLayout());

    static const float endCapZCoordinate(marlin::Global::GEAR->getEcalEndcapParameters().getExtent()[2]);
    static const unsigned int barrelSymmetryOrder(marlin::Global::GEAR->getEcalBarrelParameters().getSymmetryOrder());
    static const float barrelPhi0(marlin::Global::GEAR->getEcalBarrelParameters().getPhi0());

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
                    caloHitParameters.m_detectorRegion = (std::fabs(pCaloHit->getPosition()[2]) < endCapZCoordinate) ? pandora::BARREL : pandora::ENDCAP;
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
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, barrelSymmetryOrder, barrelPhi0, staveNumber,
                            caloHitParameters, absorberCorrection);
                    }

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_eCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
                        continue;

                    caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * pCaloHit->getEnergy();
                    caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPandora, caloHitParameters));
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
    static pandora::Pandora *pPandora = PandoraPFANewProcessor::GetPandora();

    static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout());
    static const gear::LayerLayout &barrelLayerLayout(marlin::Global::GEAR->getHcalBarrelParameters().getLayerLayout());

    static const float endCapZCoordinate(marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[2]);
    static const unsigned int barrelSymmetryOrder(marlin::Global::GEAR->getHcalBarrelParameters().getSymmetryOrder());
    static const float barrelPhi0(marlin::Global::GEAR->getHcalBarrelParameters().getPhi0());

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
                    caloHitParameters.m_detectorRegion = (std::fabs(pCaloHit->getPosition()[2]) < endCapZCoordinate) ? pandora::BARREL : pandora::ENDCAP;
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
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, barrelSymmetryOrder, barrelPhi0, barrelSymmetryOrder - int(staveNumber / 2),
                            caloHitParameters, absorberCorrection);
                    }

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = std::min(m_settings.m_hCalToHadGeV * pCaloHit->getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
                    caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPandora, caloHitParameters));
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

StatusCode CaloHitCreator::CreateLCalCaloHits(const LCEvent *const pLCEvent)
{
    static pandora::Pandora *pPandora = PandoraPFANewProcessor::GetPandora();
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

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPandora, caloHitParameters));
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
    static pandora::Pandora *pPandora = PandoraPFANewProcessor::GetPandora();
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

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPandora, caloHitParameters));
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

StatusCode CaloHitCreator::CreateMuonCaloHits(const LCEvent *const pLCEvent)
{
    static pandora::Pandora *pPandora = PandoraPFANewProcessor::GetPandora();

    static const gear::LayerLayout &endcapLayerLayout(marlin::Global::GEAR->getYokeEndcapParameters().getLayerLayout());
    static const gear::LayerLayout &barrelLayerLayout(marlin::Global::GEAR->getYokeBarrelParameters().getLayerLayout());

    static const float endCapZCoordinate(marlin::Global::GEAR->getYokeEndcapParameters().getExtent()[2]);
    static const unsigned int barrelSymmetryOrder(marlin::Global::GEAR->getYokeBarrelParameters().getSymmetryOrder());
    static const float barrelPhi0(marlin::Global::GEAR->getYokeBarrelParameters().getPhi0());

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
                    caloHitParameters.m_detectorRegion = (std::fabs(pCaloHit->getPosition()[2]) < endCapZCoordinate) ? pandora::BARREL : pandora::ENDCAP;
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
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, barrelSymmetryOrder, barrelPhi0, staveNumber,
                            caloHitParameters, absorberCorrection);
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

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPandora, caloHitParameters));
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

void CaloHitCreator::GetCommonCaloHitProperties(CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const
{
    const float *pCaloHitPosition(pCaloHit->getPosition());
    caloHitParameters.m_positionVector = pandora::CartesianVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

    caloHitParameters.m_pParentAddress = pCaloHit;

    caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
    caloHitParameters.m_time = pCaloHit->getTime();

    double interactionLengthsFromIp = 0.f;
    try
    {
        const gear::Vector3D positionIP(0,0,0);
        const float* pos = pCaloHit->getPosition();
        const gear::Vector3D positionHit(pos[0], pos[1], pos[2] );
        interactionLengthsFromIp = marlin::Global::GEAR->getDistanceProperties().getNIntlen( positionIP, positionHit);
    }
    catch( gear::Exception excpt )
    {
        ComputeInteractionLengthsFromIP(pCaloHit, interactionLengthsFromIp);
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
    // Extract geometry details
    static const float hCalBarrelOuterRadius = marlin::Global::GEAR->getHcalBarrelParameters().getExtent()[1];
    static const float hCalEndCapOuterRadius = marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[1];
    static const float hCalEndCapOuterZCoordinate = marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[3];
    static const float eCalBarrelOuterZCoordinate = marlin::Global::GEAR->getEcalBarrelParameters().getExtent()[3];
    static const float eCalEndCapInnerZCoordinate = marlin::Global::GEAR->getEcalEndcapParameters().getExtent()[2];

    // Extract layer thicknesses
    static const gear::LayerLayout &hCalBarrelLayerLayout = marlin::Global::GEAR->getHcalBarrelParameters().getLayerLayout();
    static const gear::LayerLayout &hCalEndCapLayerLayout = marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout();
    static const float hCalBarrelLayerThickness = hCalBarrelLayerLayout.getThickness(hCalBarrelLayerLayout.getNLayers() - 1);
    static const float hCalEndCapLayerThickness = hCalEndCapLayerLayout.getThickness(hCalEndCapLayerLayout.getNLayers() - 1);

    // Calo hit coordinate calculations
    static const int hCalBarrelSymmetry = marlin::Global::GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_order");
    static const float hCalBarrelPhi0 = marlin::Global::GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_phi0");

    const float barrelMaximumRadius(this->GetMaximumRadius(pCaloHit, hCalBarrelSymmetry, hCalBarrelPhi0));
    const float endCapMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_hCalEndCapInnerSymmetryOrder, m_settings.m_hCalEndCapInnerPhiCoordinate));
    const float caloHitAbsZ(std::fabs(pCaloHit->getPosition()[2]));

    // Distance from radial outer
    float radialDistanceToEdge(std::numeric_limits<float>::max());

    if (caloHitAbsZ < eCalEndCapInnerZCoordinate)
    {
        radialDistanceToEdge = (hCalBarrelOuterRadius - barrelMaximumRadius) / hCalBarrelLayerThickness;
    }
    else
    {
        radialDistanceToEdge = (hCalEndCapOuterRadius - endCapMaximumRadius) / hCalEndCapLayerThickness;
    }

    // Distance from rear of endcap outer
    float rearDistanceToEdge(std::numeric_limits<float>::max());

    if (caloHitAbsZ >= eCalEndCapInnerZCoordinate)
    {
        rearDistanceToEdge = (hCalEndCapOuterZCoordinate - caloHitAbsZ) / hCalEndCapLayerThickness;
    }
    else
    {
        const float rearDistance((eCalBarrelOuterZCoordinate - caloHitAbsZ) / hCalBarrelLayerThickness);

        if (rearDistance < m_settings.m_layersFromEdgeMaxRearDistance)
        {
            const float overlapDistance((hCalEndCapOuterRadius - endCapMaximumRadius) / hCalEndCapLayerThickness);
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



//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitCreator::ComputeInteractionLengthsFromIP(CalorimeterHit *const& pCaloHit, double& lengthInUnitsOfInteractionLength) const
{
    try
    {
        // Insert user code here ...
        static bool coordinatesAlreadyComputed = false;

        // coordinates of the sub-detectors in one quadrant
        static float rMinECalBarrel = 0.f;
        static float rMaxECalBarrel = 0.f;
        static float zMinECalBarrel = 0.f;
        static float zMaxECalBarrel = 0.f;

        static float rMinECalEndCap = 0.f;
        static float rMaxECalEndCap = 0.f;
        static float zMinECalEndCap = 0.f;
        static float zMaxECalEndCap = 0.f;

        static float rMinHCalBarrel = 0.f;
        static float rMaxHCalBarrel = 0.f;
        static float zMinHCalBarrel = 0.f;
        static float zMaxHCalBarrel = 0.f;

        static float rMinHCalEndCap = 0.f;
        static float rMaxHCalEndCap = 0.f;
        static float zMinHCalEndCap = 0.f;
        static float zMaxHCalEndCap = 0.f;

        static float rMinCoil = 0.f;
        static float rMaxCoil = 0.f;
        static float zMinCoil = 0.f;
        static float zMaxCoil = 0.f;

        static float rMinTracker = 0.f;
        static float rMaxTracker = 0.f;
        static float zMinTracker = 0.f;
        static float zMaxTracker = 0.f;

        static float rMinMuonBarrel = 0.f;
        static float rMaxMuonBarrel = 0.f;
        static float zMinMuonBarrel = 0.f;
        static float zMaxMuonBarrel = 0.f;

        static float rMinMuonEndCap = 0.f;
        static float rMaxMuonEndCap = 0.f;
        static float zMinMuonEndCap = 0.f;
        static float zMaxMuonEndCap = 0.f;


        if( !coordinatesAlreadyComputed )
        {
            const PandoraApi::Geometry::Parameters geometryParameters;

            const gear::TPCParameters &tpcParameters    = marlin::Global::GEAR->getTPCParameters();
            const gear::PadRowLayout2D &tpcPadLayout    = tpcParameters.getPadLayout();
            rMinTracker = tpcPadLayout.getPlaneExtent()[0];
            rMaxTracker = tpcPadLayout.getPlaneExtent()[1];
            zMinTracker = 0;
            zMaxTracker     = tpcParameters.getMaxDriftLength();

            const gear::GearParameters &coilParameters  = marlin::Global::GEAR->getGearParameters("CoilParameters");
            rMinCoil        = coilParameters.getDoubleVal("Coil_cryostat_inner_radius");
            rMaxCoil        = coilParameters.getDoubleVal("Coil_cryostat_outer_radius");
            zMinCoil        = 0.f;
            zMaxCoil        = coilParameters.getDoubleVal("Coil_cryostat_half_z");

            const gear::CalorimeterParameters &hCalBarrelParameters = marlin::Global::GEAR->getHcalBarrelParameters();
            const gear::CalorimeterParameters &hCalEndCapParameters = marlin::Global::GEAR->getHcalEndcapParameters();
            const gear::CalorimeterParameters &eCalBarrelParameters = marlin::Global::GEAR->getEcalBarrelParameters();
            const gear::CalorimeterParameters &eCalEndCapParameters = marlin::Global::GEAR->getEcalEndcapParameters();
            const gear::CalorimeterParameters &muonBarrelParameters = marlin::Global::GEAR->getYokeBarrelParameters();
            const gear::CalorimeterParameters &muonEndCapParameters = marlin::Global::GEAR->getYokeEndcapParameters();


            rMinECalBarrel =  eCalBarrelParameters.getExtent()[0];
            rMaxECalBarrel =  eCalBarrelParameters.getExtent()[1];
            zMinECalBarrel =  eCalBarrelParameters.getExtent()[2];
            zMaxECalBarrel =  eCalBarrelParameters.getExtent()[3];

            rMinHCalBarrel =  hCalBarrelParameters.getExtent()[0];
            rMaxHCalBarrel =  hCalBarrelParameters.getExtent()[1];
            zMinHCalBarrel =  hCalBarrelParameters.getExtent()[2];
            zMaxHCalBarrel =  hCalBarrelParameters.getExtent()[3];

            rMinECalEndCap =  eCalEndCapParameters.getExtent()[0];
            rMaxECalEndCap =  eCalEndCapParameters.getExtent()[1];
            zMinECalEndCap =  eCalEndCapParameters.getExtent()[2];
            zMaxECalEndCap =  eCalEndCapParameters.getExtent()[3];

            rMinHCalEndCap =  hCalEndCapParameters.getExtent()[0];
            rMaxHCalEndCap =  hCalEndCapParameters.getExtent()[1];
            zMinHCalEndCap =  hCalEndCapParameters.getExtent()[2];
            zMaxHCalEndCap =  hCalEndCapParameters.getExtent()[3];

            rMinMuonBarrel =  muonBarrelParameters.getExtent()[0];
            rMaxMuonBarrel =  muonBarrelParameters.getExtent()[1];
            zMinMuonBarrel =  muonBarrelParameters.getExtent()[2];
            zMaxMuonBarrel =  muonBarrelParameters.getExtent()[3];
            
            rMinMuonEndCap =  muonEndCapParameters.getExtent()[0];
            rMaxMuonEndCap =  muonEndCapParameters.getExtent()[1];
            zMinMuonEndCap =  muonEndCapParameters.getExtent()[2];
            zMaxMuonEndCap =  muonEndCapParameters.getExtent()[3];

            coordinatesAlreadyComputed = true;
        }

        const float* position = pCaloHit->getPosition();
        pandora::CartesianVector pPosition( position[0], position[1], position[2] );
        float radius, phi, z;
        pPosition.GetCylindricalCoordinates(radius, phi, z);
        pPosition.SetValues( radius, 0, fabs(z) );

        
//        lengthInUnitsOfInteractionLength += ComputePathLengthFromIPInRectangle( pPosition, rMinTracker, zMinTracker, rMaxTracker, zMaxTracker ) * m_settings.avgIntLengthTracker;
        lengthInUnitsOfInteractionLength += ComputePathLengthFromIPInRectangle( pPosition, rMinECalBarrel, zMinECalBarrel, rMaxECalBarrel, zMaxECalBarrel ) * m_settings.avgIntLengthECalBarrel;
        lengthInUnitsOfInteractionLength += ComputePathLengthFromIPInRectangle( pPosition, rMinHCalBarrel, zMinHCalBarrel, rMaxHCalBarrel, zMaxHCalBarrel ) * m_settings.avgIntLengthHCalBarrel;
        lengthInUnitsOfInteractionLength += ComputePathLengthFromIPInRectangle( pPosition, rMinCoil, zMinCoil, rMaxCoil, zMaxCoil ) * m_settings.avgIntLengthCoil;
        lengthInUnitsOfInteractionLength += ComputePathLengthFromIPInRectangle( pPosition, rMinECalEndCap, zMinECalEndCap, rMaxECalEndCap, zMaxECalEndCap ) * m_settings.avgIntLengthECalEndCap;
        lengthInUnitsOfInteractionLength += ComputePathLengthFromIPInRectangle( pPosition, rMinHCalEndCap, zMinHCalEndCap, rMaxHCalEndCap, zMaxHCalEndCap ) * m_settings.avgIntLengthHCalEndCap;

        lengthInUnitsOfInteractionLength += ComputePathLengthFromIPInRectangle( pPosition, rMinMuonBarrel, zMinMuonBarrel, rMaxMuonBarrel, zMaxMuonBarrel ) * m_settings.avgIntLengthMuonBarrel;
        lengthInUnitsOfInteractionLength += ComputePathLengthFromIPInRectangle( pPosition, rMinMuonEndCap, zMinMuonEndCap, rMaxMuonEndCap, zMaxMuonEndCap ) * m_settings.avgIntLengthMuonEndCap;
    }
    catch (gear::UnknownParameterException &e)
    {
        streamlog_out(ERROR) << "Failed to extract geometry information from gear." << std::endl;
        return STATUS_CODE_FAILURE;
    }
    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

float CaloHitCreator::ComputePathLengthFromIPInRectangle(const pandora::CartesianVector& pPosition, 
                                                         const float& rMin, const float& zMin, const float& rMax, const float& zMax) const
{
    // compute cuts with rectangle borders
    float phi, radius, z;
    pPosition.GetCylindricalCoordinates( radius, phi, z );

    float xInt[4];
    float zInt[4];
    bool valid[4];
    valid[0] = IntersectLines2D( 0, 0, radius, z,   rMin, zMin,   rMin, zMax,   (xInt[0]), (zInt[0]) ); // first edge of rectangle at rMin
    valid[1] = IntersectLines2D( 0, 0, radius, z,   rMin, zMin,   rMax, zMin,   (xInt[1]), (zInt[1]) ); // first edge of rectangle at zMin
    valid[2] = IntersectLines2D( 0, 0, radius, z,   rMax, zMax,   rMax, zMin,   (xInt[2]), (zInt[2]) ); // first edge of rectangle at rMax
    valid[3] = IntersectLines2D( 0, 0, radius, z,   rMax, zMax,   rMin, zMax,   (xInt[3]), (zInt[3]) ); // first edge of rectangle at zMax

    int indexFirstPoint = -1;
    int indexSecondPoint = -1;

    for( int i = 0; i < 4; ++i )
    {
//        std::cout << "valid["<<i<<"] " << valid[i] << "  x " << xInt[i] << "  y " << zInt[i] << std::endl;
        if( valid[i] )
        {
            if( indexFirstPoint == -1 ) // first point not yet set
                indexFirstPoint = i;
            else if( indexSecondPoint == -1 ) // second point not yet set
                indexSecondPoint = i;
            else
                std::cout << "ERROR calohitcreator problem at computing path length from IP to calohit in rectangle" << std::endl;
        }
    }

    if( indexFirstPoint == -1 )
        return 0.f;

    pandora::CartesianVector intersectionA( xInt[indexFirstPoint], 0.f, zInt[indexFirstPoint] );
    if( indexSecondPoint == -1 )
    {
        float length = pandora::CartesianVector(pPosition - intersectionA).GetMagnitude();
        return length;
    }

    pandora::CartesianVector intersectionB( xInt[indexSecondPoint], 0.f, zInt[indexSecondPoint] );
    float length = pandora::CartesianVector(intersectionA-intersectionB).GetMagnitude();
    return length;
}


//------------------------------------------------------------------------------------------------------------------------------------------

bool CaloHitCreator::IntersectLines2D( const float& lineAXStart, const float& lineAYStart, const float& lineAXEnd, const float& lineAYEnd, 
                                       const float& lineBXStart, const float& lineBYStart, const float& lineBXEnd, const float& lineBYEnd, 
                                       float& xIntersect, float& yIntersect) const
{
    float k0,k1;    // the slopes of the two lines

    const float epsilon = 1e-6;

    bool parallelToY_A = false;
    bool parallelToY_B = false;

    if ( fabs(lineAXEnd-lineAXStart) > epsilon )
        k0 = (lineAYEnd-lineAYStart)/(lineAXEnd-lineAXStart);
    else
    {
        parallelToY_A = true;
        k0 = std::numeric_limits<float>::max();   // take max float value instead of infinity
    }

    if ( fabs(lineBXEnd-lineBXStart) > epsilon )
        k1 = (lineBYEnd-lineBYStart)/(lineBXEnd-lineBXStart);
    else
    {
        parallelToY_B = true;
        k1 = std::numeric_limits<float>::max();   // take max float value instead of infinity
    }

    if( parallelToY_A && parallelToY_B )
    {
//         if( lineAXStart == lineBXStart ) 
        xIntersect = 0.f;
        yIntersect = 0.f;
        return false;
    }       

    if( parallelToY_A )
    {
        xIntersect = lineAXStart;
        yIntersect = lineBXStart+k1*(xIntersect-lineBXStart);
    }
    else if( parallelToY_B )
    {
        xIntersect = lineBXStart;
        yIntersect = lineAXStart+k0*(xIntersect-lineAXStart);
    }
    else
    {

        // compute constants
        const float b0 = -1;
        const float b1 = -1;

        const float c0 = (lineAYStart-k0*lineAXStart);
        const float c1 = (lineBYStart-k1*lineBXStart);

        // compute the inverse of the determinate

        const float inverseDet = 1/(k0*b1 - k1*b0);

        // use Kramers rule to compute xi and yi
        xIntersect=((b0*c1 - b1*c0)*inverseDet);
        yIntersect=((k1*c0 - k0*c1)*inverseDet);
    }

    // check if intersections are within end-points of lines
    // check x coordinate, line A
    if( !(  (xIntersect >= lineAXStart && xIntersect <= lineAXEnd) || (xIntersect >= lineAXEnd && xIntersect <= lineAXStart) ) )
        return false;

    // check x coordinate, line B
    if( !(  (xIntersect >= lineBXStart && xIntersect <= lineBXEnd) || (xIntersect >= lineBXEnd && xIntersect <= lineBXStart) ) )
        return false;

    // check y coordinate, line A
    if( !(  (yIntersect >= lineAYStart && yIntersect <= lineAYEnd) || (yIntersect >= lineAYEnd && yIntersect <= lineAYStart) ) )
        return false;

    // check y coordinate, line B
    if( !(  (yIntersect >= lineBYStart && yIntersect <= lineBYEnd) || (yIntersect >= lineBYEnd && yIntersect <= lineBYStart) ) )
        return false;

    return true;
} 

