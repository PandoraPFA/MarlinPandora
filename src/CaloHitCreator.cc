/**
 *  @file   PandoraPFANew/src/CaloHitCreator.cc
 * 
 *  @brief  Implementation of the calo hit creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "gear/CalorimeterParameters.h"

#include "CaloHitCreator.h"
#include "PandoraPFANewProcessor.h"

#include <cmath>

CalorimeterHitVector CaloHitCreator::m_calorimeterHitVector;

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitCreator::CreateCaloHits(const LCEvent *const pLCEvent)
{
    CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits(pLCEvent));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreateHCalCaloHits(pLCEvent));

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
            CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);

            for (int i = 0, iMax = pCaloHitCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::ECAL;
                    caloHitParameters.m_detectorRegion = (fabs(pCaloHit->getPosition()[2]) < endCapZCoordinate) ? pandora::BARREL : pandora::ENDCAP;
                    caloHitParameters.m_isInOuterSamplingLayer = false;

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
            streamlog_out(WARNING) << "Failed to extract ecal calo hit collection: " << *iter << std::endl;
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
            CellIDDecoder<CalorimeterHit> cellIdDecoder(pCaloHitCollection);

            for (int i = 0, iMax = pCaloHitCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    CalorimeterHit *pCaloHit = dynamic_cast<CalorimeterHit*>(pCaloHitCollection->getElementAt(i));

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::HCAL;
                    caloHitParameters.m_detectorRegion = (fabs(pCaloHit->getPosition()[2]) < endCapZCoordinate) ? pandora::BARREL : pandora::ENDCAP;
                    caloHitParameters.m_isInOuterSamplingLayer = (this->GetNLayersFromEdge(pCaloHit) <= m_settings.m_nOuterSamplingLayers);

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
            streamlog_out(WARNING) << "Failed to extract hcal calo hit collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetCommonCaloHitProperties(CalorimeterHit *const pCaloHit, CellIDDecoder<CalorimeterHit> &cellIdDecoder,
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
