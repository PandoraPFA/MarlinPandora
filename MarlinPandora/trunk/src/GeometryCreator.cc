/**
 *  @file   PandoraPFANew/src/GeometryCreator.cc
 * 
 *  @brief  Implementation of the geometry creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "gear/BField.h"
#include "gear/GEAR.h"
#include "gear/GearParameters.h"
#include "gear/CalorimeterParameters.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/LayerLayout.h"

#include "GeometryCreator.h"
#include "PandoraPFANewProcessor.h"

StatusCode GeometryCreator::CreateGeometry() const
{
    try
    {
        // Insert user code here ...
        static pandora::Pandora *pPandora = PandoraPFANewProcessor::GetPandora();
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
        geometryParameters.m_eCalEndCapParameters.m_innerSymmetryOrder = m_settings.m_eCalEndCapInnerSymmetryOrder;
        geometryParameters.m_eCalEndCapParameters.m_innerPhiCoordinate = m_settings.m_eCalEndCapInnerPhiCoordinate;
        geometryParameters.m_hCalEndCapParameters.m_innerSymmetryOrder = m_settings.m_hCalEndCapInnerSymmetryOrder;
        geometryParameters.m_hCalEndCapParameters.m_innerPhiCoordinate = m_settings.m_hCalEndCapInnerPhiCoordinate;
        geometryParameters.m_hCalBarrelParameters.m_outerPhiCoordinate = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_phi0");
        geometryParameters.m_hCalBarrelParameters.m_outerSymmetryOrder = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_order");

        // Additional subdetectors
        this->SetAdditionalSubDetectorParameters(geometryParameters);

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::Create(*pPandora, geometryParameters));
    }
    catch (gear::UnknownParameterException &e)
    {
        streamlog_out(ERROR) << "Failed to extract geometry information from gear." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GeometryCreator::SetDefaultSubDetectorParameters(const gear::CalorimeterParameters &inputParameters,
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

    for (int i = 0; i < layerLayout.getNLayers(); ++i)
    {
        PandoraApi::Geometry::Parameters::LayerParameters layerParameters;
        layerParameters.m_closestDistanceToIp   = layerLayout.getDistance(i) + (0.5 * (layerLayout.getThickness(i) + layerLayout.getAbsorberThickness(i)));
        layerParameters.m_nRadiationLengths     = m_settings.m_absorberRadiationLength * layerLayout.getAbsorberThickness(i);
        layerParameters.m_nInteractionLengths   = m_settings.m_absorberInteractionLength * layerLayout.getAbsorberThickness(i);
        subDetectorParameters.m_layerParametersList.push_back(layerParameters);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GeometryCreator::SetAdditionalSubDetectorParameters(PandoraApi::GeometryParameters &geometryParameters) const
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
