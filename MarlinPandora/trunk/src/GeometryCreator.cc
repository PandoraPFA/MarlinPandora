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

GeometryCreator::GeometryCreator(const Settings &settings) :
    m_settings(settings),
    m_pPandora(PandoraPFANewProcessor::GetPandora())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

GeometryCreator::~GeometryCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::CreateGeometry() const
{
    try
    {
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

        const gear::CalorimeterParameters &hCalBarrelParameters = marlin::Global::GEAR->getHcalBarrelParameters();
        const gear::CalorimeterParameters &hCalEndCapParameters = marlin::Global::GEAR->getHcalEndcapParameters();

        // Initialize settings to gear defaults
        SetDefaultSubDetectorParameters(marlin::Global::GEAR->getEcalBarrelParameters(), geometryParameters.m_eCalBarrelParameters);
        SetDefaultSubDetectorParameters(marlin::Global::GEAR->getEcalEndcapParameters(), geometryParameters.m_eCalEndCapParameters);
        SetDefaultSubDetectorParameters(hCalBarrelParameters, geometryParameters.m_hCalBarrelParameters);
        SetDefaultSubDetectorParameters(hCalEndCapParameters, geometryParameters.m_hCalEndCapParameters);
        SetDefaultSubDetectorParameters(marlin::Global::GEAR->getYokeBarrelParameters(), geometryParameters.m_muonBarrelParameters);
        SetDefaultSubDetectorParameters(marlin::Global::GEAR->getYokeEndcapParameters(), geometryParameters.m_muonEndCapParameters);

        // Additional subdetectors
        this->SetAdditionalSubDetectorParameters(geometryParameters);

        // Set positions of gaps in ILD detector and add information missing from GEAR parameters file
        if (std::string::npos != marlin::Global::GEAR->getDetectorName().find("ILD"))
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->SetILDSpecificGeometry(geometryParameters));
        }

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::Create(*m_pPandora, geometryParameters));
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(ERROR) << "Failure in marlin pandora geometry creator, gear exception: " << exception.what() << std::endl;
        throw exception;
    }

    return pandora::STATUS_CODE_SUCCESS;
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
    try
    {
        PandoraApi::Geometry::Parameters::SubDetectorParameters eCalPlugParameters;
        const gear::CalorimeterParameters &eCalPlugInputParameters = marlin::Global::GEAR->getEcalPlugParameters();
        SetDefaultSubDetectorParameters(eCalPlugInputParameters, eCalPlugParameters);
        geometryParameters.m_additionalSubDetectors["ECalPlug"] = eCalPlugParameters;
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(WARNING) << "Marlin pandora geometry creator: " << exception.what() << std::endl;
    }

    try
    {
        PandoraApi::Geometry::Parameters::SubDetectorParameters hCalRingParameters;
        const gear::CalorimeterParameters &hCalRingInputParameters = marlin::Global::GEAR->getHcalRingParameters();
        SetDefaultSubDetectorParameters(hCalRingInputParameters, hCalRingParameters);
        geometryParameters.m_additionalSubDetectors["HCalRing"] = hCalRingParameters;
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(WARNING) << "Marlin pandora geometry creator: " << exception.what() << std::endl;
    }

    try
    {
        PandoraApi::Geometry::Parameters::SubDetectorParameters lCalParameters;
        const gear::CalorimeterParameters &lCalInputParameters = marlin::Global::GEAR->getLcalParameters();
        SetDefaultSubDetectorParameters(lCalInputParameters, lCalParameters);
        geometryParameters.m_additionalSubDetectors["LCal"] = lCalParameters;
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(WARNING) << "Marlin pandora geometry creator: " << exception.what() << std::endl;
    }

    try
    {
        PandoraApi::Geometry::Parameters::SubDetectorParameters lHCalParameters;
        const gear::CalorimeterParameters &lHCalInputParameters = marlin::Global::GEAR->getLHcalParameters();
        SetDefaultSubDetectorParameters(lHCalInputParameters, lHCalParameters);
        geometryParameters.m_additionalSubDetectors["LHCal"] = lHCalParameters;
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(WARNING) << "Marlin pandora geometry creator: " << exception.what() << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::SetILDSpecificGeometry(PandoraApi::GeometryParameters &geometryParameters) const
{
    // Non-default values (and those missing from GEAR parameters file)...
    const gear::CalorimeterParameters &hCalBarrelParameters = marlin::Global::GEAR->getHcalBarrelParameters();
    geometryParameters.m_hCalBarrelParameters.m_outerPhiCoordinate = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_phi0");
    geometryParameters.m_hCalBarrelParameters.m_outerSymmetryOrder = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_order");

    geometryParameters.m_eCalEndCapParameters.m_innerSymmetryOrder = m_settings.m_eCalEndCapInnerSymmetryOrder;
    geometryParameters.m_eCalEndCapParameters.m_innerPhiCoordinate = m_settings.m_eCalEndCapInnerPhiCoordinate;
    geometryParameters.m_eCalEndCapParameters.m_outerSymmetryOrder = m_settings.m_eCalEndCapOuterSymmetryOrder;
    geometryParameters.m_eCalEndCapParameters.m_outerPhiCoordinate = m_settings.m_eCalEndCapOuterPhiCoordinate;

    geometryParameters.m_hCalEndCapParameters.m_innerSymmetryOrder = m_settings.m_hCalEndCapInnerSymmetryOrder;
    geometryParameters.m_hCalEndCapParameters.m_innerPhiCoordinate = m_settings.m_hCalEndCapInnerPhiCoordinate;
    geometryParameters.m_hCalEndCapParameters.m_outerSymmetryOrder = m_settings.m_hCalEndCapOuterSymmetryOrder;
    geometryParameters.m_hCalEndCapParameters.m_outerPhiCoordinate = m_settings.m_hCalEndCapOuterPhiCoordinate;

    // Gaps in detector active material
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalBarrelBoxGaps(geometryParameters));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalEndCapBoxGaps(geometryParameters));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalBarrelConcentricGaps(geometryParameters));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::CreateHCalBarrelBoxGaps(PandoraApi::GeometryParameters &geometryParameters) const
{
    const std::string detectorName(marlin::Global::GEAR->getDetectorName());

    const gear::CalorimeterParameters &hCalBarrelParameters = marlin::Global::GEAR->getHcalBarrelParameters();
    const unsigned int innerSymmetryOrder(hCalBarrelParameters.getSymmetryOrder());
    const unsigned int outerSymmetryOrder(hCalBarrelParameters.getIntVal("Hcal_outer_polygon_order"));

    if ((0 == innerSymmetryOrder) || (2 != outerSymmetryOrder / innerSymmetryOrder))
    {
        streamlog_out(ERROR) << " Detector " << detectorName << " doesn't conform to expected ILD-specific geometry" << std::endl;
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    const float innerRadius(hCalBarrelParameters.getExtent()[0]);
    const float outerRadius(hCalBarrelParameters.getExtent()[1]);
    const float outerZ(hCalBarrelParameters.getExtent()[3]);
    const float phi0(hCalBarrelParameters.getPhi0());

    const float staveGap(hCalBarrelParameters.getDoubleVal("Hcal_stave_gaps"));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(innerSymmetryOrder, phi0, innerRadius, outerRadius,
        -outerZ, outerZ, staveGap, geometryParameters));

    static const float pi(std::acos(-1.));
    const float outerPseudoPhi0(pi / static_cast<float>(innerSymmetryOrder));
    const float cosOuterPseudoPhi0(std::cos(outerPseudoPhi0));

    if ((0 == outerPseudoPhi0) || (0.f == cosOuterPseudoPhi0))
    {
        streamlog_out(ERROR) << " Detector " << detectorName << " doesn't conform to expected ILD-specific geometry" << std::endl;
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    const float middleStaveGap(hCalBarrelParameters.getDoubleVal("Hcal_middle_stave_gaps"));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(innerSymmetryOrder, outerPseudoPhi0,
        innerRadius / cosOuterPseudoPhi0, outerRadius, -outerZ, outerZ, middleStaveGap, geometryParameters));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::CreateHCalEndCapBoxGaps(PandoraApi::GeometryParameters &geometryParameters) const
{
    const gear::CalorimeterParameters &hCalEndCapParameters = marlin::Global::GEAR->getHcalEndcapParameters();

    const float staveGap(hCalEndCapParameters.getDoubleVal("Hcal_stave_gaps"));
    const float innerRadius(hCalEndCapParameters.getExtent()[0]);
    const float outerRadius(hCalEndCapParameters.getExtent()[1]);
    const float innerZ(hCalEndCapParameters.getExtent()[2]);
    const float outerZ(hCalEndCapParameters.getExtent()[3]);

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(m_settings.m_hCalEndCapInnerSymmetryOrder,
        m_settings.m_hCalEndCapInnerPhiCoordinate, innerRadius, outerRadius, innerZ, outerZ, staveGap, geometryParameters,
        pandora::CartesianVector(-innerRadius, 0, 0)));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(m_settings.m_hCalEndCapInnerSymmetryOrder,
        m_settings.m_hCalEndCapInnerPhiCoordinate, innerRadius, outerRadius, -outerZ, -innerZ, staveGap, geometryParameters,
        pandora::CartesianVector(innerRadius, 0, 0)));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::CreateHCalBarrelConcentricGaps(PandoraApi::GeometryParameters &geometryParameters) const
{
    const gear::CalorimeterParameters &hCalBarrelParameters = marlin::Global::GEAR->getHcalBarrelParameters();
    const float gapWidth(hCalBarrelParameters.getDoubleVal("Hcal_stave_gaps"));

    PandoraApi::ConcentricGap::Parameters gapParameters;

    gapParameters.m_minZCoordinate = -0.5f * gapWidth;
    gapParameters.m_maxZCoordinate =  0.5f * gapWidth;
    gapParameters.m_innerRCoordinate = hCalBarrelParameters.getExtent()[0];
    gapParameters.m_innerPhiCoordinate = hCalBarrelParameters.getPhi0();
    gapParameters.m_innerSymmetryOrder = hCalBarrelParameters.getSymmetryOrder();
    gapParameters.m_outerRCoordinate = hCalBarrelParameters.getExtent()[1];
    gapParameters.m_outerPhiCoordinate = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_phi0");
    gapParameters.m_outerSymmetryOrder = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_order");

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ConcentricGap::Create(*m_pPandora, gapParameters));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::CreateRegularBoxGaps(unsigned int symmetryOrder, float phi0, float innerRadius, float outerRadius,
    float minZ, float maxZ, float gapWidth, PandoraApi::GeometryParameters &geometryParameters, pandora::CartesianVector vertexOffset) const
{
    const pandora::CartesianVector basicGapVertex(pandora::CartesianVector(-0.5f * gapWidth, innerRadius, minZ) + vertexOffset);
    const pandora::CartesianVector basicSide1(gapWidth, 0, 0);
    const pandora::CartesianVector basicSide2(0, outerRadius - innerRadius, 0);
    const pandora::CartesianVector basicSide3(0, 0, maxZ - minZ);

    for (unsigned int i = 0; i < symmetryOrder; ++i)
    {
        static const float pi(std::acos(-1.));

        const float phi = phi0 + (2. * pi * static_cast<float>(i) / static_cast<float>(symmetryOrder));
        const float sinPhi(std::sin(phi));
        const float cosPhi(std::cos(phi));

        PandoraApi::BoxGap::Parameters gapParameters;

        gapParameters.m_vertex = pandora::CartesianVector(cosPhi * basicGapVertex.GetX() + sinPhi * basicGapVertex.GetY(),
            -sinPhi * basicGapVertex.GetX() + cosPhi * basicGapVertex.GetY(), basicGapVertex.GetZ());
        gapParameters.m_side1 = pandora::CartesianVector(cosPhi * basicSide1.GetX() + sinPhi * basicSide1.GetY(),
            -sinPhi * basicSide1.GetX() + cosPhi * basicSide1.GetY(), basicSide1.GetZ());
        gapParameters.m_side2 = pandora::CartesianVector(cosPhi * basicSide2.GetX() + sinPhi * basicSide2.GetY(),
            -sinPhi * basicSide2.GetX() + cosPhi * basicSide2.GetY(), basicSide2.GetZ());
        gapParameters.m_side3 = pandora::CartesianVector(cosPhi * basicSide3.GetX() + sinPhi * basicSide3.GetY(),
            -sinPhi * basicSide3.GetX() + cosPhi * basicSide3.GetY(), basicSide3.GetZ());

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::BoxGap::Create(*m_pPandora, gapParameters));
    }

    return pandora::STATUS_CODE_SUCCESS;
}
