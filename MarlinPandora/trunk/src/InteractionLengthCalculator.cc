/**
 *  @file   PandoraPFANew/src/InteractionLengthCalculator.cc
 * 
 *  @brief  Implementation of the calo hit creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "gear/CalorimeterParameters.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"

#include "InteractionLengthCalculator.h"

#include <limits>

bool InteractionLengthCalculator::m_instanceFlag = false;
InteractionLengthCalculator* InteractionLengthCalculator::m_pInteractionLengthCalculator = NULL;

float InteractionLengthCalculator::Settings::m_avgIntLengthTracker = 0.f;
float InteractionLengthCalculator::Settings::m_avgIntLengthCoil = 0.f;
float InteractionLengthCalculator::Settings::m_avgIntLengthECalBarrel = 0.f;
float InteractionLengthCalculator::Settings::m_avgIntLengthHCalBarrel = 0.f;
float InteractionLengthCalculator::Settings::m_avgIntLengthECalEndCap = 0.f;
float InteractionLengthCalculator::Settings::m_avgIntLengthHCalEndCap = 0.f;
float InteractionLengthCalculator::Settings::m_avgIntLengthMuonBarrel = 0.f;
float InteractionLengthCalculator::Settings::m_avgIntLengthMuonEndCap = 0.f;

//------------------------------------------------------------------------------------------------------------------------------------------

InteractionLengthCalculator *InteractionLengthCalculator::GetInstance()
{
    if (!m_instanceFlag)
    {
        m_pInteractionLengthCalculator = new InteractionLengthCalculator();
        m_instanceFlag = true;
    }

    return m_pInteractionLengthCalculator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

InteractionLengthCalculator::InteractionLengthCalculator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

InteractionLengthCalculator::~InteractionLengthCalculator()
{
    m_instanceFlag = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float InteractionLengthCalculator::GetNInteractionLengthsFromIP(const EVENT::CalorimeterHit *const pCaloHit)
{
    try
    {
        static bool initialized = false;

        // Coordinates of the sub-detectors in one quadrant
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

        if (!initialized)
        {
            const PandoraApi::Geometry::Parameters geometryParameters;

            const gear::TPCParameters &tpcParameters = marlin::Global::GEAR->getTPCParameters();
            const gear::PadRowLayout2D &tpcPadLayout = tpcParameters.getPadLayout();
            rMinTracker = tpcPadLayout.getPlaneExtent()[0];
            rMaxTracker = tpcPadLayout.getPlaneExtent()[1];
            zMinTracker = 0;
            zMaxTracker = tpcParameters.getMaxDriftLength();

            const gear::GearParameters &coilParameters = marlin::Global::GEAR->getGearParameters("CoilParameters");
            rMinCoil = coilParameters.getDoubleVal("Coil_cryostat_inner_radius");
            rMaxCoil = coilParameters.getDoubleVal("Coil_cryostat_outer_radius");
            zMinCoil = 0.f;
            zMaxCoil = coilParameters.getDoubleVal("Coil_cryostat_half_z");

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

            initialized = true;
        }

        const float *pPosition = pCaloHit->getPosition();
        pandora::CartesianVector positionVector(pPosition[0], pPosition[1], pPosition[2]);

        float radius, phi, z;
        positionVector.GetCylindricalCoordinates(radius, phi, z);
        positionVector.SetValues(radius, 0.f, std::fabs(z));

        float nInteractionLengths(0.f);

        nInteractionLengths += ComputePathLengthFromIPInRectangle(positionVector, rMinECalBarrel, zMinECalBarrel, rMaxECalBarrel, zMaxECalBarrel) * Settings::m_avgIntLengthECalBarrel;
        nInteractionLengths += ComputePathLengthFromIPInRectangle(positionVector, rMinHCalBarrel, zMinHCalBarrel, rMaxHCalBarrel, zMaxHCalBarrel) * Settings::m_avgIntLengthHCalBarrel;

        nInteractionLengths += ComputePathLengthFromIPInRectangle(positionVector, rMinCoil, zMinCoil, rMaxCoil, zMaxCoil) * Settings::m_avgIntLengthCoil;

        nInteractionLengths += ComputePathLengthFromIPInRectangle(positionVector, rMinECalEndCap, zMinECalEndCap, rMaxECalEndCap, zMaxECalEndCap) * Settings::m_avgIntLengthECalEndCap;
        nInteractionLengths += ComputePathLengthFromIPInRectangle(positionVector, rMinHCalEndCap, zMinHCalEndCap, rMaxHCalEndCap, zMaxHCalEndCap) * Settings::m_avgIntLengthHCalEndCap;

        nInteractionLengths += ComputePathLengthFromIPInRectangle(positionVector, rMinMuonBarrel, zMinMuonBarrel, rMaxMuonBarrel, zMaxMuonBarrel) * Settings::m_avgIntLengthMuonBarrel;
        nInteractionLengths += ComputePathLengthFromIPInRectangle(positionVector, rMinMuonEndCap, zMinMuonEndCap, rMaxMuonEndCap, zMaxMuonEndCap) * Settings::m_avgIntLengthMuonEndCap;

        return nInteractionLengths;
    }
    catch (gear::Exception &exception)
    {
        streamlog_out(ERROR) << "InteractionLengthCalculator: failed to extract gear geometry information" << std::endl;
        throw exception;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float InteractionLengthCalculator::ComputePathLengthFromIPInRectangle(const pandora::CartesianVector &position, float rMin, float zMin,
    float rMax, float zMax)
{
    // compute cuts with rectangle borders
    float phi, radius, z;
    position.GetCylindricalCoordinates(radius, phi, z);

    bool valid[4];
    float xInt[4], zInt[4];

    valid[0] = IntersectLines2D(0.f, 0.f, radius, z, rMin, zMin, rMin, zMax, xInt[0], zInt[0]); // first edge of rectangle at rMin
    valid[1] = IntersectLines2D(0.f, 0.f, radius, z, rMin, zMin, rMax, zMin, xInt[1], zInt[1]); // first edge of rectangle at zMin
    valid[2] = IntersectLines2D(0.f, 0.f, radius, z, rMax, zMax, rMax, zMin, xInt[2], zInt[2]); // first edge of rectangle at rMax
    valid[3] = IntersectLines2D(0.f, 0.f, radius, z, rMax, zMax, rMin, zMax, xInt[3], zInt[3]); // first edge of rectangle at zMax

    int indexFirstPoint = -1;
    int indexSecondPoint = -1;

    for (int i = 0; i < 4; ++i)
    {
        if (valid[i])
        {
            if (indexFirstPoint == -1)
            {
                indexFirstPoint = i;
            }
            else if (indexSecondPoint == -1)
            {
                indexSecondPoint = i;
            }
            else
            {
                std::cout << "ERROR calohitcreator problem at computing path length from IP to calohit in rectangle" << std::endl;
            }
        }
    }

    if (indexFirstPoint == -1)
        return 0.f;

    pandora::CartesianVector intersectionA( xInt[indexFirstPoint], 0.f, zInt[indexFirstPoint] );

    if (indexSecondPoint == -1)
    {
        const float length(pandora::CartesianVector(position - intersectionA).GetMagnitude());
        return length;
    }

    pandora::CartesianVector intersectionB(xInt[indexSecondPoint], 0.f, zInt[indexSecondPoint]);

    const float length(pandora::CartesianVector(intersectionA - intersectionB).GetMagnitude());
    return length;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool InteractionLengthCalculator::IntersectLines2D(float lineAXStart, float lineAYStart, float lineAXEnd, float lineAYEnd,
    float lineBXStart, float lineBYStart, float lineBXEnd, float lineBYEnd, float &xIntersect, float &yIntersect)
{
    // Slopes of the two lines, take max float value instead of infinity
    float k0(std::numeric_limits<float>::max()), k1(std::numeric_limits<float>::max());

    bool parallelToY_A = false;
    bool parallelToY_B = false;

    if (std::fabs(lineAXEnd-lineAXStart) > std::numeric_limits<float>::epsilon())
    {
        k0 = (lineAYEnd - lineAYStart) / (lineAXEnd - lineAXStart);
    }
    else
    {
        parallelToY_A = true;
    }

    if (std::fabs(lineBXEnd - lineBXStart) > std::numeric_limits<float>::epsilon())
    {
        k1 = (lineBYEnd - lineBYStart) / (lineBXEnd - lineBXStart);
    }
    else
    {
        parallelToY_B = true;
    }

    if (parallelToY_A && parallelToY_B)
    {
        xIntersect = 0.f;
        yIntersect = 0.f;

        return false;
    }

    if (parallelToY_A)
    {
        xIntersect = lineAXStart;
        yIntersect = lineBXStart + k1 * (xIntersect - lineBXStart);
    }
    else if (parallelToY_B)
    {
        xIntersect = lineBXStart;
        yIntersect = lineAXStart + k0 * (xIntersect - lineAXStart);
    }
    else
    {
        const float b0 = -1;
        const float b1 = -1;

        const float c0 = (lineAYStart - k0 * lineAXStart);
        const float c1 = (lineBYStart - k1 * lineBXStart);

        const float determinant(k0 * b1 - k1 * b0);

        if (0.f == determinant)
        {
            std::cout << "ERROR zero determinant in interaction length calculator" << std::endl;
            return false;
        }

        // use Kramers rule to compute xi and yi
        xIntersect=((b0 * c1 - b1 * c0) / determinant);
        yIntersect=((k1 * c0 - k0 * c1) / determinant);
    }

    // check if intersections are within end-points of lines
    // check x coordinate, line A
    if (!((xIntersect >= lineAXStart && xIntersect <= lineAXEnd) || (xIntersect >= lineAXEnd && xIntersect <= lineAXStart)))
        return false;

    // check x coordinate, line B
    if (!((xIntersect >= lineBXStart && xIntersect <= lineBXEnd) || (xIntersect >= lineBXEnd && xIntersect <= lineBXStart)))
        return false;

    // check y coordinate, line A
    if (!((yIntersect >= lineAYStart && yIntersect <= lineAYEnd) || (yIntersect >= lineAYEnd && yIntersect <= lineAYStart)))
        return false;

    // check y coordinate, line B
    if (!((yIntersect >= lineBYStart && yIntersect <= lineBYEnd) || (yIntersect >= lineBYEnd && yIntersect <= lineBYStart)))
        return false;

    return true;
} 
