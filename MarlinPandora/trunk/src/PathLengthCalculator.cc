/**
 *  @file   MarlinPandora/src/PathLengthCalculator.cc
 * 
 *  @brief  Implementation of the calo hit creator class.
 * 
 *  $Log: $
 */

#include "marlin/Processor.h"

#include "PathLengthCalculator.h"

#include <limits>

PathLengthCalculator *PathLengthCalculator::GetInstance()
{
    if (!m_instanceFlag)
    {
        m_pPathLengthCalculator = new PathLengthCalculator();
        m_instanceFlag = true;
    }

    return m_pPathLengthCalculator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PathLengthCalculator::PathLengthCalculator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PathLengthCalculator::~PathLengthCalculator()
{
    m_instanceFlag = false;
    m_isGeometryInitialized = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PathLengthCalculator::GetPathLengths(const EVENT::CalorimeterHit *const pCaloHit, float &nRadiationLengthsFromIp,
    float &nInteractionLengthsFromIp)
{
    if (!m_isGeometryInitialized)
        PathLengthCalculator::InitializeGeometry();

    const pandora::CartesianVector positionVector(pCaloHit->getPosition()[0], pCaloHit->getPosition()[1], pCaloHit->getPosition()[2]);
    const float lineLength(positionVector.GetMagnitude());

    const float eCalBarrelPathLength(lineLength * GetFractionInSubDetector(positionVector, m_eCalBarrelParameters));
    const float hCalBarrelPathLength(lineLength * GetFractionInSubDetector(positionVector, m_hCalBarrelParameters));
    const float muonBarrelPathLength(lineLength * GetFractionInSubDetector(positionVector, m_muonBarrelParameters));
    const float eCalEndCapPathLength(lineLength * GetFractionInSubDetector(positionVector, m_eCalEndCapParameters));
    const float hCalEndCapPathLength(lineLength * GetFractionInSubDetector(positionVector, m_hCalEndCapParameters));
    const float muonEndCapPathLength(lineLength * GetFractionInSubDetector(positionVector, m_muonEndCapParameters));
    const float coilPathLength(lineLength * GetFractionInSubDetector(positionVector, m_coilParameters));

    nRadiationLengthsFromIp = (eCalBarrelPathLength * Settings::m_avgRadLengthECalBarrel) +
        (hCalBarrelPathLength * Settings::m_avgRadLengthHCalBarrel) +
        (muonBarrelPathLength * Settings::m_avgRadLengthMuonBarrel) +
        (eCalEndCapPathLength * Settings::m_avgRadLengthECalEndCap) +
        (hCalEndCapPathLength * Settings::m_avgRadLengthHCalEndCap) +
        (muonEndCapPathLength * Settings::m_avgRadLengthMuonEndCap) +
        (coilPathLength * Settings::m_avgRadLengthCoil);

    nInteractionLengthsFromIp = (eCalBarrelPathLength * Settings::m_avgIntLengthECalBarrel) +
        (hCalBarrelPathLength * Settings::m_avgIntLengthHCalBarrel) +
        (muonBarrelPathLength * Settings::m_avgIntLengthMuonBarrel) +
        (eCalEndCapPathLength * Settings::m_avgIntLengthECalEndCap) +
        (hCalEndCapPathLength * Settings::m_avgIntLengthHCalEndCap) +
        (muonEndCapPathLength * Settings::m_avgIntLengthMuonEndCap) +
        (coilPathLength * Settings::m_avgIntLengthCoil);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PathLengthCalculator::GetFractionInSubDetector(const pandora::CartesianVector &positionVector, const SubDetectorParameters &subDetectorParameters)
{
    const float rCoordinate(std::sqrt(positionVector.GetX() * positionVector.GetX() + positionVector.GetY() * positionVector.GetY()));
    const float zCoordinate(std::fabs(positionVector.GetZ()));

    if ((rCoordinate < subDetectorParameters.GetInnerRCoordinate()) && (zCoordinate < subDetectorParameters.GetInnerZCoordinate()))
        return 0.f;

    // Use r over z ratios to determine where line enters/exits the subdetector
    const float rOverZ((zCoordinate == 0.f) ? std::numeric_limits<float>::max() : rCoordinate / zCoordinate);

    const float innerROverZ((subDetectorParameters.GetInnerZCoordinate() == 0.f) ? std::numeric_limits<float>::max() :
        subDetectorParameters.GetInnerRCoordinate() / subDetectorParameters.GetInnerZCoordinate());

    const float outerROverZ((subDetectorParameters.GetOuterZCoordinate() == 0.f) ? std::numeric_limits<float>::max() :
        subDetectorParameters.GetOuterRCoordinate() / subDetectorParameters.GetOuterZCoordinate());

    // Find point at which line enters subdetector
    float innerFraction(0.f);

    if (rOverZ <= innerROverZ)
    {
        innerFraction = GetLengthFraction(positionVector, subDetectorParameters.GetInnerRCoordinate(), subDetectorParameters.GetInnerRNormalVectors());
    }
    else
    {
        innerFraction = GetLengthFraction(positionVector, subDetectorParameters.GetInnerZCoordinate(), subDetectorParameters.GetZNormalVectors());
    }

    if ((innerFraction <= 0.f) || (innerFraction >= 1.f))
        return 0.f;

    // Point at which line exits subdetector
    float outerFraction(0.f);

    if (rOverZ >= outerROverZ)
    {
        outerFraction = GetLengthFraction(positionVector, subDetectorParameters.GetOuterRCoordinate(), subDetectorParameters.GetOuterRNormalVectors());
    }
    else
    {
        outerFraction = GetLengthFraction(positionVector, subDetectorParameters.GetOuterZCoordinate(), subDetectorParameters.GetZNormalVectors());
    }

    if (outerFraction < innerFraction)
        return 0.f;

    return (std::min(1.f, outerFraction) - innerFraction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PathLengthCalculator::GetLengthFraction(const pandora::CartesianVector &positionVector, const float closestDistanceToIp,
    const NormalVectorList &normalVectorList)
{
    // Deal with cylindrical case
    if (normalVectorList.empty())
    {
        const float radius(std::sqrt(positionVector.GetX() * positionVector.GetX() + positionVector.GetY() * positionVector.GetY()));
        return closestDistanceToIp / radius;
    }

    // Deal with regular polygon case
    float maxDotProduct(0.f);

    for (NormalVectorList::const_iterator iter = normalVectorList.begin(), iterEnd = normalVectorList.end(); iter != iterEnd; ++iter)
    {
        const float dotProduct(positionVector.GetDotProduct(*iter));

        if (dotProduct > maxDotProduct)
            maxDotProduct = dotProduct;
    }

    if (maxDotProduct <= 0.f)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    return closestDistanceToIp / maxDotProduct;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PathLengthCalculator::InitializeGeometry()
{
    try
    {
        m_eCalBarrelParameters.Initialize(pandora::GeometryHelper::GetECalBarrelParameters());
        m_hCalBarrelParameters.Initialize(pandora::GeometryHelper::GetHCalBarrelParameters());
        m_muonBarrelParameters.Initialize(pandora::GeometryHelper::GetMuonBarrelParameters());
        m_eCalEndCapParameters.Initialize(pandora::GeometryHelper::GetECalEndCapParameters());
        m_hCalEndCapParameters.Initialize(pandora::GeometryHelper::GetHCalEndCapParameters());
        m_muonEndCapParameters.Initialize(pandora::GeometryHelper::GetMuonEndCapParameters());
        m_coilParameters.Initialize_Cylinder(pandora::GeometryHelper::GetCoilInnerRadius(), pandora::GeometryHelper::GetCoilOuterRadius(),
            0.f, pandora::GeometryHelper::GetCoilZExtent());

        m_isGeometryInitialized = true;
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        streamlog_out(ERROR) << "PathLengthCalculator: Failed to initialize geometry: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void PathLengthCalculator::SubDetectorParameters::Initialize(const pandora::GeometryHelper::SubDetectorParameters &subDetectorParameters)
{
    m_innerRCoordinate = subDetectorParameters.GetInnerRCoordinate();
    m_outerRCoordinate = subDetectorParameters.GetOuterRCoordinate();
    m_innerZCoordinate = subDetectorParameters.GetInnerZCoordinate();
    m_outerZCoordinate = subDetectorParameters.GetOuterZCoordinate();

    this->GetPolygonNormalVectors(subDetectorParameters.GetInnerSymmetryOrder(), subDetectorParameters.GetInnerPhiCoordinate(), m_innerRNormalVectors);
    this->GetPolygonNormalVectors(subDetectorParameters.GetOuterSymmetryOrder(), subDetectorParameters.GetOuterPhiCoordinate(), m_outerRNormalVectors);
    m_zNormalVectors.push_back(pandora::CartesianVector(0, 0, 1));
    m_zNormalVectors.push_back(pandora::CartesianVector(0, 0, -1));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PathLengthCalculator::SubDetectorParameters::Initialize_Cylinder(const float innerRCoordinate, const float outerRCoordinate,
    const float innerZCoordinate, const float outerZCoordinate)
{
    m_innerRCoordinate = innerRCoordinate;
    m_outerRCoordinate = outerRCoordinate;
    m_innerZCoordinate = innerZCoordinate;
    m_outerZCoordinate = outerZCoordinate;

    m_zNormalVectors.push_back(pandora::CartesianVector(0, 0, 1));
    m_zNormalVectors.push_back(pandora::CartesianVector(0, 0, -1));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PathLengthCalculator::SubDetectorParameters::GetPolygonNormalVectors(const unsigned int symmetry, const float phi0, NormalVectorList &normalVectorList) const
{
    for (unsigned int iSymmetry = 0; iSymmetry < symmetry; ++iSymmetry)
    {
        static const float pi(std::acos(-1.));
        const float phi = phi0 + (2. * pi * static_cast<float>(iSymmetry) / static_cast<float>(symmetry));
        normalVectorList.push_back(pandora::CartesianVector(std::sin(phi), std::cos(phi), 0));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PathLengthCalculator::m_instanceFlag = false;
bool PathLengthCalculator::m_isGeometryInitialized = false;
PathLengthCalculator* PathLengthCalculator::m_pPathLengthCalculator = NULL;

PathLengthCalculator::SubDetectorParameters PathLengthCalculator::m_eCalBarrelParameters;
PathLengthCalculator::SubDetectorParameters PathLengthCalculator::m_hCalBarrelParameters;
PathLengthCalculator::SubDetectorParameters PathLengthCalculator::m_muonBarrelParameters;
PathLengthCalculator::SubDetectorParameters PathLengthCalculator::m_eCalEndCapParameters;
PathLengthCalculator::SubDetectorParameters PathLengthCalculator::m_hCalEndCapParameters;
PathLengthCalculator::SubDetectorParameters PathLengthCalculator::m_muonEndCapParameters;
PathLengthCalculator::SubDetectorParameters PathLengthCalculator::m_coilParameters;

float PathLengthCalculator::Settings::m_avgRadLengthCoil = 0.f;
float PathLengthCalculator::Settings::m_avgRadLengthECalBarrel = 0.f;
float PathLengthCalculator::Settings::m_avgRadLengthHCalBarrel = 0.f;
float PathLengthCalculator::Settings::m_avgRadLengthECalEndCap = 0.f;
float PathLengthCalculator::Settings::m_avgRadLengthHCalEndCap = 0.f;
float PathLengthCalculator::Settings::m_avgRadLengthMuonBarrel = 0.f;
float PathLengthCalculator::Settings::m_avgRadLengthMuonEndCap = 0.f;

float PathLengthCalculator::Settings::m_avgIntLengthCoil = 0.f;
float PathLengthCalculator::Settings::m_avgIntLengthECalBarrel = 0.f;
float PathLengthCalculator::Settings::m_avgIntLengthHCalBarrel = 0.f;
float PathLengthCalculator::Settings::m_avgIntLengthECalEndCap = 0.f;
float PathLengthCalculator::Settings::m_avgIntLengthHCalEndCap = 0.f;
float PathLengthCalculator::Settings::m_avgIntLengthMuonBarrel = 0.f;
float PathLengthCalculator::Settings::m_avgIntLengthMuonEndCap = 0.f;
