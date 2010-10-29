/**
 *  @file   PandoraPFANew/src/SimpleBFieldCalculator.cc
 * 
 *  @brief  Implementation of the simple bfield calculator class.
 * 
 *  $Log: $
 */

#include "SimpleBFieldCalculator.h"

#include <cmath>

float SimpleBFieldCalculator::m_innerBField = 4.f;
float SimpleBFieldCalculator::m_muonBarrelBField = 1.5f;
float SimpleBFieldCalculator::m_muonEndCapBField = 0.01f;

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleBFieldCalculator::Initialize(const pandora::GeometryHelper *const pGeometryHelper)
{
    m_muonEndCapInnerZ = pGeometryHelper->GetMuonEndCapParameters().GetInnerZCoordinate();
    m_coilMidPointR = (0.5f * (pGeometryHelper->GetCoilInnerRadius() + pGeometryHelper->GetCoilOuterRadius()));
};

//------------------------------------------------------------------------------------------------------------------------------------------

float SimpleBFieldCalculator::GetBField(const pandora::CartesianVector &positionVector) const
{
    if (std::fabs(positionVector.GetZ()) >= m_muonEndCapInnerZ)
        return m_muonEndCapBField;

    if (std::sqrt(positionVector.GetX() * positionVector.GetX() + positionVector.GetY() * positionVector.GetY()) >= m_coilMidPointR)
        return m_muonBarrelBField;

    return m_innerBField;
};
