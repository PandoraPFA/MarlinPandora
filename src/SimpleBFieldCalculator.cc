/**
 *  @file   MarlinPandora/src/SimpleBFieldCalculator.cc
 * 
 *  @brief  Implementation of the simple bfield calculator class.
 * 
 *  $Log: $
 */

#include "Helpers/GeometryHelper.h"

#include "SimpleBFieldCalculator.h"

#include <cmath>

float SimpleBFieldCalculator::m_innerBField = 4.f;
float SimpleBFieldCalculator::m_muonBarrelBField = -1.5f;
float SimpleBFieldCalculator::m_muonEndCapBField = 0.01f;

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleBFieldCalculator::InitializeGeometry()
{
    m_muonEndCapInnerZ = pandora::GeometryHelper::GetMuonEndCapParameters().GetInnerZCoordinate();
    m_coilMidPointR = (0.5f * (pandora::GeometryHelper::GetCoilInnerRadius() + pandora::GeometryHelper::GetCoilOuterRadius()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float SimpleBFieldCalculator::GetBField(const pandora::CartesianVector &positionVector) const
{
    if (std::fabs(positionVector.GetZ()) >= m_muonEndCapInnerZ)
        return m_muonEndCapBField;

    if (std::sqrt(positionVector.GetX() * positionVector.GetX() + positionVector.GetY() * positionVector.GetY()) >= m_coilMidPointR)
        return m_muonBarrelBField;

    return m_innerBField;
}
