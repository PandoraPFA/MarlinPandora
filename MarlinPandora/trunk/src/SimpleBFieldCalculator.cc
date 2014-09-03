/**
 *  @file   MarlinPandora/src/SimpleBFieldCalculator.cc
 * 
 *  @brief  Implementation of the simple bfield calculator class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "SimpleBFieldCalculator.h"

// TODO static make parameters configurable via Marlin processor
SimpleBFieldCalculator::SimpleBFieldCalculator() :
    m_innerBField(4.f),
    m_muonBarrelBField(-1.5f),
    m_muonEndCapBField(0.01f),
    m_muonEndCapInnerZ(std::numeric_limits<float>::max()),
    m_coilMidPointR(std::numeric_limits<float>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SimpleBFieldCalculator::~SimpleBFieldCalculator()
{
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

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode SimpleBFieldCalculator::Initialize()
{
// TODO read detector name from settings
    try
    {
        m_muonEndCapInnerZ = this->GetPandora().GetGeometry()->GetSubDetector("MuonEndCap").GetInnerZCoordinate();
        m_coilMidPointR = (0.5f * (this->GetPandora().GetGeometry()->GetSubDetector("Coil").GetInnerRCoordinate() +
            this->GetPandora().GetGeometry()->GetSubDetector("Coil").GetOuterRCoordinate()));
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        std::cout << "SimpleBFieldCalculator: Unable to extract Muon EndCap and Coil geometry." << std::endl;
        return statusCodeException.GetStatusCode();
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode SimpleBFieldCalculator::ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/)
{
    return pandora::STATUS_CODE_SUCCESS;
}
