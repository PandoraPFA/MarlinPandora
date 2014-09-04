/**
 *  @file   MarlinPandora/src/BFieldPlugin.cc
 * 
 *  @brief  Implementation of the bfield plugin class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "BFieldPlugin.h"

BFieldPlugin::BFieldPlugin(const Settings &settings) :
    m_innerBField(settings.m_innerBField),
    m_muonBarrelBField(settings.m_muonBarrelBField),
    m_muonEndCapBField(settings.m_muonEndCapBField),
    m_muonEndCapInnerZ(std::numeric_limits<float>::max()),
    m_coilMidPointR(std::numeric_limits<float>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float BFieldPlugin::GetBField(const pandora::CartesianVector &positionVector) const
{
    if (std::fabs(positionVector.GetZ()) >= m_muonEndCapInnerZ)
        return m_muonEndCapBField;

    if (std::sqrt(positionVector.GetX() * positionVector.GetX() + positionVector.GetY() * positionVector.GetY()) >= m_coilMidPointR)
        return m_muonBarrelBField;

    return m_innerBField;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode BFieldPlugin::Initialize()
{
    try
    {
        m_muonEndCapInnerZ = this->GetPandora().GetGeometry()->GetSubDetector(pandora::MUON_ENDCAP).GetInnerZCoordinate();
        m_coilMidPointR = (0.5f * (this->GetPandora().GetGeometry()->GetSubDetector(pandora::COIL).GetInnerRCoordinate() +
            this->GetPandora().GetGeometry()->GetSubDetector(pandora::COIL).GetOuterRCoordinate()));
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        std::cout << "BFieldPlugin: Unable to extract Muon EndCap and Coil geometry." << std::endl;
        return statusCodeException.GetStatusCode();
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode BFieldPlugin::ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/)
{
    return pandora::STATUS_CODE_SUCCESS;
}
