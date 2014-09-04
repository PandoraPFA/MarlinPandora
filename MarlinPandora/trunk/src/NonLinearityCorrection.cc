/**
 *  @file   MarlinPandora/src/NonLinearityCorrection.cc
 * 
 *  @brief  Implementation of the non linearity correction class.
 * 
 *  $Log: $
 */

#include "NonLinearityCorrection.h"

NonLinearityCorrection::NonLinearityCorrection(const Settings &settings) :
    m_inputEnergyCorrectionPoints(settings.m_inputEnergyCorrectionPoints)
{
    const unsigned int nEnergyBins(m_inputEnergyCorrectionPoints.size());

    if (nEnergyBins != settings.m_outputEnergyCorrectionPoints.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    for (unsigned int i = 0; i < nEnergyBins; ++i)
    {
        const float inputEnergy(m_inputEnergyCorrectionPoints.at(i));
        const float outputEnergy(settings.m_outputEnergyCorrectionPoints.at(i));

        if (std::fabs(inputEnergy) < std::numeric_limits<float>::epsilon())
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

        m_energyCorrections.push_back(outputEnergy / inputEnergy);
    }

    if (nEnergyBins != m_energyCorrections.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode NonLinearityCorrection::MakeEnergyCorrections(const pandora::Cluster *const pCluster, float &correctedEnergy) const
{
    const unsigned int nEnergyBins(m_energyCorrections.size());
    unsigned int index(nEnergyBins);

    for (unsigned int i = 0; i < nEnergyBins; ++i)
    {
        if (correctedEnergy < m_inputEnergyCorrectionPoints.at(i))
        {
            index = i;
            break;
        }
    }

    float correction(1.f);

    if ((0 == index) || (nEnergyBins == index))
    {
        correction = m_energyCorrections.at(std::min(index, nEnergyBins - 1));
    }
    else
    {
        const float lowCorrection(m_energyCorrections.at(index - 1)), highCorrection(m_energyCorrections.at(index));
        const float lowEnergy(m_inputEnergyCorrectionPoints.at(index - 1)), highEnergy(m_inputEnergyCorrectionPoints.at(index));
        correction = lowCorrection + (correctedEnergy - lowEnergy) * (highCorrection - lowCorrection) / (highEnergy - lowEnergy);
    }

    correctedEnergy *= correction;

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode NonLinearityCorrection::ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/)
{
    return pandora::STATUS_CODE_SUCCESS;
}
