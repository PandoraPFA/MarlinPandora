/**
 *  @file   MarlinPandora/include/NonLinearityCorrection.h
 * 
 *  @brief  Header file for the non linearity correction class.
 * 
 *  $Log: $
 */

#ifndef NON_LINEARITY_CORRECTION_H
#define NON_LINEARITY_CORRECTION_H 1

#include "Plugins/EnergyCorrectionsPlugin.h"

/**
 *   @brief  Correct cluster energy to account for non-linearities in calibration
 */
class NonLinearityCorrection : public pandora::EnergyCorrectionPlugin
{
public:
    typedef std::vector<float> FloatVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        FloatVector     m_inputEnergyCorrectionPoints;      ///< The input energy points for energy correction
        FloatVector     m_outputEnergyCorrectionPoints;     ///< The output energy points for energy correction
    };

    /**
     *  @brief  Default constructor
     * 
     *  @param  settings the settings
     */
    NonLinearityCorrection(const Settings &settings);

    pandora::StatusCode MakeEnergyCorrections(const pandora::Cluster *const pCluster, float &correctedEnergy) const;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    FloatVector         m_inputEnergyCorrectionPoints;      ///< The input energy points for energy correction
    FloatVector         m_energyCorrections;                ///< The energy correction factors
};

#endif // #ifndef NON_LINEARITY_CORRECTION_H
