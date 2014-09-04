/**
 *  @file   MarlinPandora/include/BFieldPlugin.h
 * 
 *  @brief  Header file for the bfield plugin class.
 * 
 *  $Log: $
 */

#ifndef BFIELD_PLUGIN_H
#define BFIELD_PLUGIN_H 1

#include "Plugins/BFieldPlugin.h"

/**
 *  @brief  BFieldPlugin class
 */
class BFieldPlugin : public pandora::BFieldPlugin
{
public:
    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        float   m_innerBField;          ///< The bfield in the main tracker, ecal and hcal, units Tesla
        float   m_muonBarrelBField;     ///< The bfield in the muon barrel, units Tesla
        float   m_muonEndCapBField;     ///< The bfield in the muon endcap, units Tesla
    };

    /**
     *  @brief  Default constructor
     * 
     *  @param  settings the settings
     */
    BFieldPlugin(const Settings &settings);

    float GetBField(const pandora::CartesianVector &positionVector) const;

private:
    pandora::StatusCode Initialize();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float   m_innerBField;              ///< The bfield in the main tracker, ecal and hcal, units Tesla
    float   m_muonBarrelBField;         ///< The bfield in the muon barrel, units Tesla
    float   m_muonEndCapBField;         ///< The bfield in the muon endcap, units Tesla

    float   m_muonEndCapInnerZ;         ///< The muon endcap inner z coordinate, units mm
    float   m_coilMidPointR;            ///< The r coordinate at the coil midpoint, units mm
};

#endif // #ifndef BFIELD_PLUGIN_H
