/**
 *  @file   MarlinPandora/include/SimpleBFieldCalculator.h
 * 
 *  @brief  Header file for the simple bfield calculator class.
 * 
 *  $Log: $
 */

#ifndef SIMPLE_BFIELD_CALCULATOR_H
#define SIMPLE_BFIELD_CALCULATOR_H 1

#include "Plugins/BFieldPlugin.h"

/**
 *  @brief  SimpleBFieldCalculator class
 */
class SimpleBFieldCalculator : public pandora::BFieldPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    SimpleBFieldCalculator();

    /**
     *  @brief  Destructor
     */
    ~SimpleBFieldCalculator();

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

#endif // #ifndef SIMPLE_BFIELD_CALCULATOR_H
