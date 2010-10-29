/**
 *  @file   PandoraPFANew/include/SimpleBFieldCalculator.h
 * 
 *  @brief  Header file for the simple bfield calculator class.
 * 
 *  $Log: $
 */

#ifndef SIMPLE_BFIELD_CALCULATOR_H
#define SIMPLE_BFIELD_CALCULATOR_H 1

#include "Utilities/BFieldCalculator.h"

/**
 *  @brief  SimpleBFieldCalculator class
 */
class SimpleBFieldCalculator : public pandora::BFieldCalculator
{
public:
    static float        m_innerBField;              ///< The bfield in the main tracker, ecal and hcal, units Tesla
    static float        m_muonBarrelBField;         ///< The bfield in the muon barrel, units Tesla
    static float        m_muonEndCapBField;         ///< The bfield in the muon endcap, units Tesla

private:
    void Initialize(const pandora::GeometryHelper *const pGeometryHelper);
    float GetBField(const pandora::CartesianVector &positionVector) const;

    float               m_muonEndCapInnerZ;         ///< The muon endcap inner z coordinate, units mm
    float               m_coilMidPointR;            ///< The r coordinate at the coil midpoint, units mm
};

#endif // #ifndef SIMPLE_BFIELD_CALCULATOR_H
