/**
 *  @file   PandoraPFANew/include/GeometryCreator.h
 * 
 *  @brief  Header file for the geometry creator class.
 * 
 *  $Log: $
 */

#ifndef GEOMETRY_CREATOR_H
#define GEOMETRY_CREATOR_H 1

#include "gear/CalorimeterParameters.h"

#include "Api/PandoraApi.h"

/**
 *  @brief  GeometryCreator class
 */
class GeometryCreator
{
public:
    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        float           m_absorberRadiationLength;              ///< The absorber radiation length
        float           m_absorberInteractionLength;            ///< The absorber interaction length

        int             m_eCalEndCapInnerSymmetryOrder;         ///< ECal end cap inner symmetry order (missing from ILD00 gear file)
        float           m_eCalEndCapInnerPhiCoordinate;         ///< ECal end cap inner phi coordinate (missing from ILD00 gear file)
        int             m_hCalEndCapInnerSymmetryOrder;         ///< HCal end cap inner symmetry order (missing from ILD00 gear file)
        float           m_hCalEndCapInnerPhiCoordinate;         ///< HCal end cap inner phi coordinate (missing from ILD00 gear file)
    };

    /**
     *  @brief  Create geometry, insert user code here
     */
    StatusCode CreateGeometry() const;

    Settings                m_settings;         ///< The settings

private:
    /**
     *  @brief  Set sub detector parameters to their gear default values
     * 
     *  @param  inputParameters input parameters, from gear
     *  @param  subDetectorParameters the sub detector parameters
     */
    void SetDefaultSubDetectorParameters(const gear::CalorimeterParameters &inputParameters,
        PandoraApi::GeometryParameters::SubDetectorParameters &subDetectorParameters) const;

    /**
     *  @brief  Set additional sub detector parameters
     * 
     *  @param  geometryParameters the pandora geometry parameters
     */
    void SetAdditionalSubDetectorParameters(PandoraApi::GeometryParameters &geometryParameters) const;

    /**
     *  @brief  Set positions of gaps in ILD detector and add information missing from GEAR parameters file
     * 
     *  @param  geometryParameters the pandora geometry parameters
     */
    StatusCode SetILDSpecificGeometry(PandoraApi::GeometryParameters &geometryParameters) const;

    /**
     *  @brief  Specify positions of hcal barrel box gaps - ILD specific
     * 
     *  @param  geometryParameters the pandora geometry parameters
     */
    StatusCode CreateHCalBarrelBoxGaps(PandoraApi::GeometryParameters &geometryParameters) const;

    /**
     *  @brief  Specify positions of hcal barrel concentric polygon gaps - ILD specific
     * 
     *  @param  geometryParameters the pandora geometry parameters
     */
    StatusCode CreateHCalBarrelConcentricGaps(PandoraApi::GeometryParameters &geometryParameters) const;

    /**
     *  @brief  Create box gaps at regular positions on polygonal prism, oriented along main z axis - ILD specific
     * 
     *  @param  symmetryOrder the pandora geometry parameters
     *  @param  phi0
     *  @param  innerRadius
     *  @param  outerRadius
     *  @param  outerZ
     *  @param  gapWidth
     *  @param  geometryParameters the pandora geometry parameters
     */
    StatusCode CreateRegularBoxGaps(unsigned int symmetryOrder, float phi0, float innerRadius, float outerRadius, float outerZ,
        float gapWidth, PandoraApi::GeometryParameters &geometryParameters) const;
};

#endif // #ifndef GEOMETRY_CREATOR_H
