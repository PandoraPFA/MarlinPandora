/**
 *  @file   MarlinPandora/include/GeometryCreator.h
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

        int             m_eCalEndCapInnerSymmetryOrder;         ///< ECal end cap inner symmetry order (missing from ILD gear files)
        float           m_eCalEndCapInnerPhiCoordinate;         ///< ECal end cap inner phi coordinate (missing from ILD gear files)
        int             m_eCalEndCapOuterSymmetryOrder;         ///< ECal end cap outer symmetry order (missing from ILD gear files)
        float           m_eCalEndCapOuterPhiCoordinate;         ///< ECal end cap outer phi coordinate (missing from ILD gear files)

        int             m_hCalEndCapInnerSymmetryOrder;         ///< HCal end cap inner symmetry order (missing from ILD gear files)
        float           m_hCalEndCapInnerPhiCoordinate;         ///< HCal end cap inner phi coordinate (missing from ILD gear files)
        int             m_hCalEndCapOuterSymmetryOrder;         ///< HCal end cap outer symmetry order (missing from ILD gear files)
        float           m_hCalEndCapOuterPhiCoordinate;         ///< HCal end cap outer phi coordinate (missing from ILD gear files)
        
        int             m_hCalRingInnerSymmetryOrder;           ///< HCal ring inner symmetry order (missing from ILD gear files)
        float           m_hCalRingInnerPhiCoordinate;           ///< HCal ring inner phi coordinate (missing from ILD gear files)
        int             m_hCalRingOuterSymmetryOrder;           ///< HCal ring outer symmetry order (missing from ILD gear files)
        float           m_hCalRingOuterPhiCoordinate;           ///< HCal ring outer phi coordinate (missing from ILD gear files)
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     */
     GeometryCreator(const Settings &settings);

    /**
     *  @brief  Destructor
     */
     ~GeometryCreator();

    /**
     *  @brief  Create geometry
     */
    pandora::StatusCode CreateGeometry() const;

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
    pandora::StatusCode SetILDSpecificGeometry(PandoraApi::GeometryParameters &geometryParameters) const;

    /**
     *  @brief  Add information missing from GEAR parameters file for ILD SDHCAL detector
     * 
     *  @param  geometryParameters the pandora geometry parameters
     */
    pandora::StatusCode SetILD_SDHCALSpecificGeometry(PandoraApi::GeometryParameters &geometryParameters) const;

    /**
     *  @brief  Specify positions of hcal barrel box gaps - ILD specific
     */
    pandora::StatusCode CreateHCalBarrelBoxGaps() const;

    /**
     *  @brief  Specify positions of hcal end cap box gaps - ILD specific
     */
    pandora::StatusCode CreateHCalEndCapBoxGaps() const;

    /**
     *  @brief  Specify positions of hcal barrel concentric polygon gaps - ILD specific
     */
    pandora::StatusCode CreateHCalBarrelConcentricGaps() const;

    /**
     *  @brief  Create box gaps at regular positions on polygonal prism, oriented along main z axis - ILD specific
     * 
     *  @param  symmetryOrder the pandora geometry parameters
     *  @param  phi0 the phi coordinate
     *  @param  innerRadius the inner r coordinate
     *  @param  outerRadius the outer r coordinate
     *  @param  minZ the minimum z coordinate
     *  @param  maxZ the maximum z coordinate
     *  @param  gapWidth the gap width
     *  @param  vertexOffset position offset for vertex that doesn't point back to origin of xy plane
     */
    pandora::StatusCode CreateRegularBoxGaps(unsigned int symmetryOrder, float phi0, float innerRadius, float outerRadius, float minZ,
        float maxZ, float gapWidth, pandora::CartesianVector vertexOffset = pandora::CartesianVector(0, 0, 0)) const;

    const Settings          m_settings;                     ///< The geometry creator settings
    const pandora::Pandora *m_pPandora;                     ///< Address of the pandora object to create the geometry
};

#endif // #ifndef GEOMETRY_CREATOR_H
