/**
 *  @file   MarlinPandora/include/PathLengthCalculator.h
 * 
 *  @brief  Header file for the path length calculator class.
 * 
 *  $Log: $
 */

#ifndef PATH_LENGTH_CALCULATOR_H
#define PATH_LENGTH_CALCULATOR_H 1

#include "EVENT/CalorimeterHit.h"

#include "Api/PandoraApi.h"

#include "Helpers/GeometryHelper.h"

/**
 *  @brief  PathLengthCalculator class
 */
class PathLengthCalculator
{
public:
    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        static float        m_avgRadLengthCoil;                     ///< Average radiation length per mm in the coil
        static float        m_avgRadLengthECalBarrel;               ///< Average radiation length per mm in the ECal barrel
        static float        m_avgRadLengthHCalBarrel;               ///< Average radiation length per mm in the HCal barrel
        static float        m_avgRadLengthECalEndCap;               ///< Average radiation length per mm in the ECal endcap
        static float        m_avgRadLengthHCalEndCap;               ///< Average radiation length per mm in the HCal endcap
        static float        m_avgRadLengthMuonBarrel;               ///< Average radiation length per mm in the Muon barrel
        static float        m_avgRadLengthMuonEndCap;               ///< Average radiation length per mm in the Muon endcap

        static float        m_avgIntLengthCoil;                     ///< Average interaction length per mm in the coil
        static float        m_avgIntLengthECalBarrel;               ///< Average interaction length per mm in the ECal barrel
        static float        m_avgIntLengthHCalBarrel;               ///< Average interaction length per mm in the HCal barrel
        static float        m_avgIntLengthECalEndCap;               ///< Average interaction length per mm in the ECal endcap
        static float        m_avgIntLengthHCalEndCap;               ///< Average interaction length per mm in the HCal endcap
        static float        m_avgIntLengthMuonBarrel;               ///< Average interaction length per mm in the Muon barrel
        static float        m_avgIntLengthMuonEndCap;               ///< Average interaction length per mm in the Muon endcap
    };

    typedef std::vector<pandora::CartesianVector> NormalVectorList;

    /**
     *  @brief  SubDetectorParameters class
     */
    class SubDetectorParameters
    {
    public:
        /**
         *  @brief  Initialize the subdetector parameters
         *
         *  @param  subDetectorParameters the pandora geometry helper subdetector parameters
         */
        void Initialize(const pandora::GeometryHelper::SubDetectorParameters &subDetectorParameters);

        /**
         *  @brief  Initialize the subdetector parameters for a cylindrical subdetector
         *
         *  @param  innerRCoordinate the inner r coordinate
         *  @param  outerRCoordinate the outer r coordinate
         *  @param  innerZCoordinate the inner z coordinate
         *  @param  outerZCoordinate the outer z coordinate
         */
        void Initialize_Cylinder(const float innerRCoordinate, const float outerRCoordinate, const float innerZCoordinate,
            const float outerZCoordinate);

        /**
         *  @brief  Get the subdetector inner r coordinate
         * 
         *  @return The subdetector inner r coordinate
         */
        float GetInnerRCoordinate() const;

        /**
         *  @brief  Get the subdetector outer r coordinate
         * 
         *  @return The subdetector outer r coordinate
         */
        float GetOuterRCoordinate() const;

        /**
         *  @brief  Get the subdetector inner z coordinate
         * 
         *  @return The subdetector inner z coordinate
         */
        float GetInnerZCoordinate() const;

        /**
         *  @brief  Get the subdetector outer z coordinate
         * 
         *  @return The subdetector outer z coordinate
         */
        float GetOuterZCoordinate() const;

        /**
         *  @brief  Get the subdetector inner r normal vectors, empty for cylinders
         * 
         *  @return The subdetector inner r normal vectors
         */
        const NormalVectorList &GetInnerRNormalVectors() const;

        /**
         *  @brief  Get the subdetector outer r normal vectors, empty for cylinders
         * 
         *  @return The subdetector outer r normal vectors
         */
        const NormalVectorList &GetOuterRNormalVectors() const;

        /**
         *  @brief  Get the subdetector z normal vectors
         * 
         *  @return The subdetector z normal vectors
         */
        const NormalVectorList &GetZNormalVectors() const;

    private:
        /**
         *  @brief  Get list of vectors, normal to sides of a specified regular polygon
         * 
         *  @param  symmetry the polgon order of symmetry
         *  @param  phi0 the polygon phi0 coordinate
         *  @param  normalVectorList to receive the normal vector list
         */
        void GetPolygonNormalVectors(const unsigned int symmetry, const float phi0, NormalVectorList &normalVectorList) const;

        float                       m_innerRCoordinate;             ///< The subdetector inner r coordinate
        float                       m_outerRCoordinate;             ///< The subdetector outer r coordinate
        float                       m_innerZCoordinate;             ///< The subdetector inner z coordinate
        float                       m_outerZCoordinate;             ///< The subdetector outer z coordinate
        NormalVectorList            m_innerRNormalVectors;          ///< The subdetector inner r normal vectors; leave empty for cylinders
        NormalVectorList            m_outerRNormalVectors;          ///< The subdetector outer r normal vectors; leave empty for cylinders
        NormalVectorList            m_zNormalVectors;               ///< The subdetector z normal vectors
    };

    /**
     *  @brief  Get the interaction length calculator singleton
     */
    static PathLengthCalculator *GetInstance();

    /**
     *  @brief  Destructor
     */
    ~PathLengthCalculator();

    /**
     *  @brief  Get the path length from the ip to the position of a calorimeter hit in units of interaction lengths
     * 
     *  @param  pCaloHit address of the calorimeter hit
     *  @param  nRadiationLengthsFromIp to receive the path length in radiation lengths
     *  @param  nInteractionLengthsFromIp to receive the path length in interaction lengths
     */
    static void GetPathLengths(const EVENT::CalorimeterHit *const pCaloHit, float &nRadiationLengthsFromIp, float &nInteractionLengthsFromIp);

private:
    /**
     *  @brief  Constructor
     */
    PathLengthCalculator();

    /**
     *  @brief  Get the fraction of a straight line path (from the ip to a specified position vector) represented by its traversal
     *          through a specified subdetector
     * 
     *  @param  positionVector the position vector
     *  @param  subDetectorParameters the subdetector parameters
     * 
     *  @return The fraction of the straight line path represented by its traversal through the subdetector
     */
    static float GetFractionInSubDetector(const pandora::CartesianVector &positionVector, const SubDetectorParameters &subDetectorParameters);

    /**
     *  @brief  Get the fraction of a straight line path (from the ip to a specified position vector) represented by the projection
     *          onto a vector of length "closestDistanceToIp" and direction normal to the most appropriate side of a regular polygon.
     * 
     *  @param  positionVector the position vector
     *  @param  closestDistanceToIp the closest distance to the ip
     *  @param  normalVectorList the list of unit vectors normal to sides of a regular polygon
     * 
     *  @return The fraction of the straight line path represented by the projection
     */
    static float GetLengthFraction(const pandora::CartesianVector &positionVector, const float closestDistanceToIp,
        const NormalVectorList &normalVectorList);

    /**
     *  @brief  Initialize the path length calculator geometry
     */
    static void InitializeGeometry();

    static bool                     m_instanceFlag;                 ///< The path length calculator instance flag
    static bool                     m_isGeometryInitialized;        ///< Whether the path length calculator geometry has been initialized
    static PathLengthCalculator    *m_pPathLengthCalculator;        ///< The path length calculator instance

    static SubDetectorParameters    m_eCalBarrelParameters;         ///< The ecal barrel parameters
    static SubDetectorParameters    m_hCalBarrelParameters;         ///< The hcal barrel parameters
    static SubDetectorParameters    m_muonBarrelParameters;         ///< The muon barrel parameters
    static SubDetectorParameters    m_eCalEndCapParameters;         ///< The ecal endcap parameters
    static SubDetectorParameters    m_hCalEndCapParameters;         ///< The hcal endcap parameters
    static SubDetectorParameters    m_muonEndCapParameters;         ///< The muon endcap parameters
    static SubDetectorParameters    m_coilParameters;               ///< The coil parameters
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float PathLengthCalculator::SubDetectorParameters::GetInnerRCoordinate() const
{
    return m_innerRCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PathLengthCalculator::SubDetectorParameters::GetOuterRCoordinate() const
{
    return m_outerRCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PathLengthCalculator::SubDetectorParameters::GetInnerZCoordinate() const
{
    return m_innerZCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PathLengthCalculator::SubDetectorParameters::GetOuterZCoordinate() const
{
    return m_outerZCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const PathLengthCalculator::NormalVectorList &PathLengthCalculator::SubDetectorParameters::GetInnerRNormalVectors() const
{
    return m_innerRNormalVectors;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const PathLengthCalculator::NormalVectorList &PathLengthCalculator::SubDetectorParameters::GetOuterRNormalVectors() const
{
    return m_outerRNormalVectors;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const PathLengthCalculator::NormalVectorList &PathLengthCalculator::SubDetectorParameters::GetZNormalVectors() const
{
    return m_zNormalVectors;
}

#endif // #ifndef PATH_LENGTH_CALCULATOR_H
