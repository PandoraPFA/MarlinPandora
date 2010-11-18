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
     *  @brief  Compute the path length of the intersection of the line from the IP to the position of a CalorimeterHit with a rectangle
     * 
     *  @param  position position of the calorimeter-hit
     *  @param  rMin minimum radius coordinate of the rectangle (assuming zylindrical coordinates)
     *  @param  zMin minimum z-position coordinate of the rectangle
     *  @param  rMax maximum radius coordinate of the rectangle (assuming zylindrical coordinates)
     *  @param  zMax maximum z-position coordinate of the rectangle
     * 
     *  @return length of the line segment within the rectangle
     */
    static float ComputePathLengthFromIPInRectangle(const pandora::CartesianVector &position, float rMin, float zMin, float rMax, float zMax);

    /**
     *  @brief  Compute the intersection of two lines
     * 
     *  @param  lineAXStart x coordinate of the start of the first line
     *  @param  lineAYStart y coordinate of the start of the first line
     *  @param  lineAXEnd x coordinate of the end of the first line
     *  @param  lineAYEnd y coordinate of the end of the first line
     *  @param  lineBXStart x coordinate of the start of the second line
     *  @param  lineBYStart y coordinate of the start of the second line
     *  @param  lineBXEnd x coordinate of the end of the second line
     *  @param  lineBYEnd y coordinate of the end of the second line
     *  @param  xIntersect to receive the x coordinate of the intersection point
     *  @param  yIntersect to receive the y coordinate of the intersection point
     * 
     *  @return true if the lines are intersecting within their respective start and end points
     */
    static bool IntersectLines2D(float lineAXStart, float lineAYStart, float lineAXEnd, float lineAYEnd, float lineBXStart,
        float lineBYStart, float lineBXEnd, float lineBYEnd, float &xIntersect, float &yIntersect);

    static bool                     m_instanceFlag;                ///< The path length calculator instance flag
    static PathLengthCalculator    *m_pPathLengthCalculator;       ///< The path length calculator instance
};

#endif // #ifndef PATH_LENGTH_CALCULATOR_H
