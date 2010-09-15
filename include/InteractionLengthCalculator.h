/**
 *  @file   PandoraPFANew/include/InteractionLengthCalculator.h
 * 
 *  @brief  Header file for the interaction length calculator class.
 * 
 *  $Log: $
 */

#ifndef INTERACTION_LENGTH_CALCULATOR_H
#define INTERACTION_LENGTH_CALCULATOR_H 1

#include "EVENT/CalorimeterHit.h"

#include "Api/PandoraApi.h"

/**
 *  @brief  InteractionLengthCalculator class
 */
class InteractionLengthCalculator
{
public:
    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        static float        m_avgIntLengthTracker;                  ///< Average interaction length per mm in the tracker
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
    static InteractionLengthCalculator *GetInstance();

    /**
     *  @brief  Destructor
     */
    ~InteractionLengthCalculator();

    /**
     *  @brief  Compute the path length from the IP to the position of a CalorimeterHit in units of interaction lengths
     * 
     *  @param  pCaloHit Calorimeter hit
     *  @param  nInteractionLengths to receive the length from the IP to the position of the calorimeter hit in units of interaction length
     */
    static StatusCode ComputeInteractionLengthsFromIP(const EVENT::CalorimeterHit *const pCaloHit, float &nInteractionLengths);

private:
    /**
     *  @brief  Constructor
     */
    InteractionLengthCalculator();

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

    static bool                             m_instanceFlag;                         ///< The interaction length calculator instance flag
    static InteractionLengthCalculator     *m_pInteractionLengthCalculator;         ///< The interaction length calculator instance
};

#endif // #ifndef INTERACTION_LENGTH_CALCULATOR_H
