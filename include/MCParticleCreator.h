/**
 *  @file   PandoraPFANew/include/MCParticleCreator.h
 * 
 *  @brief  Header file for the mc particle creator class.
 * 
 *  $Log: $
 */

#ifndef MC_PARTICLE_CREATOR_H
#define MC_PARTICLE_CREATOR_H 1

#include "EVENT/LCEvent.h"

#include "Api/PandoraApi.h"

using namespace EVENT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MCParticleCreator class
 */
class MCParticleCreator
{
public:
    typedef std::vector<std::string> StringVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        StringVector    m_mcParticleCollections;                ///< The mc particle collections
        StringVector    m_lcCaloHitRelationCollections;         ///< The SimCaloHit to CaloHit particle relations
        StringVector    m_lcTrackRelationCollections;           ///< The SimTrackerHit to TrackerHit particle relations
    };

    /**
     *  @brief  Create MCParticles, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */    
    StatusCode CreateMCParticles(const LCEvent *const pLCEvent) const;

    /**
     *  @brief  Create Track to mc particle relationships
     *
     *  @param  pLCEvent the lcio event
     */
    StatusCode CreateTrackToMCParticleRelationships(const LCEvent *const pLCEvent) const;

    /**
     *  @brief  Create calo hit to mc particle relationships
     *
     *  @param  pLCEvent the lcio event
     */
    StatusCode CreateCaloHitToMCParticleRelationships(const LCEvent *const pLCEvent) const;

    Settings                m_settings;         ///< The settings
};

#endif // #ifndef MC_PARTICLE_CREATOR_H
