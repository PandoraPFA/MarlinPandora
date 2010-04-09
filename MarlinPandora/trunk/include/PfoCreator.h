/**
 *  @file   PandoraPFANew/include/PfoCreator.h
 * 
 *  @brief  Header file for the pfo creator class.
 * 
 *  $Log: $
 */

#ifndef PFO_CREATOR_H
#define PFO_CREATOR_H 1

#include "EVENT/LCEvent.h"

#include "Api/PandoraApi.h"

using namespace EVENT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  PfoCreator class
 */
class PfoCreator
{
public:
    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        std::string     m_clusterCollectionName;                ///< The name of the cluster output collection
        std::string     m_pfoCollectionName;                    ///< The name of the PFO output collection
    };

    /**
     *  @brief  Create particle flow objects, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */    
    StatusCode CreateParticleFlowObjects(LCEvent *pLCEvent);

    Settings                m_settings;         ///< The settings
};

#endif // #ifndef PFO_CREATOR_H
