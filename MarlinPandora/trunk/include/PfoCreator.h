/**
 *  @file   MarlinPandora/include/PfoCreator.h
 * 
 *  @brief  Header file for the pfo creator class.
 * 
 *  $Log: $
 */

#ifndef PFO_CREATOR_H
#define PFO_CREATOR_H 1

#include "EVENT/LCEvent.h"

#include "Api/PandoraApi.h"

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
        std::string     m_pfoCollectionName;                    ///< The name of the pfo output collection
        std::string     m_reclusterMonitoringCollectionName;    ///< The name of the recluster monitoring collection
        std::string     m_reclusterRelationsCollectionName;     ///< The name of the collection relating recluster monitoring to tracks
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     */
     PfoCreator(const Settings &settings);

    /**
     *  @brief  Destructor
     */
     ~PfoCreator();

    /**
     *  @brief  Create particle flow objects
     * 
     *  @param  pLCEvent the lcio event
     */    
    pandora::StatusCode CreateParticleFlowObjects(EVENT::LCEvent *pLCEvent);

private:
    const Settings          m_settings;                         ///< The pfo creator settings
    const pandora::Pandora *m_pPandora;                         ///< Address of the pandora object from which to extract the pfos
};

#endif // #ifndef PFO_CREATOR_H
