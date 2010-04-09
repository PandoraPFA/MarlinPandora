/**
 *  @file   PandoraPFANew/include/PandoraPFANewProcessor.h
 * 
 *  @brief  Header file for the pandora pfa new processor class.
 * 
 *  $Log: $
 */

#ifndef PANDORA_PFA_NEW_PROCESSOR_H
#define PANDORA_PFA_NEW_PROCESSOR_H 1

#include "marlin/Processor.h"

#include "CaloHitCreator.h"
#include "GeometryCreator.h"
#include "MCParticleCreator.h"
#include "PfoCreator.h"
#include "TrackCreator.h"

namespace pandora {class Pandora;}

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  PandoraPFANewProcessor class
 */
class PandoraPFANewProcessor : public marlin::Processor
{
public:
    typedef std::vector<std::string> StringVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        std::string     m_pandoraSettingsXmlFile;           ///< The pandora settings xml file
        int             m_nEventsToSkip;                    ///< Number of events to skip at start of reconstruction job
    };

    /**
     *  @brief  Default constructor
     */
    PandoraPFANewProcessor();

    /**
     *  @brief  Create new processor
     */
    virtual Processor *newProcessor();

    /**
     *  @brief  Initialize, called at startup
     */
    virtual void init();

    /**
     *  @brief  Process run header
     *
     *  @param  pLCRunHeader the lc run header
     */
    virtual void processRunHeader(lcio::LCRunHeader *pLCRunHeader);

    /**
     *  @brief  Process event, main entry point
     *
     *  @param  pLCEvent the lc event
     */
    virtual void processEvent(lcio::LCEvent *pLCEvent);

    /**
     *  @brief  Checks for event
     *
     *  @param  pLCEvent the lc event
     */
    virtual void check(lcio::LCEvent *pLCEvent);

    /**
     *  @brief  End, called at shutdown
     */
    virtual void end();

    /**
     *  @brief  Get address of the pandora instance
     * 
     *  @return address of the pandora instance
     */
    static pandora::Pandora *GetPandora();

private:
    /**
     *  @brief  Register user algorithm factories, insert user code here
     */
    StatusCode RegisterUserAlgorithmFactories() const;

    /**
     *  @brief  Process steering file parameters, insert user code here
     */
    void ProcessSteeringFile();

    /**
     *  @brief  Copy some steering parameters between settings objects
     */
    void FinaliseSteeringParameters();

    static pandora::Pandora    *m_pPandora;                 ///< Address of the pandora instance
    Settings                    m_settings;                 ///< The settings for the pandora pfa new processor
    std::string                 m_detectorName;             ///< The detector name
    unsigned int                m_nRun;                     ///< The run number
    unsigned int                m_nEvent;                   ///< The event number

    CaloHitCreator              m_caloHitCreator;           ///< The calo hit creator
    GeometryCreator             m_geometryCreator;          ///< The geometry creator
    MCParticleCreator           m_mcParticleCreator;        ///< The mc particle creator
    PfoCreator                  m_pfoCreator;               ///< The pfo creator
    TrackCreator                m_trackCreator;             ///< The track creator
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline marlin::Processor *PandoraPFANewProcessor::newProcessor()
{
    return new PandoraPFANewProcessor;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Pandora *PandoraPFANewProcessor::GetPandora()
{
    if (NULL == m_pPandora)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return m_pPandora;
}

#endif // #ifndef PANDORA_PFA_NEW_PROCESSOR_H
