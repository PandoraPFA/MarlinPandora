/**
 *  @file   MarlinPandora/include/PandoraPFANewProcessor.h
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
    typedef std::vector<float> FloatVector;
    typedef std::vector<std::string> StringVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Settings();

        std::string     m_pandoraSettingsXmlFile;           ///< The pandora settings xml file

        float           m_innerBField;                      ///< The bfield in the main tracker, ecal and hcal, units Tesla
        float           m_muonBarrelBField;                 ///< The bfield in the muon barrel, units Tesla
        float           m_muonEndCapBField;                 ///< The bfield in the muon endcap, units Tesla

        FloatVector     m_inputEnergyCorrectionPoints;      ///< The input energy points for non-linearity energy correction
        FloatVector     m_outputEnergyCorrectionPoints;     ///< The output energy points for non-linearity energy correction
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
    virtual void processEvent(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  Checks for event
     *
     *  @param  pLCEvent the lc event
     */
    virtual void check(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  End, called at shutdown
     */
    virtual void end();

    /**
     *  @brief  Get address of the pandora instance
     * 
     *  @return address of the pandora instance
     */
    const pandora::Pandora *GetPandora() const;

    /**
     *  @brief  Get address of the current lcio event
     * 
     *  @param  pPandora address of the relevant pandora instance
     * 
     *  @return address of the current lcio event
     */
    static const EVENT::LCEvent *GetCurrentEvent(const pandora::Pandora *const pPandora);

private:
    /**
     *  @brief  Register user algorithm factories, energy correction functions and particle id functions,
     *          insert user code here
     */
    pandora::StatusCode RegisterUserComponents() const;

    /**
     *  @brief  Process steering file parameters, insert user code here
     */
    void ProcessSteeringFile();

    /**
     *  @brief  Copy some steering parameters between settings objects
     */
    void FinaliseSteeringParameters();

    /**
     *  @brief  Reset the pandora pfa new processor
     */
    void Reset();

    pandora::Pandora                   *m_pPandora;                         ///< Address of the pandora instance
    CaloHitCreator                     *m_pCaloHitCreator;                  ///< The calo hit creator
    GeometryCreator                    *m_pGeometryCreator;                 ///< The geometry creator
    TrackCreator                       *m_pTrackCreator;                    ///< The track creator
    MCParticleCreator                  *m_pMCParticleCreator;               ///< The mc particle creator
    PfoCreator                         *m_pPfoCreator;                      ///< The pfo creator

    Settings                            m_settings;                         ///< The settings for the pandora pfa new processor
    CaloHitCreator::Settings            m_caloHitCreatorSettings;           ///< The calo hit creator settings
    GeometryCreator::Settings           m_geometryCreatorSettings;          ///< The geometry creator settings
    MCParticleCreator::Settings         m_mcParticleCreatorSettings;        ///< The mc particle creator settings
    TrackCreator::Settings              m_trackCreatorSettings;             ///< The track creator settings
    PfoCreator::Settings                m_pfoCreatorSettings;               ///< The pfo creator settings

    typedef std::map<const pandora::Pandora *, EVENT::LCEvent *> PandoraToLCEventMap;
    static PandoraToLCEventMap          m_pandoraToLCEventMap;              ///< The pandora to lc event map
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline marlin::Processor *PandoraPFANewProcessor::newProcessor()
{
    return new PandoraPFANewProcessor;
}

#endif // #ifndef PANDORA_PFA_NEW_PROCESSOR_H
