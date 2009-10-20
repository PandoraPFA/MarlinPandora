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

#include "Pandora/Pandora.h"

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
        std::string     m_pandoraSettingsXmlFile;       ///< The pandora settings xml file

        StringVector    m_trackCollections;             ///< The reconstructed track collections
        StringVector    m_v0VertexCollections;          ///< The v0 vertex collections
        StringVector    m_caloHitCollections;           ///< The calorimeter hit collections
        StringVector    m_mcParticleCollections;        ///< The mc particle collections
        StringVector    m_lcCaloHitRelationCollections; ///< The SimCaloHit to CaloHit particle relations
        StringVector    m_lcTrackRelationCollections;   ///< The SimTrackerHit to TrackerHit particle relations

        std::string     m_particleCollectionName;       ///< The name of the PFO output collection

        float           m_absorberRadiationLength;      ///< The absorber radation length
        float           m_absorberInteractionLength;    ///< The absorber interaction length

        float           m_eCalToMip;                    ///< The calibration from deposited ECal energy to mip
        float           m_hCalToMip;                    ///< The calibration from deposited HCal energy to mip
        float           m_eCalMipThreshold;             ///< Threshold for creating calo hits in the ECal, units mip
        float           m_hCalMipThreshold;             ///< Threshold for creating calo hits in the HCal, units mip

        float           m_eCalToEMGeV;                  ///< The calibration from deposited ECal energy to EM energy
        float           m_hCalToEMGeV;                  ///< The calibration from deposited HCal energy to EM energy
        float           m_eCalToHadGeV;                 ///< The calibration from deposited ECal energy to hadronic energy
        float           m_hCalToHadGeV;                 ///< The calibration from deposited HCal energy to hadronic energy
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

private:
    /**
     *  @brief  Create geometry, insert user code here
     */
    StatusCode CreateGeometry() const;

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
     *  @brief  Register user algorithm factories, insert user code here
     */
    StatusCode RegisterUserAlgorithmFactories() const;

    /**
     *  @brief  Create tracks, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */    
    StatusCode CreateTracks(const LCEvent *const pLCEvent);

    /**
     *  @brief  Create MCParticles, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */    
    StatusCode CreateMCParticles(const LCEvent *const pLCEvent) const;

    /**
     *  @brief  Create calo hits, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */    
    StatusCode CreateCaloHits(const LCEvent *const pLCEvent);

    /**
     *  @brief  Create calo hit to mc particle relationships
     *
     *  @param  pLCEvent the lcio event
     */
    StatusCode CreateCaloHitToMCParticleRelationships(const LCEvent *const pLCEvent) const;

    /**
     *  @brief  Create Track to mc particle relationships
     *
     *  @param  pLCEvent the lcio event
     */
    StatusCode CreateTrackToMCParticleRelationships(const LCEvent *const pLCEvent) const;

    /**
     *  @brief  Process particle flow objects, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */    
    StatusCode ProcessParticleFlowObjects( LCEvent * pLCEvent);

    /**
     *  @brief  Process steering file parameters, insert user code here
     */
    void ProcessSteeringFile();

    typedef std::vector<CalorimeterHit *> CalorimeterHitVector;
    typedef std::vector<Track *>          TrackVector;

    pandora::Pandora            m_pandora;                  ///< The pandora instance
    Settings                    m_settings;                 ///< The settings for the pandora pfa new processor
    std::string                 m_detectorName;             ///< The detector name
    unsigned int                m_nRun;                     ///< The run number
    unsigned int                m_nEvent;                   ///< The event number
    CalorimeterHitVector        m_calorimeterHitVector;     ///< The calorimeter hit vector
    TrackVector                 m_trackVector;              ///< The track vector
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline marlin::Processor *PandoraPFANewProcessor::newProcessor()
{
    return new PandoraPFANewProcessor;
}

#endif // #ifndef PANDORA_PFA_NEW_PROCESSOR_H
