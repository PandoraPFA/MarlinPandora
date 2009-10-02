/**
 *  @file   PandoraPFANew/include/PandoraPFANewProcessor.h
 * 
 *  @brief  Header file for the pandora pfa new processor class.
 * 
 *  $Log: $
 */
#ifndef PANDORA_PFA_NEW_PROCESSOR_H
#define PANDORA_PFA_NEW_PROCESSOR_H 1

#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include "EVENT/LCCollection.h"

#include <algorithm>

#include "marlin/Processor.h"

#include "Pandora/Pandora.h"

#include "Test/TestMCManager.h"



/* // for sorting a vector of pairs with an algorithm by "second" */
/*
/* template< typename T1, typename T2 > */
/*    inline bool less_than_second( const std::pair<T1,T2>& b1, const std::pair<T1,T2>& b2 ){ */
/*    return b1.second < b2.second; */
/* } */



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
        StringVector    m_lcRelationCollections;        ///< The caloHit to MC particle relations

        float           m_absorberRadiationLength;      ///< The absorber radation length
        float           m_absorberInteractionLength;    ///< The absorber interaction length
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
    StatusCode CreateTracks(const LCEvent *const pLCEvent) const;

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
    StatusCode CreateCaloHits(const LCEvent *const pLCEvent) const;

    /**
     *  @brief  Process particle flow objects, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */    
    StatusCode ProcessParticleFlowObjects(const LCEvent *const pLCEvent) const;

    /**
     *  @brief  Process steering file parameters, insert user code here
     */
    void ProcessSteeringFile();


    /**
     *  @brief  Get the MCParticles which contributed to a CaloHit ordered by energy
     *          (the MCParticle .at(0) is the one with the highest contribution)
     *
     *  @param  pLCEvent the lcio event
     *  @param  pCaloHit points to CalorimeterHit of which the MCParticles are requested
     *  @param  pMCParticles is a reference to a vector of pairs which gets the pointers 
     *          to the MCParticles together with their relative energy contributions to the hit
     */
    void GetCaloHitMCParticles( const LCEvent *const pLCEvent,
                                CalorimeterHit* pCaloHit, 
                                std::vector<std::pair<MCParticle*,double> >& pMcParticles ) const;

    

    pandora::Pandora    m_pandora;          ///< The pandora instance
    Settings            m_settings;         ///< The settings for the pandora pfa new processor
    std::string         m_detectorName;     ///< The detector name
    unsigned int        m_nRun;             ///< The run number
    unsigned int        m_nEvent;           ///< The event number
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline marlin::Processor *PandoraPFANewProcessor::newProcessor()
{
    return new PandoraPFANewProcessor;
}

#endif // #ifndef PANDORA_PFA_NEW_PROCESSOR_H
