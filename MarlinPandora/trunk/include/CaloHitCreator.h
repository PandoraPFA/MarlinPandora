/**
 *  @file   PandoraPFANew/include/CaloHitCreator.h
 * 
 *  @brief  Header file for the calo hit creator class.
 * 
 *  $Log: $
 */

#ifndef CALO_HIT_CREATOR_H
#define CALO_HIT_CREATOR_H 1

#include "EVENT/CalorimeterHit.h"

#include "gear/LayerLayout.h"

#include "Api/PandoraApi.h"

typedef std::vector<CalorimeterHit *> CalorimeterHitVector;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  CaloHitCreator class
 */
class CaloHitCreator
{
public:
    typedef std::vector<std::string> StringVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        StringVector    m_eCalCaloHitCollections;               ///< The ecal calorimeter hit collections
        StringVector    m_hCalCaloHitCollections;               ///< The hcal calorimeter hit collections
        StringVector    m_lCalCaloHitCollections;               ///< The lcal calorimeter hit collections
        StringVector    m_lHCalCaloHitCollections;              ///< The lhcal calorimeter hit collections
        StringVector    m_muonCaloHitCollections;               ///< The muon calorimeter hit collections

        float           m_absorberRadiationLength;              ///< The absorber radiation length
        float           m_absorberInteractionLength;            ///< The absorber interaction length
        float           m_eCalToMip;                            ///< The calibration from deposited ECal energy to mip
        float           m_hCalToMip;                            ///< The calibration from deposited HCal energy to mip
        float           m_muonToMip;                            ///< The calibration from deposited Muon energy to mip
        float           m_eCalMipThreshold;                     ///< Threshold for creating calo hits in the ECal, units mip
        float           m_hCalMipThreshold;                     ///< Threshold for creating calo hits in the HCal, units mip
        float           m_muonMipThreshold;                     ///< Threshold for creating calo hits in the HCal, units mip

        float           m_eCalToEMGeV;                          ///< The calibration from deposited ECal energy to EM energy
        float           m_hCalToEMGeV;                          ///< The calibration from deposited HCal energy to EM energy
        float           m_eCalToHadGeV;                         ///< The calibration from deposited ECal energy to hadronic energy
        float           m_hCalToHadGeV;                         ///< The calibration from deposited HCal energy to hadronic energy
        int             m_muonDigitalHits;                      ///< Muon hits are treated as digital (energy from hit count)
        float           m_muonHitEnergy;                        ///< The energy for a digital muon calorimeter hit, units GeV

        float           m_maxHCalHitHadronicEnergy;             ///< The maximum hadronic energy allowed for a single hcal hit
        int             m_nOuterSamplingLayers;                 ///< Number of layers from edge for hit to be flagged as an outer layer hit
        float           m_layersFromEdgeMaxRearDistance;        ///< Maximum number of layers from candidate outer layer hit to rear of detector

        int             m_hCalEndCapInnerSymmetryOrder;         ///< HCal end cap inner symmetry order (missing from ILD00 gear file)
        float           m_hCalEndCapInnerPhiCoordinate;         ///< HCal end cap inner phi coordinate (missing from ILD00 gear file)
    };

    /**
     *  @brief  Create calo hits, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */    
    StatusCode CreateCaloHits(const LCEvent *const pLCEvent);

    /**
     *  @brief  Get the calorimeter hit vector
     * 
     *  @return The calorimeter hit vector
     */
    static const CalorimeterHitVector &GetCalorimeterHitVector();

    /**
     *  @brief  Reset the calo hit creator
     */
    void Reset();

    Settings                m_settings;         ///< The settings

private:
    /**
     *  @brief  Create ecal calo hits
     * 
     *  @param  pLCEvent the lcio event
     */    
    StatusCode CreateECalCaloHits(const LCEvent *const pLCEvent);

    /**
     *  @brief  Create hcal calo hits
     * 
     *  @param  pLCEvent the lcio event
     */
    StatusCode CreateHCalCaloHits(const LCEvent *const pLCEvent);

    /**
     *  @brief  Create lcal calo hits
     * 
     *  @param  pLCEvent the lcio event
     */    
    StatusCode CreateLCalCaloHits(const LCEvent *const pLCEvent);

    /**
     *  @brief  Create lhcal calo hits
     * 
     *  @param  pLCEvent the lcio event
     */
    StatusCode CreateLHCalCaloHits(const LCEvent *const pLCEvent);

    /**
     *  @brief  Create muon calo hits
     * 
     *  @param  pLCEvent the lcio event
     */
    StatusCode CreateMuonCaloHits(const LCEvent *const pLCEvent);

    /**
     *  @brief  Get common calo hit properties: position, parent address, input energy and time
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  caloHitParameters the calo hit parameters to populate
     */
    void GetCommonCaloHitProperties(CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const;

    /**
     *  @brief  Get end cap specific calo hit properties: cell size, absorber radiation and interaction lengths, normal vector
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  layerLayout the gear end cap layer layout
     *  @param  caloHitParameters the calo hit parameters to populate
     *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
     */
    void GetEndCapCaloHitProperties(CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
        PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const;

    /**
     *  @brief  Get barrel specific calo hit properties: cell size, absorber radiation and interaction lengths, normal vector
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  layerLayout the gear barrel layer layout
     *  @param  barrelSymmetryOrder the barrel order of symmetry
     *  @param  barrelPhi0 the barrel orientation
     *  @param  staveNumber the stave number
     *  @param  caloHitParameters the calo hit parameters to populate
     *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
     */
    void GetBarrelCaloHitProperties(CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout, unsigned int barrelSymmetryOrder,
        float barrelPhi0, unsigned int staveNumber, PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const;

    /**
     *  @brief  Get number of active layers from position of a calo hit to the edge of the detector
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     */
    int GetNLayersFromEdge(CalorimeterHit *const pCaloHit) const;

    /**
     *  @brief  Get the maximum radius of a calo hit in a polygonal detector structure
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  symmetryOrder the symmetry order
     *  @param  phi0 the angular orientation
     * 
     *  @return the maximum radius
     */
    float GetMaximumRadius(CalorimeterHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const;

    static CalorimeterHitVector        m_calorimeterHitVector;     ///< The calorimeter hit vector
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CalorimeterHitVector &CaloHitCreator::GetCalorimeterHitVector()
{
    return m_calorimeterHitVector;
}


//------------------------------------------------------------------------------------------------------------------------------------------

inline void CaloHitCreator::Reset()
{
    m_calorimeterHitVector.clear();
}

#endif // #ifndef CALO_HIT_CREATOR_H
