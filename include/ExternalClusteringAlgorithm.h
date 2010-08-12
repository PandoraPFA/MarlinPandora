/**
 *  @file   MarlinPandora/include/Algorithms/ExternalClusteringAlgorithm.h
 * 
 *  @brief  Header file for the external clustering algorithm class.
 * 
 *  $Log: $
 */
#ifndef EXTERNAL_CLUSTERING_ALGORITHM_H
#define EXTERNAL_CLUSTERING_ALGORITHM_H 1

#include "Algorithms/Algorithm.h"

/**
 *  @brief  ExternalClusteringAlgorithm class
 */
class ExternalClusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        Algorithm *CreateAlgorithm() const;
    };

private:
    StatusCode Run();
    StatusCode ReadSettings(const TiXmlHandle xmlHandle);

    typedef std::map<const void *, pandora::CaloHit *> ParentAddressToCaloHitMap;

    std::string     m_externalClusterCollectionName;        ///< The collection name for the external clusters
    bool            m_flagClustersAsPhotons;                ///< Whether to automatically flag new clusters as fixed photons
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ExternalClusteringAlgorithm::Factory::CreateAlgorithm() const
{
    return new ExternalClusteringAlgorithm();
}

#endif // #ifndef EXTERNAL_CLUSTERING_ALGORITHM_H
