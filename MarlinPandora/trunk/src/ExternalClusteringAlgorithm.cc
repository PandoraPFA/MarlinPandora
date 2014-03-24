/**
 *  @file   MarlinPandora/src/ExternalClusteringAlgorithm.cc
 * 
 *  @brief  Implementation of the external clustering algorithm class.
 * 
 *  $Log: $
 */

#include "EVENT/LCCollection.h"
#include "EVENT/LCEvent.h"
#include "EVENT/Cluster.h"

#include "ExternalClusteringAlgorithm.h"
#include "PandoraPFANewProcessor.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

StatusCode ExternalClusteringAlgorithm::Run()
{
    try
    {
        const CaloHitList *pCaloHitList = NULL;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

        if (pCaloHitList->empty())
            return STATUS_CODE_SUCCESS;

        // Get external photon cluster collection
        const EVENT::LCEvent *const pLCEvent(PandoraPFANewProcessor::GetCurrentEvent());

        const EVENT::LCCollection *pExternalClusterCollection = pLCEvent->getCollection(m_externalClusterCollectionName);
        const unsigned int nExternalClusters(pExternalClusterCollection->getNumberOfElements());

        if (0 == nExternalClusters)
            return STATUS_CODE_SUCCESS;

        // Populate pandora parent address to calo hit map
        ParentAddressToCaloHitMap parentAddressToCaloHitMap;

        for (CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end(); hitIter != hitIterEnd; ++hitIter)
        {
            pandora::CaloHit *pCaloHit = *hitIter;
            parentAddressToCaloHitMap.insert(ParentAddressToCaloHitMap::value_type(pCaloHit->GetParentCaloHitAddress(), pCaloHit));
        }

        // Recreate external clusters within the pandora framework
        for (unsigned int iCluster = 0; iCluster < nExternalClusters; ++iCluster)
        {
            EVENT::Cluster *pExternalCluster = dynamic_cast<EVENT::Cluster *>(pExternalClusterCollection->getElementAt(iCluster));
            const CalorimeterHitVec &calorimeterHitVec(pExternalCluster->getCalorimeterHits());

            if (NULL == pExternalCluster)
                throw EVENT::Exception("Collection type mismatch");

            pandora::Cluster *pPandoraCluster = NULL;

            for (CalorimeterHitVec::const_iterator iter = calorimeterHitVec.begin(), iterEnd = calorimeterHitVec.end(); iter != iterEnd; ++iter)
            {
                ParentAddressToCaloHitMap::const_iterator pandoraCaloHitIter = parentAddressToCaloHitMap.find(*iter);

                if (parentAddressToCaloHitMap.end() == pandoraCaloHitIter)
                {
                    continue;
                }

                pandora::CaloHit *pPandoraCaloHit = pandoraCaloHitIter->second;

                if (NULL == pPandoraCluster)
                {
                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, pPandoraCaloHit, pPandoraCluster));

                    if (m_flagClustersAsPhotons)
                        pPandoraCluster->SetIsFixedPhotonFlag(true);
                }
                else
                {
                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pPandoraCluster, pPandoraCaloHit));
                }
            }
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        return statusCodeException.GetStatusCode();
    }
    catch (EVENT::Exception &exception)
    {
        std::cout << "ExternalClusteringAlgorithm failure: " << exception.what() << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ExternalClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ExternalClusterCollectionName", m_externalClusterCollectionName));

    m_flagClustersAsPhotons = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FlagClustersAsPhotons", m_flagClustersAsPhotons));

    return STATUS_CODE_SUCCESS;
}
