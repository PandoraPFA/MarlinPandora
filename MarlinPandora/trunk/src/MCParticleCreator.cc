/**
 *  @file   PandoraPFANew/src/MCParticleCreator.cc
 * 
 *  @brief  Implementation of the mc particle creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/SimCalorimeterHit.h"

#include "UTIL/LCRelationNavigator.h"

#include "gear/BField.h"

#include "CaloHitCreator.h"
#include "MCParticleCreator.h"
#include "PandoraPFANewProcessor.h"
#include "TrackCreator.h"

#include <cmath>
#include <limits>

StatusCode MCParticleCreator::CreateMCParticles(const LCEvent *const pLCEvent) const
{
    // Insert user code here ...
    static pandora::Pandora *pPandora = PandoraPFANewProcessor::GetPandora();

    for (StringVector::const_iterator iter = m_settings.m_mcParticleCollections.begin(), iterEnd = m_settings.m_mcParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pMCParticleCollection = pLCEvent->getCollection(*iter);

            for (int i = 0, iMax = pMCParticleCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    MCParticle *pMcParticle = dynamic_cast<MCParticle*>(pMCParticleCollection->getElementAt(i));

                    double innerRadius = 0.;
                    double outerRadius = 0.;
                    pandora::CartesianVector momentum(pMcParticle->getMomentum()[0], pMcParticle->getMomentum()[1], pMcParticle->getMomentum()[2]);

                    for(int i = 0; i < 3; ++i)
                    {
                        innerRadius += pow(pMcParticle->getVertex()[i], 2);
                        outerRadius += pow(pMcParticle->getEndpoint()[i], 2);
                    }

                    innerRadius = std::sqrt(innerRadius);
                    outerRadius = std::sqrt(outerRadius);
         
                    PandoraApi::MCParticle::Parameters mcParticleParameters;
                    mcParticleParameters.m_energy = pMcParticle->getEnergy();
                    mcParticleParameters.m_particleId = pMcParticle->getPDG();
                    mcParticleParameters.m_momentum = momentum;
                    mcParticleParameters.m_innerRadius = innerRadius;
                    mcParticleParameters.m_outerRadius = outerRadius;
                    mcParticleParameters.m_pParentAddress = pMcParticle;

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPandora, mcParticleParameters));

                    // Create parent-daughter relationships
                    for(MCParticleVec::const_iterator itDaughter = pMcParticle->getDaughters().begin(),
                        itDaughterEnd = pMcParticle->getDaughters().end(); itDaughter != itDaughterEnd; ++itDaughter)
                    {
                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*pPandora, pMcParticle,
                            *itDaughter));
                    }
                }
                catch (StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract MCParticle: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract MCParticle, unrecognised exception" << std::endl;
                }
            }
        }
        catch (...)
        {
            streamlog_out(WARNING) << "Failed to extract MCParticles collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCParticleCreator::CreateTrackToMCParticleRelationships(const LCEvent *const pLCEvent) const
{
    static pandora::Pandora *pPandora = PandoraPFANewProcessor::GetPandora();
    static const float bField(marlin::Global::GEAR->getBField().at(gear::Vector3D(0., 0., 0.)).z());

    for (StringVector::const_iterator iter = m_settings.m_lcTrackRelationCollections.begin(), iterEnd = m_settings.m_lcTrackRelationCollections.end();
         iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pMCRelationCollection = pLCEvent->getCollection(*iter);
            LCRelationNavigator navigate(pMCRelationCollection);

            const TrackVector &trackVector = TrackCreator::GetTrackVector();

            for (TrackVector::const_iterator trackIter = trackVector.begin(), trackIterEnd = trackVector.end();
                trackIter != trackIterEnd; ++trackIter)
            {
                try
                {
                    Track *pTrack = *trackIter;
                    const LCObjectVec &objectVec = navigate.getRelatedToObjects(*trackIter);

                    // Get reconstructed momentum at dca
                    HelixClass helixFit;
                    helixFit.Initialize_Canonical(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(),
                        pTrack->getTanLambda(), bField);

                    const float recoMomentum(pandora::CartesianVector(helixFit.getMomentum()[0], helixFit.getMomentum()[1],
                        helixFit.getMomentum()[2]).GetMagnitude());

                    // Use momentum magnitude to identify best mc particle
                    MCParticle *pBestMCParticle = NULL;
                    float bestDeltaMomentum(std::numeric_limits<float>::max());

                    for (LCObjectVec::const_iterator itRel = objectVec.begin(), itRelEnd = objectVec.end(); itRel != itRelEnd; ++itRel)
                    {
                        MCParticle *pMCParticle = NULL;
                        pMCParticle = dynamic_cast<MCParticle *>(*itRel);

                        if (NULL == pMCParticle)
                            continue;

                        const float trueMomentum(pandora::CartesianVector(pMCParticle->getMomentum()[0], pMCParticle->getMomentum()[1],
                            pMCParticle->getMomentum()[2]).GetMagnitude());

                        const float deltaMomentum(std::fabs(recoMomentum - trueMomentum));

                        if (deltaMomentum < bestDeltaMomentum)
                        {
                            pBestMCParticle = pMCParticle;
                            bestDeltaMomentum = deltaMomentum;
                        }
                    }

                    if (NULL == pBestMCParticle)
                        continue;

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackToMCParticleRelationship(*pPandora, pTrack,
                        pBestMCParticle));
                }
                catch (StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract track to mc particle relationship: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract track to mc particle relationship, unrecognised exception" << std::endl;
                }
            }
        }
        catch(...)
        {
            streamlog_out(WARNING) << "Failed to extract track to mc particle relationships collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCParticleCreator::CreateCaloHitToMCParticleRelationships(const LCEvent *const pLCEvent) const
{
    static pandora::Pandora *pPandora = PandoraPFANewProcessor::GetPandora();

    typedef std::map<MCParticle *, float> MCParticleToEnergyWeightMap;
    MCParticleToEnergyWeightMap mcParticleToEnergyWeightMap;

    for (StringVector::const_iterator iter = m_settings.m_lcCaloHitRelationCollections.begin(), iterEnd = m_settings.m_lcCaloHitRelationCollections.end();
         iter != iterEnd; ++iter)
    {
        try
        {
            const LCCollection *pMCRelationCollection = pLCEvent->getCollection(*iter);
            LCRelationNavigator navigate(pMCRelationCollection);

            const CalorimeterHitVector &calorimeterHitVector = CaloHitCreator::GetCalorimeterHitVector();

            for (CalorimeterHitVector::const_iterator caloHitIter = calorimeterHitVector.begin(),
                caloHitIterEnd = calorimeterHitVector.end(); caloHitIter != caloHitIterEnd; ++caloHitIter)
            {
                try
                {
                    mcParticleToEnergyWeightMap.clear();
                    const LCObjectVec &objectVec = navigate.getRelatedToObjects(*caloHitIter);

                    for (LCObjectVec::const_iterator itRel = objectVec.begin(), itRelEnd = objectVec.end(); itRel != itRelEnd; ++itRel)
                    {
                        SimCalorimeterHit *pSimHit = dynamic_cast<SimCalorimeterHit *>(*itRel);

                        if (NULL == pSimHit)
                            continue;

                        for (int iCont = 0, iEnd = pSimHit->getNMCContributions(); iCont < iEnd; ++iCont)
                        {
                            mcParticleToEnergyWeightMap[pSimHit->getParticleCont(iCont)] += pSimHit->getEnergyCont(iCont);
                        }
                    }

                    for (MCParticleToEnergyWeightMap::const_iterator mcParticleIter = mcParticleToEnergyWeightMap.begin(),
                        mcParticleIterEnd = mcParticleToEnergyWeightMap.end(); mcParticleIter != mcParticleIterEnd; ++mcParticleIter)
                    {
                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(*pPandora,
                            *caloHitIter, mcParticleIter->first, mcParticleIter->second));
                    }
                }
                catch (StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract calo hit to mc particle relationship: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract calo hit to mc particle relationship, unrecognised exception" << std::endl;
                }
            }
        }
        catch(...)
        {
            streamlog_out(WARNING) << "Failed to extract calo hit to mc particle relationships collection: " << *iter << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}
