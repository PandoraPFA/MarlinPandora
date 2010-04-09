/**
 *  @file   PandoraPFANew/include/TrackCreator.h
 * 
 *  @brief  Header file for the track creator class.
 * 
 *  $Log: $
 */

#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H 1

#include "EVENT/LCEvent.h"
#include "EVENT/Track.h"

#include "HelixClass.h"

#include "Api/PandoraApi.h"

using namespace EVENT;

typedef std::vector<Track *> TrackVector;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TrackCreator class
 */
class TrackCreator
{
public:
    typedef std::vector<std::string> StringVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        StringVector    m_trackCollections;                     ///< The reconstructed track collections
        StringVector    m_v0VertexCollections;                  ///< The v0 vertex collections

        int             m_minTrackHits;                         ///< Track quality cut: the minimum number of track hits
        int             m_maxTrackHits;                         ///< Track quality cut: the maximum number of track hits
        int             m_nHitsForHelixFits;                    ///< The number of hits to be used in helix fits at start/end of tracks

        int             m_useEndTrackHelixForECalProjection;    ///< Use end track fit or full track helix for ECal projection

        float           m_d0TrackCut;                           ///< Track d0 cut used to determine whether track can be used to form pfo
        float           m_z0TrackCut;                           ///< Track z0 cut used to determine whether track can be used to form pfo
        float           m_maxTpcInnerRDistance;                 ///< Track cut on distance from tpc inner r to id whether track can form pfo

        int             m_usingNonVertexTracks;                 ///< Whether can form pfos from tracks that don't start at vertex
        int             m_usingUnmatchedNonVertexTracks;        ///< Whether can form pfos from unmatched tracks that don't start at vertex
        int             m_usingUnmatchedVertexTracks;           ///< Whether can form pfos from unmatched tracks that start at vertex
        float           m_unmatchedVertexTrackMaxEnergy;        ///< Maximum energy for unmatched vertex track

        float           m_d0UnmatchedVertexTrackCut;            ///< d0 cut used to determine whether unmatched vertex track can form pfo
        float           m_z0UnmatchedVertexTrackCut;            ///< z0 cut used to determine whether unmatched vertex track can form pfo
        float           m_zCutForNonVertexTracks;               ///< Non vtx track z cut to determine whether track can be used to form pfo

        int             m_reachesECalNTpcHits;                  ///< Minimum number of tpc hits to consider track as reaching ecal
        float           m_reachesECalTpcOuterDistance;          ///< Max distance from track to tpc r max to id whether track reaches ecal
        float           m_reachesECalTpcZMaxDistance;           ///< Max distance from track to tpc z max to id whether track reaches ecal
        float           m_curvatureToMomentumFactor;            ///< Constant relating track curvature in b field to momentum
    };

    /**
     *  @brief  Create tracks, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */
    StatusCode CreateTracks(const LCEvent *const pLCEvent);

    /**
     *  @brief  Create associations between tracks, V0s, kinks, etc
     * 
     *  @param  pLCEvent the lcio event
     */
    StatusCode CreateTrackAssociations(const LCEvent *const pLCEvent) const;

    /**
     *  @brief  Get the track vector
     * 
     *  @return The track vector
     */
    static const TrackVector &GetTrackVector();

    Settings                m_settings;         ///< The settings

private:
    /**
     *  @brief  Decide whether track reaches the ecal surface
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    void TrackReachesECAL(const Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Perform helix fits to calculate track parameters: momentum at dca, start and end track states
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    void FitTrackHelices(const Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Determine whether a track can be used to form a pfo under the following conditions:
     *          1) if the track proves to be associated with a cluster, OR
     *          2) if the track proves to have no cluster associations
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    void DefineTrackPfoUsage(const Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Project helix to the surface of the ecal
     * 
     *  @param  pHelix helix fit to be projected to ecal surface
     *  @param  referencePoint helix reference point
     *  @param  signPz sign w.r.t. increasing z direction
     */
    pandora::TrackState GetECalProjection(HelixClass *const pHelix, float referencePoint[3], int signPz) const;

    static TrackVector      m_trackVector;      ///< The track vector
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TrackVector &TrackCreator::GetTrackVector()
{
    return m_trackVector;
}

#endif // #ifndef TRACK_CREATOR_H
