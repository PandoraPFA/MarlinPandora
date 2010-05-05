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
typedef std::set<const Track *> TrackList;

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
        StringVector    m_kinkVertexCollections;                ///< The kink vertex collections

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

        int             m_shouldFormTrackRelationships;         ///< Whether to form pandora track relationships using v0 and kink info
    };

    /**
     *  @brief  Create associations between tracks, V0s, kinks, etc
     * 
     *  @param  pLCEvent the lcio event
     */
    StatusCode CreateTrackAssociations(const LCEvent *const pLCEvent);

    /**
     *  @brief  Create tracks, insert user code here
     * 
     *  @param  pLCEvent the lcio event
     */
    StatusCode CreateTracks(const LCEvent *const pLCEvent) const;

    /**
     *  @brief  Get the track vector
     * 
     *  @return The track vector
     */
    static const TrackVector &GetTrackVector();

    /**
     *  @brief  Reset the track creator
     */
    void Reset();

    Settings                m_settings;         ///< The settings

private:
    /**
     *  @brief  Extract kink information from specified lcio collections
     * 
     *  @param  pLCEvent the lcio event
     */
    StatusCode ExtractKinks(const LCEvent *const pLCEvent);

    /**
     *  @brief  Extract v0 information from specified lcio collections
     * 
     *  @param  pLCEvent the lcio event
     */
    StatusCode ExtractV0s(const LCEvent *const pLCEvent);

    /**
     *  @brief  Whether the track vertex conflicts with previously provided relationship information
     * 
     *  @param  trackVec the vector of tracks associated with the vertex
     */
    bool IsConflictingRelationship(const TrackVec &trackVec) const;

    /**
     *  @brief  Whether a track is a v0 track
     * 
     *  @param  pTrack the lcio track
     * 
     *  @return boolean
     */
    bool IsV0(const Track *const pTrack) const;

    /**
     *  @brief  Whether a track is a parent track
     * 
     *  @param  pTrack the lcio track
     * 
     *  @return boolean
     */
    bool IsParent(const Track *const pTrack) const;

    /**
     *  @brief  Whether a track is a daughter track
     * 
     *  @param  pTrack the lcio track
     * 
     *  @return boolean
     */
    bool IsDaughter(const Track *const pTrack) const;

    /**
     *  @brief  Perform helix fits to calculate track parameters: momentum at dca, start and end track states
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    void FitTrackHelices(const Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Decide whether track reaches the ecal surface
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    void TrackReachesECAL(const Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

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

    /**
     *  @brief  Whether track passes the quality cuts required in order to be used to form a pfo
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     *  @param  rInner the track inner radius
     * 
     *  @return boolean
     */
    bool PassesQualityCuts(const Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters, const float rInner) const;

    static TrackVector      m_trackVector;      ///< The track vector

    TrackList               m_v0TrackList;      ///< The list of v0 tracks
    TrackList               m_parentTrackList;  ///< The list of parent tracks
    TrackList               m_daughterTrackList;///< The list of daughter tracks
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TrackVector &TrackCreator::GetTrackVector()
{
    return m_trackVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackCreator::Reset()
{
    m_trackVector.clear();
    m_v0TrackList.clear();
    m_parentTrackList.clear();
    m_daughterTrackList.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackCreator::IsV0(const Track *const pTrack) const
{
    return (m_v0TrackList.end() != m_v0TrackList.find(pTrack));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackCreator::IsParent(const Track *const pTrack) const
{
    return (m_parentTrackList.end() != m_parentTrackList.find(pTrack));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackCreator::IsDaughter(const Track *const pTrack) const
{
    return (m_daughterTrackList.end() != m_daughterTrackList.find(pTrack));
}

#endif // #ifndef TRACK_CREATOR_H
