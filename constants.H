#pragma once

#include <string>
#include <unordered_map>

/// IMC
static const int NUM_IMC_CHANNELS = 6; // includes channels on all IMCs in a socket. 2 IMCs with each have 3 channels, making up to 6 in total.
static const int NUM_IMC_COUNTERS = 5; // 0-3 are the 4 programmable counters, 4 is the fixed-function DCLK counter

/// general
static const int CACHE_LINE_SIZE = 64;
static const int NUM_SOCKETS = 2;
static const int NUM_CHA_BOXES = 28;

/// base
static const long CHA_MSR_PMON_CTRL_BASE = 0x0E01L;
static const long CHA_MSR_PMON_CTR_BASE = 0x0E08L;
static const int CHA_BASE = 0x10;

/// mesh traffic
// block / data
static const unsigned int LEFT_BL_READ = 0x004003AB;     /// horizontal_bl_ring
static const unsigned int RIGHT_BL_READ = 0x00400CAB;    /// horizontal_bl_ring
static const unsigned int UP_BL_READ = 0x004003AA;       /// vertical_bl_ring
static const unsigned int DOWN_BL_READ = 0x00400CAA;     /// vertical_bl_ring
static const unsigned int VERT_BL_ALL_READ = 0x00400FAA; /// vertical_bl_ring

// addresses / snooping?
static const unsigned int LEFT_AD_READ = 0x004003A7;  /// horizontal_ad_ring
static const unsigned int RIGHT_AD_READ = 0x00400CA7; /// horizontal_ad_ring
static const unsigned int UP_AD_READ = 0x004003A6;    /// vertical_ad_ring
static const unsigned int DOWN_AD_READ = 0x00400CA6;  /// vertical_ad_ring

// acknowledgements
static const unsigned int LEFT_AK_READ = 0x004003A9;  /// horizontal_ak_ring
static const unsigned int RIGHT_AK_READ = 0x00400CA9; /// horizontal_ak_ring
static const unsigned int UP_AK_READ = 0x004003A8;    /// vertical_ak_ring
static const unsigned int DOWN_AK_READ = 0x00400CA8;  /// vertical_ak_ring

// invalidations
static const unsigned int LEFT_IV_READ = 0x004001AD;  /// horizontal_iv_ring
static const unsigned int RIGHT_IV_READ = 0x004004AD; /// horizontal_iv_ring
static const unsigned int UP_IV_READ = 0x004001AC;    /// vertical_iv_ring
static const unsigned int DOWN_IV_READ = 0x004004AC;  /// vertical_iv_ring
// static const unsigned int DOWN_IV_ALL  = 0x004004AC; /// vertical_iv_ring

/// hitme cache
// lookup:  Counts Number of times HitMe Cache is accessed.
static const unsigned int HITME_LOOKUP_READ = 0x0040015E;
static const unsigned int HITME_LOOKUP_WRITE = 0x0040025E;
static const unsigned int HITME_LOOKUP_ALL = 0x0040035E;

// hit:   Counts Number of Hits in HitMe Cache
static const unsigned int HITME_HIT_EX_RDS = 0x0040015F;
static const unsigned int HITME_HIT_SHARED_OWNREQ = 0x0040045F;
static const unsigned int HITME_HIT_WBMTOE = 0x0040085F;
static const unsigned int HITME_HIT_WBMTOI_OR_S = 0x0040105F;
static const unsigned int HITME_HIT_ALL = 0x00401D5F;

// miss:  Counts Number of Misses in HitMe Cache
static const unsigned int HITME_MISS_SHARED_RDINVOWN = 0x00402060;
static const unsigned int HITME_MISS_NOTSHARED_RDINVOWN = 0x00404060;
static const unsigned int HITME_MISS_READ_OR_INV = 0x00408060;
static const unsigned int HITME_MISS_ALL = 0x0040E060;

// update:  Counts the number of Allocate/Update to HitMe Cache
static const unsigned int HITME_UPDATE_DEALLOCATE_RSPFWDI_LOC = 0x00400161;
static const unsigned int HITME_UPDATE_RSPFWDI_REM = 0x00400261;
static const unsigned int HITME_UPDATE_SHARED = 0x00400461;
static const unsigned int HITME_UPDATE_RDINVOWN = 0x00400861;
static const unsigned int HITME_UPDATE_DEALLOCATE = 0x00401061;
static const unsigned int HITME_UPDATE_ALL = 0x00401F61;

///// cache

// llc lookups:  Cache and Snoop Filter Lookups

// Definition: Counts the number of times the LLC was accessed - this includes code,
// data, prefetches and hints coming from L2. This has numerous filters available. Note
// the non-standard filtering equation. This event will count requests that lookup the
// cache multiple times with multiple increments. One must **************ALWAYS************** set umask bit 0
// AND SELECT A STATE OR STATES TO MATCH. --> Otherwise, the event will count nothing. <-- CHA-
// Filter0[24:21,17] bits correspond to [FMESI] state.
// • NOTE: ***Bit 0 of the umask must ***ALWAYS*** be set for this event***. This allows us to match
// against a given state (or states) as programmed in the Cn_MSR_PMON_BOX_FIL-
// TER0.state field bitmask. 0 = I (miss), 4 = S, 5 = E, 6 = M, 7 = F. For example, if
// you wanted to monitor F and S hits, you could set 00001001b in the 8-bit state field.
// To monitor any lookup, set the field to 0x1F. Extra note - it may be a little confusing
// for customers of earlier products. With the CBo and HA functionality combined, it's
// possible to also measure Snoop Filter lookups with bits 1-3 of the FILTER0.state field

/// LLC_LOOKUP may be filtered by the cacheline state (using CHA filter registers).
static const unsigned int LLC_ANY_LOOKUP = 0x00401134;
static const unsigned int LLC_LOCAL_LOOKUP = 0x00403134;
static const unsigned int LLC_REMOTE_LOOKUP = 0x00409134;
static const unsigned int LLC_DATA_READ_LOOKUP = 0x00400334;
static const unsigned int LLC_WRITE_LOOKUP = 0x00400534;
static const unsigned int LLC_REMOTE_SNOOP_LOOKUP = 0x00400934; /// mccalping made these 0x0043, but this implies making reserved bit 1. why???

static const unsigned int DIR_LOOKUP_SNP_NOSNP = 0x00400353;    /// umask is SNP | NO_SNP
static const unsigned int DIR_LOOKUP_SNP = 0x00400153;          /// umask is SNP
static const unsigned int DIR_LOOKUP_NOSNP = 0x00400253;        /// umask is NO_SNP

static const unsigned int XSNP_RESP_EVICT_RSP_HITFSE = 0x00408132;
static const unsigned int SF_EVICTION_S_STATE = 0x0040043d;
static const unsigned int SF_EVICTION_M_STATE = 0x0040013d;
static const unsigned int SF_EVICTION_E_STATE = 0x0040023d;

/// filter
static const unsigned int FILTER0_OFF = 0;
static const unsigned int FILTER0_ALL_LLC = (1 << 24) | (1 << 23) | (1 << 22) | (1 << 21) | (1 << 17); /// for filtering only and all LLC events (FMESI).
static const unsigned int FILTER0_ALL_SF = (1 << 20) | (1 << 19) | (1 << 18);                          /// for filtering only and all SF_EVICTIONs.
static const unsigned int FILTER0_ALL_COMBINED_LLC_SF = FILTER0_ALL_LLC | FILTER0_ALL_SF;              /// for filtering only and all SF_EVICTIONs.
// static const unsigned int FILTER0_ALL_SF = 0x0000003B; /// for filtering only LLC events.
static const unsigned int FILTER1_OFF = 0x0000003B; /// 3B essentially turns off this filter.

static std::unordered_map<unsigned int, std::string> descriptions {
    {                       LEFT_BL_READ,                        "LEFT_BL_READ"},
    {                      RIGHT_BL_READ,                       "RIGHT_BL_READ"},

    {                         UP_BL_READ,                          "UP_BL_READ"},
    {                       DOWN_BL_READ,                        "DOWN_BL_READ"},
    {                   VERT_BL_ALL_READ,                    "VERT_BL_ALL_READ"},

    {                       LEFT_AD_READ,                        "LEFT_AD_READ"},
    {                      RIGHT_AD_READ,                       "RIGHT_AD_READ"},
    {                         UP_AD_READ,                          "UP_AD_READ"},
    {                       DOWN_AD_READ,                        "DOWN_AD_READ"},

    {                       LEFT_AK_READ,                        "LEFT_AK_READ"},
    {                      RIGHT_AK_READ,                       "RIGHT_AK_READ"},
    {                         UP_AK_READ,                          "UP_AK_READ"},
    {                       DOWN_AK_READ,                        "DOWN_AK_READ"},

    {                       LEFT_IV_READ,                        "LEFT_IV_READ"},
    {                      RIGHT_IV_READ,                       "RIGHT_IV_READ"},
    {                         UP_IV_READ,                          "UP_IV_READ"},
    {                       DOWN_IV_READ,                        "DOWN_IV_READ"},

    {                  HITME_LOOKUP_READ,                   "HITME_LOOKUP_READ"},
    {                 HITME_LOOKUP_WRITE,                  "HITME_LOOKUP_WRITE"},
    {                   HITME_LOOKUP_ALL,                    "HITME_LOOKUP_ALL"},

    {                   HITME_HIT_EX_RDS,                    "HITME_HIT_EX_RDS"},
    {            HITME_HIT_SHARED_OWNREQ,             "HITME_HIT_SHARED_OWNREQ"},
    {                   HITME_HIT_WBMTOE,                    "HITME_HIT_WBMTOE"},
    {              HITME_HIT_WBMTOI_OR_S,               "HITME_HIT_WBMTOI_OR_S"},
    {                      HITME_HIT_ALL,                       "HITME_HIT_ALL"},

    {         HITME_MISS_SHARED_RDINVOWN,          "HITME_MISS_SHARED_RDINVOWN"},
    {      HITME_MISS_NOTSHARED_RDINVOWN,       "HITME_MISS_NOTSHARED_RDINVOWN"},
    {             HITME_MISS_READ_OR_INV,              "HITME_MISS_READ_OR_INV"},
    {                     HITME_MISS_ALL,                      "HITME_MISS_ALL"},

    {HITME_UPDATE_DEALLOCATE_RSPFWDI_LOC, "HITME_UPDATE_DEALLOCATE_RSPFWDI_LOC"},
    {           HITME_UPDATE_RSPFWDI_REM,            "HITME_UPDATE_RSPFWDI_REM"},
    {                HITME_UPDATE_SHARED,                 "HITME_UPDATE_SHARED"},
    {              HITME_UPDATE_RDINVOWN,               "HITME_UPDATE_RDINVOWN"},
    {            HITME_UPDATE_DEALLOCATE,             "HITME_UPDATE_DEALLOCATE"},
    {                   HITME_UPDATE_ALL,                    "HITME_UPDATE_ALL"},

    {                     LLC_ANY_LOOKUP,                      "LLC_ANY_LOOKUP"},
    {                   LLC_LOCAL_LOOKUP,                    "LLC_LOCAL_LOOKUP"},
    {                  LLC_REMOTE_LOOKUP,                   "LLC_REMOTE_LOOKUP"},
    {            LLC_REMOTE_SNOOP_LOOKUP,             "LLC_REMOTE_SNOOP_LOOKUP"},
    {               LLC_DATA_READ_LOOKUP,                "LLC_DATA_READ_LOOKUP"},
    {                   LLC_WRITE_LOOKUP,                    "LLC_WRITE_LOOKUP"},

    {               DIR_LOOKUP_SNP_NOSNP,                "DIR_LOOKUP_SNP_NOSNP"},
    {                     DIR_LOOKUP_SNP,                      "DIR_LOOKUP_SNP"},
    {                   DIR_LOOKUP_NOSNP,                    "DIR_LOOKUP_NOSNP"},

    {         XSNP_RESP_EVICT_RSP_HITFSE,          "XSNP_RESP_EVICT_RSP_HITFSE"},
    {                SF_EVICTION_S_STATE,                 "SF_EVICTION_S_STATE"},
    {                SF_EVICTION_M_STATE,                 "SF_EVICTION_M_STATE"},
    {                SF_EVICTION_E_STATE,                 "SF_EVICTION_E_STATE"},

    {                        FILTER0_OFF,                         "FILTER0_OFF"},
    {                    FILTER0_ALL_LLC,                     "FILTER0_ALL_LLC"},
    {                     FILTER0_ALL_SF,                      "FILTER0_ALL_SF"},
    {        FILTER0_ALL_COMBINED_LLC_SF,         "FILTER0_ALL_COMBINED_LLC_SF"},
    {                        FILTER1_OFF,                         "FILTER1_OFF"},
};