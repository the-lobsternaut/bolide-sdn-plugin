#ifndef BOLIDE_TYPES_H
#define BOLIDE_TYPES_H

/**
 * Bolide / Fireball Detection Plugin Types
 * ==========================================
 *
 * Consumes data from two NASA sources:
 *
 * 1. CNEOS Fireball API (JPL)
 *    https://ssd-api.jpl.nasa.gov/fireball.api
 *    Historical fireball events: date, location, altitude, energy,
 *    impact energy, entry velocity components
 *
 * 2. NASA Bolide Detection Service (GLM)
 *    https://neo-bolide.ndc.nasa.gov/service/event/public
 *    Near-real-time bolide detections from Geostationary Lightning
 *    Mapper sensors (GLM-16, GLM-17, GLM-18, EUMETSAT LI)
 *    Includes: lat/lon, duration, brightness, energy, ground track,
 *    NetCDF source files, trajectory charts, confidence scores
 *
 * Output: $BOL FlatBuffer-aligned binary records
 *
 * File identifier: $BOL
 * Wire format:
 *   Header (16 bytes): magic[4]="$BOL", version(u32), source(u32), count(u32)
 *   N × BolideRecord (176 bytes each)
 *
 * Source enum: 0=CNEOS_FIREBALL, 1=NASA_GLM, 2=EUMETSAT_LI, 3=COMBINED
 */

#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <array>

namespace bolide {

// ============================================================================
// Constants
// ============================================================================

static constexpr char BOL_FILE_ID[4] = {'$', 'B', 'O', 'L'};
static constexpr uint32_t BOL_VERSION = 1;

static constexpr double DEG2RAD = M_PI / 180.0;
static constexpr double R_EARTH = 6371000.0;  // m

// ============================================================================
// Enums
// ============================================================================

enum class DataSource : uint32_t {
    CNEOS_FIREBALL = 0,   // JPL fireball database
    NASA_GLM       = 1,   // NOAA GLM-16/17/18 via NASA Bolide Service
    EUMETSAT_LI    = 2,   // EUMETSAT Lightning Imager (MTG-I)
    COMBINED       = 3,   // Merged from multiple sources
};

enum class BrightnessClass : uint8_t {
    UNKNOWN   = 0,
    FAINT     = 1,   // barely detectable
    MODERATE  = 2,   // clearly visible
    BRIGHT    = 3,   // saturating detector
    SUPERBOL  = 4,   // >0.1 kt impact energy
};

enum class GroundTrackType : uint8_t {
    UNKNOWN      = 0,
    SINGLE_PIXEL = 1,   // point detection
    MULTI_PIXEL  = 2,   // resolved ground track
    EXTENDED     = 3,   // long trajectory (>100km)
};

enum class ConfidenceLevel : uint8_t {
    UNKNOWN  = 0,
    LOW      = 1,   // may be lightning or artifact
    MEDIUM   = 2,   // likely bolide
    HIGH     = 3,   // confirmed bolide
    VERIFIED = 4,   // cross-validated with multiple sensors
};

// ============================================================================
// Core Record Types
// ============================================================================

#pragma pack(push, 1)

/// Wire format header (16 bytes)
struct BolHeader {
    char     magic[4];    // "$BOL"
    uint32_t version;     // BOL_VERSION
    uint32_t source;      // DataSource enum
    uint32_t count;       // number of records
};
static_assert(sizeof(BolHeader) == 16, "BolHeader must be 16 bytes");

/// Bolide event record (176 bytes, packed)
struct BolideRecord {
    // Time
    double   epoch_s;          // Unix timestamp of peak brightness [s]

    // Location at peak brightness
    double   lat_deg;          // latitude [deg] (negative = South)
    double   lon_deg;          // longitude [deg] (negative = West)
    double   alt_km;           // altitude above geoid [km] (NaN if unknown)

    // Entry velocity (Earth-centered, km/s) — from CNEOS
    double   vx_kms;           // X component [km/s]
    double   vy_kms;           // Y component [km/s]
    double   vz_kms;           // Z component [km/s]
    double   speed_kms;        // total speed [km/s] (derived)

    // Energy
    double   radiated_energy_J; // total optical radiated energy [J]
    double   impact_energy_kt;  // estimated total impact energy [kilotons TNT]
    double   brightness_W;      // peak brightness [W] (from GLM)

    // Ground track extent
    double   lat_delta_deg;    // latitude extent [deg]
    double   lon_delta_deg;    // longitude extent [deg]
    double   duration_s;       // event duration [s]

    // Trajectory (if resolved)
    double   trajectory_len_km; // ground track length [km]
    double   azimuth_deg;       // trajectory azimuth [deg from N]
    double   entry_angle_deg;   // entry angle from horizontal [deg]

    // Classification
    uint8_t  source;            // DataSource enum
    uint8_t  brightness_class;  // BrightnessClass enum
    uint8_t  ground_track;      // GroundTrackType enum
    uint8_t  confidence;        // ConfidenceLevel enum

    // Assessment
    float    ml_score;          // ML classification score [0-1]

    // Sensor
    char     detector[8];       // e.g., "GLM-16", "GLM-18", "LI"

    uint8_t  reserved[4];
};
static_assert(sizeof(BolideRecord) == 156, "BolideRecord size mismatch");

#pragma pack(pop)

// ============================================================================
// Parsed Event (in-memory, from JSON API responses)
// ============================================================================

struct FireballEvent {
    // CNEOS fields
    std::string date;          // "YYYY-MM-DD hh:mm:ss"
    double lat = NAN, lon = NAN;
    char   lat_dir = ' ', lon_dir = ' ';
    double alt_km = NAN;
    double energy_10J = 0;     // ×10^10 joules
    double impact_e_kt = 0;    // kilotons
    double vx = 0, vy = 0, vz = 0; // km/s

    // Derived
    double epoch_s = 0;
    double speed_kms = 0;

    /// Convert lat/lon directions to signed degrees
    double signedLat() const {
        return (lat_dir == 'S') ? -std::abs(lat) : std::abs(lat);
    }
    double signedLon() const {
        return (lon_dir == 'W') ? -std::abs(lon) : std::abs(lon);
    }
};

struct GLMBolideEvent {
    std::string id;
    std::string status;         // "PUBLISHED", "SUBMITTED", etc.
    double epoch_s = 0;         // Unix timestamp ms / 1000
    double lat = 0, lon = 0;
    double lat_delta = 0, lon_delta = 0;
    double duration_s = 0;
    double brightness_W = 0;
    double energy_J = 0;
    std::string detector;       // "GLM-16", "GLM-18", "LI"
    std::string brightness_cat; // "Faint", "Bright", etc.
    std::string ground_track_cat; // "Single-Pixel", "Multi-Pixel"
    float  ml_score = 0;
    bool   is_manual = false;
};

// ============================================================================
// Serialization Helpers
// ============================================================================

inline std::vector<uint8_t> serializeBolides(const std::vector<BolideRecord>& records,
                                               DataSource source = DataSource::COMBINED) {
    size_t size = sizeof(BolHeader) + records.size() * sizeof(BolideRecord);
    std::vector<uint8_t> buf(size);

    BolHeader hdr;
    std::memcpy(hdr.magic, BOL_FILE_ID, 4);
    hdr.version = BOL_VERSION;
    hdr.source = static_cast<uint32_t>(source);
    hdr.count = static_cast<uint32_t>(records.size());
    std::memcpy(buf.data(), &hdr, sizeof(BolHeader));

    if (!records.empty()) {
        std::memcpy(buf.data() + sizeof(BolHeader),
                    records.data(),
                    records.size() * sizeof(BolideRecord));
    }
    return buf;
}

inline bool deserializeBolides(const uint8_t* data, size_t len,
                                BolHeader& hdr,
                                std::vector<BolideRecord>& records) {
    if (len < sizeof(BolHeader)) return false;
    std::memcpy(&hdr, data, sizeof(BolHeader));
    if (std::memcmp(hdr.magic, BOL_FILE_ID, 4) != 0) return false;

    size_t expected = sizeof(BolHeader) + hdr.count * sizeof(BolideRecord);
    if (len < expected) return false;

    records.resize(hdr.count);
    if (hdr.count > 0) {
        std::memcpy(records.data(),
                    data + sizeof(BolHeader),
                    hdr.count * sizeof(BolideRecord));
    }
    return true;
}

// ============================================================================
// Conversion: FireballEvent → BolideRecord
// ============================================================================

inline BolideRecord fireballToRecord(const FireballEvent& fb) {
    BolideRecord r{};
    r.epoch_s = fb.epoch_s;
    r.lat_deg = fb.signedLat();
    r.lon_deg = fb.signedLon();
    r.alt_km = fb.alt_km;
    r.vx_kms = fb.vx;
    r.vy_kms = fb.vy;
    r.vz_kms = fb.vz;
    r.speed_kms = std::sqrt(fb.vx*fb.vx + fb.vy*fb.vy + fb.vz*fb.vz);
    r.radiated_energy_J = fb.energy_10J * 1e10;
    r.impact_energy_kt = fb.impact_e_kt;
    r.brightness_W = 0; // not available from CNEOS
    r.lat_delta_deg = 0;
    r.lon_delta_deg = 0;
    r.duration_s = 0;
    r.trajectory_len_km = 0;
    r.azimuth_deg = NAN;
    r.entry_angle_deg = NAN;
    r.source = static_cast<uint8_t>(DataSource::CNEOS_FIREBALL);

    // Classify by impact energy
    if (fb.impact_e_kt >= 0.1)
        r.brightness_class = static_cast<uint8_t>(BrightnessClass::SUPERBOL);
    else if (fb.impact_e_kt >= 0.01)
        r.brightness_class = static_cast<uint8_t>(BrightnessClass::BRIGHT);
    else if (fb.impact_e_kt >= 0.001)
        r.brightness_class = static_cast<uint8_t>(BrightnessClass::MODERATE);
    else
        r.brightness_class = static_cast<uint8_t>(BrightnessClass::FAINT);

    r.ground_track = static_cast<uint8_t>(GroundTrackType::UNKNOWN);
    r.confidence = static_cast<uint8_t>(ConfidenceLevel::VERIFIED); // CNEOS = verified
    r.ml_score = 1.0f;
    std::strncpy(r.detector, "CNEOS", 7);
    r.detector[7] = '\0';
    return r;
}

inline BolideRecord glmToRecord(const GLMBolideEvent& glm) {
    BolideRecord r{};
    r.epoch_s = glm.epoch_s;
    r.lat_deg = glm.lat;
    r.lon_deg = glm.lon;
    r.alt_km = NAN; // GLM doesn't report altitude
    r.vx_kms = 0; r.vy_kms = 0; r.vz_kms = 0;
    r.speed_kms = 0;
    r.radiated_energy_J = glm.energy_J;
    r.impact_energy_kt = 0; // can be estimated
    r.brightness_W = glm.brightness_W;
    r.lat_delta_deg = glm.lat_delta;
    r.lon_delta_deg = glm.lon_delta;
    r.duration_s = glm.duration_s;

    // Estimate ground track length from lat/lon delta
    double dlat_m = glm.lat_delta * DEG2RAD * R_EARTH;
    double dlon_m = glm.lon_delta * DEG2RAD * R_EARTH * std::cos(glm.lat * DEG2RAD);
    r.trajectory_len_km = std::sqrt(dlat_m*dlat_m + dlon_m*dlon_m) / 1000.0;

    r.azimuth_deg = (r.trajectory_len_km > 0.1)
                    ? std::atan2(dlon_m, dlat_m) / DEG2RAD : NAN;
    r.entry_angle_deg = NAN;

    r.source = static_cast<uint8_t>(
        (glm.detector.find("LI") != std::string::npos)
        ? DataSource::EUMETSAT_LI : DataSource::NASA_GLM);

    // Brightness classification
    if (glm.brightness_cat == "Bright")
        r.brightness_class = static_cast<uint8_t>(BrightnessClass::BRIGHT);
    else if (glm.brightness_cat == "Moderate")
        r.brightness_class = static_cast<uint8_t>(BrightnessClass::MODERATE);
    else
        r.brightness_class = static_cast<uint8_t>(BrightnessClass::FAINT);

    // Ground track
    if (glm.ground_track_cat == "Multi-Pixel")
        r.ground_track = static_cast<uint8_t>(GroundTrackType::MULTI_PIXEL);
    else if (glm.ground_track_cat == "Single-Pixel")
        r.ground_track = static_cast<uint8_t>(GroundTrackType::SINGLE_PIXEL);
    else
        r.ground_track = static_cast<uint8_t>(GroundTrackType::UNKNOWN);

    // Confidence from ML score
    if (glm.ml_score >= 0.99)
        r.confidence = static_cast<uint8_t>(ConfidenceLevel::HIGH);
    else if (glm.ml_score >= 0.95)
        r.confidence = static_cast<uint8_t>(ConfidenceLevel::MEDIUM);
    else
        r.confidence = static_cast<uint8_t>(ConfidenceLevel::LOW);

    r.ml_score = glm.ml_score;
    std::strncpy(r.detector, glm.detector.c_str(), 7);
    r.detector[7] = '\0';
    return r;
}

// ============================================================================
// Analysis Helpers
// ============================================================================

/// Estimate pre-atmospheric mass from impact energy (kg)
/// E_kt ≈ 0.5 * m * v^2 / (4.184e12)  →  m = 2 * E_kt * 4.184e12 / v^2
inline double estimateMass(double impact_e_kt, double speed_kms) {
    if (speed_kms <= 0 || impact_e_kt <= 0) return 0;
    double v_ms = speed_kms * 1000.0;
    return 2.0 * impact_e_kt * 4.184e12 / (v_ms * v_ms);
}

/// Estimate diameter from mass assuming spherical stony body (ρ = 3500 kg/m³)
inline double estimateDiameter(double mass_kg, double density = 3500.0) {
    if (mass_kg <= 0) return 0;
    double volume = mass_kg / density;
    return 2.0 * std::cbrt(3.0 * volume / (4.0 * M_PI));
}

/// Check if event is within a geographic bounding box
inline bool inBBox(const BolideRecord& r,
                    double latMin, double latMax,
                    double lonMin, double lonMax) {
    return r.lat_deg >= latMin && r.lat_deg <= latMax &&
           r.lon_deg >= lonMin && r.lon_deg <= lonMax;
}

/// Filter records by minimum impact energy
inline std::vector<BolideRecord> filterByEnergy(
    const std::vector<BolideRecord>& records, double minEnergy_kt) {
    std::vector<BolideRecord> out;
    for (const auto& r : records) {
        if (r.impact_energy_kt >= minEnergy_kt) out.push_back(r);
    }
    return out;
}

/// Filter records by time range
inline std::vector<BolideRecord> filterByTime(
    const std::vector<BolideRecord>& records,
    double epochMin, double epochMax) {
    std::vector<BolideRecord> out;
    for (const auto& r : records) {
        if (r.epoch_s >= epochMin && r.epoch_s <= epochMax) out.push_back(r);
    }
    return out;
}

}  // namespace bolide

#endif  // BOLIDE_TYPES_H
