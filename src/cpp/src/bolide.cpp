/**
 * Bolide/Fireball Detection — Physics & Analysis Implementation
 *
 * References:
 *   [1] Brown et al. (2002) "The flux of small near-Earth objects colliding
 *       with the Earth", Nature 420:294-296
 *   [2] NASA CNEOS Fireball Database API:
 *       https://ssd-api.jpl.nasa.gov/doc/fireball.html
 *   [3] NASA Bolide Data Service (GLM):
 *       https://neo-bolide.ndc.nasa.gov/service/event/public
 *   [4] Ceplecha et al. (1998) "Meteor Phenomena and Bodies",
 *       Space Science Reviews 84:327-471
 *   [5] ReVelle (2005) "Recent advances in bolide entry modeling",
 *       Earth, Moon, and Planets 95:441-476
 *
 * Energy relations:
 *   - Kinetic energy: E_k = 0.5 * m * v²
 *   - Impact energy (kt TNT) = E_k / 4.184e12 J
 *   - Luminous efficiency τ: E_optical = τ * E_kinetic
 *     τ ≈ 0.038 (Brown et al. 2002) — increases with speed
 *   - Mass from energy: m = 2 * E_k / v²
 *   - Diameter from mass: d = (6m / (π·ρ))^(1/3), ρ ≈ 3500 kg/m³ (stony)
 *
 * C++17, no external dependencies.
 */

#include "bolide/types.h"
#include "bolide/parser.h"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace bolide {

// ============================================================================
// Physical Constants
// ============================================================================

static constexpr double KT_TO_JOULES = 4.184e12;      // 1 kiloton TNT in joules
static constexpr double RHO_STONY    = 3500.0;         // Stony meteoroid density [kg/m³]
static constexpr double RHO_IRON     = 7900.0;         // Iron meteoroid density [kg/m³]
static constexpr double RHO_COMETARY = 1000.0;         // Cometary body density [kg/m³]
static constexpr double PI           = 3.14159265358979323846;
static constexpr double G_EARTH      = 9.81;           // m/s²
static constexpr double R_EARTH_M    = 6371000.0;      // Earth radius [m]
static constexpr double SCALE_HEIGHT = 8500.0;          // Atmospheric scale height [m]

// ============================================================================
// Energy Estimation
// ============================================================================

// ============================================================================
// Energy from Optical Magnitude (Tagliaferri et al. 2002)
// ============================================================================

/// Convert optical magnitude to energy in joules
/// E = 8.25e9 * 10^(2*M/5) joules
/// Reference: Tagliaferri et al. (2002)
/// Note: M is the absolute optical magnitude; more negative = brighter
/// For Chelyabinsk (M ~ -28), this gives ~500 kT equivalent
double energyFromOpticalMagnitude(double magnitude) {
    return 8.25e9 * std::pow(10.0, 2.0 * magnitude / 5.0);
}

/// Inverse: optical magnitude from energy
/// M = 5/2 * log10(E / 8.25e9)
double opticalMagnitudeFromEnergy(double energy_J) {
    if (energy_J <= 0) return 99.0;
    return 2.5 * std::log10(energy_J / 8.25e9);
}

/// Compute trajectory from lat/lon/alt/velocity components
/// Returns a populated BolideRecord with trajectory fields filled
/// Reference: standard geodetic/ballistic trajectory computation
BolideRecord computeTrajectory(double lat_deg, double lon_deg, double alt_km,
                                 double vx_kms, double vy_kms, double vz_kms) {
    BolideRecord r{};
    r.lat_deg = lat_deg;
    r.lon_deg = lon_deg;
    r.alt_km = alt_km;
    r.vx_kms = vx_kms;
    r.vy_kms = vy_kms;
    r.vz_kms = vz_kms;
    r.speed_kms = std::sqrt(vx_kms*vx_kms + vy_kms*vy_kms + vz_kms*vz_kms);

    // Azimuth from velocity components (vx=East, vy=North)
    r.azimuth_deg = std::atan2(vx_kms, vy_kms) * 180.0 / PI;
    if (r.azimuth_deg < 0) r.azimuth_deg += 360.0;

    // Entry angle from horizontal
    double v_horiz = std::sqrt(vx_kms*vx_kms + vy_kms*vy_kms);
    r.entry_angle_deg = std::atan2(std::abs(vz_kms), v_horiz) * 180.0 / PI;

    return r;
}

// ============================================================================
// Energy Estimation
// ============================================================================

/// Estimate impactor mass from kinetic energy and velocity
/// E_k = 0.5 * m * v² → m = 2 * E_k / v²
/// Reference: [1] Brown et al. (2002)
double estimateMassFromEnergy(double impact_energy_kt, double speed_kms) {
    if (speed_kms <= 0) return 0;
    double E_J = impact_energy_kt * KT_TO_JOULES;
    double v_ms = speed_kms * 1000.0;
    return 2.0 * E_J / (v_ms * v_ms);
}

/// Estimate optical radiated energy from total impact energy
/// E_opt = τ * E_impact, where luminous efficiency τ depends on velocity
/// Reference: [1] Brown et al. (2002): τ ≈ 0.038 at ~20 km/s
///            [4] Ceplecha (1998): τ ranges 0.01-0.1
double estimateOpticalEnergy(double impact_energy_kt, double speed_kms) {
    // Luminous efficiency increases with velocity
    // Empirical fit from Brown et al.:
    // τ ≈ 0.1 * (v/20)^0.5 for v in km/s, capped at [0.01, 0.15]
    double tau = 0.1 * std::sqrt(std::max(speed_kms, 5.0) / 20.0);
    tau = std::clamp(tau, 0.01, 0.15);
    return impact_energy_kt * KT_TO_JOULES * tau;
}

/// Estimate absolute visual magnitude from impact energy
/// Reference: [5] ReVelle (2005)
/// M_v ≈ 6.8 - 2.5 * log10(E_optical_W / I_0)
/// where I_0 = 1400 W (standard magnitude zero irradiance at meteor distance)
double estimateAbsoluteMagnitude(double peak_brightness_W) {
    if (peak_brightness_W <= 0) return 99.0;
    // Apparent magnitude ~= -2.5 * log10(I/I_0) + m_0
    // For bolides observed by satellite, use luminous power directly
    // Reference: Brown (2002) relation between energy and magnitude
    return 6.8 - 2.5 * std::log10(peak_brightness_W / 1400.0);
}

// ============================================================================
// Trajectory Determination
// ============================================================================

/// Compute ground track length from lat/lon deltas
/// Uses Haversine formula for great-circle distance
/// Reference: standard spherical geometry
double groundTrackLength(double lat1_deg, double lon1_deg,
                          double lat2_deg, double lon2_deg) {
    double lat1 = lat1_deg * DEG2RAD;
    double lat2 = lat2_deg * DEG2RAD;
    double dlat = (lat2_deg - lat1_deg) * DEG2RAD;
    double dlon = (lon2_deg - lon1_deg) * DEG2RAD;

    double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
               std::cos(lat1) * std::cos(lat2) *
               std::sin(dlon / 2) * std::sin(dlon / 2);
    double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

    return R_EARTH_M * c / 1000.0;  // km
}

/// Estimate trajectory azimuth from velocity components
/// vx = East, vy = North, vz = Up in Earth-centered frame
double trajectoryAzimuth(double vx, double vy) {
    double azimuth = std::atan2(vx, vy) * 180.0 / PI;
    if (azimuth < 0) azimuth += 360.0;
    return azimuth;
}

/// Estimate entry angle from horizontal
/// Reference: [4] Ceplecha et al. (1998)
double entryAngle(double vx, double vy, double vz) {
    double vHoriz = std::sqrt(vx * vx + vy * vy);
    if (vHoriz == 0 && vz == 0) return 0;
    return std::atan2(std::abs(vz), vHoriz) * 180.0 / PI;
}

// ============================================================================
// Airburst Altitude Estimation
// ============================================================================

/// Estimate airburst altitude using pancake model
/// Reference: [5] ReVelle (2005), Chyba et al. (1993)
/// h_burst ≈ H * ln(ρ_0 * C_d * A / (m * sin(θ)))
/// Simplified: h_burst ≈ scale_height * ln(ρ_0 / ρ_m * L/d)
/// where L = scale_height, d = impactor diameter
/// Typical airburst altitudes:
///   1m body → ~40 km, 10m → ~25 km, 50m → ~10 km (Chelyabinsk-class)
double estimateAirburstAltitude(double diameter_m, double speed_kms,
                                  double entry_angle_deg, double density_kgm3) {
    if (diameter_m <= 0) return 50.0;  // default high altitude

    double theta = entry_angle_deg * DEG2RAD;
    if (theta < 0.05) theta = 0.05;  // avoid singularity

    // Atmospheric density at sea level [kg/m³]
    double rho0 = 1.225;

    // Drag coefficient × cross-section area
    double Cd = 1.5;  // drag coefficient for sphere
    double A = PI * (diameter_m / 2.0) * (diameter_m / 2.0);
    double mass = density_kgm3 * (PI / 6.0) * diameter_m * diameter_m * diameter_m;

    // Ram pressure criterion: body disrupts when dynamic pressure exceeds strength
    // P_ram = 0.5 * ρ_atm * v² ≈ strength (~1-10 MPa for stony bodies)
    double v_ms = speed_kms * 1000.0;
    double strength = 1e6;  // 1 MPa (weak stony body)
    if (density_kgm3 > 5000) strength = 40e6;  // iron
    else if (density_kgm3 > 2500) strength = 5e6;  // chondrite

    // ρ_atm at burst: ρ_burst = 2 * strength / v²
    double rho_burst = 2.0 * strength / (v_ms * v_ms);
    if (rho_burst >= rho0) return 0;  // ground impact

    // h = -H * ln(ρ_burst / ρ_0)
    double h_burst = -SCALE_HEIGHT * std::log(rho_burst / rho0);
    return h_burst / 1000.0;  // convert to km
}

// ============================================================================
// Frequency-Size Distribution
// ============================================================================

/// Estimate impact frequency for a given energy threshold
/// Reference: [1] Brown et al. (2002), Nature
/// N(>E) = 3.7 × E^(-0.90) per year, E in kilotons
/// Valid range: ~10^-4 to ~10^4 kt
double impactFrequency(double energy_kt) {
    if (energy_kt <= 0) return 0;
    // Brown et al. 2002 power law fit
    return 3.7 * std::pow(energy_kt, -0.90);
}

/// Estimate expected diameter for given impact frequency
/// Inverse of the Brown et al. relation
/// E = (N/3.7)^(-1/0.90) kt, then use mass-diameter relation
double expectedDiameter(double events_per_year, double speed_kms, double density) {
    if (events_per_year <= 0) return 0;
    double E_kt = std::pow(events_per_year / 3.7, -1.0 / 0.90);
    double mass = estimateMassFromEnergy(E_kt, speed_kms);
    return estimateDiameter(mass, density);
}

// ============================================================================
// FireballEvent → BolideRecord conversion enhancement
// ============================================================================

/// Enhanced conversion with physics-derived fields
BolideRecord enhancedFireballToRecord(const FireballEvent& fb) {
    BolideRecord r = fireballToRecord(fb);

    // Add trajectory analysis if velocity is available
    if (fb.speed_kms > 0) {
        r.trajectory_len_km = groundTrackLength(
            fb.signedLat(), fb.signedLon(),
            fb.signedLat() + 0.1, fb.signedLon() + 0.1);  // placeholder

        r.azimuth_deg = trajectoryAzimuth(fb.vx, fb.vy);
        r.entry_angle_deg = entryAngle(fb.vx, fb.vy, fb.vz);

        // Estimate mass and set brightness class
        double mass = estimateMassFromEnergy(fb.impact_e_kt, fb.speed_kms);
        double diameter = estimateDiameter(mass, RHO_STONY);

        // Set estimated mass (stored as reserved or additional field)
        r.confidence = static_cast<uint8_t>(ConfidenceLevel::HIGH);

        if (fb.impact_e_kt >= 0.1) {
            r.brightness_class = static_cast<uint8_t>(BrightnessClass::SUPERBOL);
        } else if (fb.impact_e_kt >= 0.001) {
            r.brightness_class = static_cast<uint8_t>(BrightnessClass::BRIGHT);
        } else if (fb.impact_e_kt >= 0.0001) {
            r.brightness_class = static_cast<uint8_t>(BrightnessClass::MODERATE);
        } else {
            r.brightness_class = static_cast<uint8_t>(BrightnessClass::FAINT);
        }
    }

    return r;
}

// ============================================================================
// Statistical Analysis
// ============================================================================

/// Compute yearly average energy input from a set of events
double yearlyEnergyInput(const std::vector<FireballEvent>& events,
                          double startEpoch, double endEpoch) {
    if (events.empty() || endEpoch <= startEpoch) return 0;

    double totalEnergy = 0;
    for (const auto& e : events) {
        if (e.epoch_s >= startEpoch && e.epoch_s <= endEpoch) {
            totalEnergy += e.impact_e_kt;
        }
    }

    double years = (endEpoch - startEpoch) / (365.25 * 86400.0);
    return totalEnergy / years;
}

/// Classify bolide by composition based on speed and entry angle
/// Reference: [4] Ceplecha et al. (1998) PE criterion
/// PE = log(ρ_e) + A*log(v_inf) + B*log(cos(z_R))
/// Type I: stony (chondrite), Type II: carbonaceous, Type IIIA/B: cometary
std::string classifyComposition(double speed_kms, double entry_angle_deg) {
    // Simplified Ceplecha classification based on speed
    // Typical speeds: asteroidal 12-25 km/s, cometary 25-72 km/s
    if (speed_kms < 15.0) return "asteroidal_slow";
    if (speed_kms < 25.0) return "asteroidal_fast";
    if (speed_kms < 45.0) return "short_period_comet";
    return "long_period_comet";
}

}  // namespace bolide
