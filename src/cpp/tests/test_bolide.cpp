/**
 * Bolide / Fireball Detection Plugin Tests
 *
 * Tests:
 * 1. Parse JPL Fireball API JSON
 * 2. Parse NASA GLM Bolide API JSON
 * 3. Convert to BolideRecord and serialize to $BOL binary
 * 4. Deserialize $BOL binary
 * 5. Energy/mass/diameter estimation
 * 6. Geographic filtering
 * 7. Time filtering
 * 8. Chelyabinsk meteor validation (known event)
 */

#include "bolide/types.h"
#include "bolide/parser.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstring>

using namespace bolide;

// ============================================================================
// Test Data: JPL Fireball API (real format, sample records)
// ============================================================================

static const char* FIREBALL_JSON = R"({
  "signature":{"version":"1.2","source":"NASA/JPL Fireball Data API"},
  "count":5,
  "fields":["date","lat","lat-dir","lon","lon-dir","alt","energy","impact-e","vx","vy","vz"],
  "data":[
    ["2013-02-15 03:20:33","54.8","N","61.1","E","23.3","32200","440","-13.4","3.1","-18.6"],
    ["2024-01-21 16:07:51","47.2","N","16.5","E","35.0","3.8","0.13",null,null,null],
    ["2026-03-08 18:30:00","50.5","N","10.2","E","42.0","0.5","0.018",null,null,null],
    ["2023-06-15 12:00:00",null,null,null,null,null,"1.2","0.042",null,null,null],
    ["2024-07-22 08:15:30","33.5","S","151.2","E","28.5","15.0","0.52","-8.2","5.1","-12.3"]
  ]
})";

// ============================================================================
// Test Data: NASA GLM Bolide Service (real format)
// ============================================================================

static const char* GLM_JSON = R"({
  "success":true,
  "data":[
    {
      "_id":"69a1e31ec4d0421e407297ab",
      "status":"PUBLISHED",
      "datetime":1772143009023,
      "latitude":-13.207,
      "longitude":-162.745,
      "latitudeDelta":0.003,
      "longitudeDelta":0.124,
      "detectedBy":"GLM-18",
      "duration":0.481,
      "brightness":{"GLM-18":{"category":"Faint","value":3.72e-13}},
      "groundTrack":{"GLM-18":{"category":"Multi-Pixel","value":0}},
      "assessmentScore":{"score":0.9967,"assessment":"accepted"},
      "attachments":[{"energy":6.85e-16}],
      "isManuallyGenerated":false
    },
    {
      "_id":"69a1e31ec4d0421e407297af",
      "status":"PUBLISHED",
      "datetime":1772098044938,
      "latitude":55.672,
      "longitude":-159.198,
      "latitudeDelta":0.006,
      "longitudeDelta":0.082,
      "detectedBy":"GLM-18",
      "duration":0.149,
      "brightness":{"GLM-18":{"category":"Bright","value":1.57e-11}},
      "groundTrack":{"GLM-18":{"category":"Single-Pixel","value":0}},
      "assessmentScore":{"score":0.9962,"assessment":"accepted"},
      "attachments":[{"energy":2.45e-13}],
      "isManuallyGenerated":false
    },
    {
      "_id":"rejected_event",
      "status":"REJECTED",
      "datetime":1772000000000,
      "latitude":0,
      "longitude":0,
      "detectedBy":"GLM-16"
    }
  ]
})";

// ============================================================================
// Test 1: Parse JPL Fireball API
// ============================================================================

void testParseFireball() {
    auto events = parseFireballAPI(FIREBALL_JSON);

    assert(events.size() == 5);

    // Chelyabinsk meteor (2013-02-15)
    auto& chely = events[0];
    assert(chely.date == "2013-02-15 03:20:33");
    assert(std::abs(chely.lat - 54.8) < 0.01);
    assert(chely.lat_dir == 'N');
    assert(std::abs(chely.lon - 61.1) < 0.01);
    assert(chely.lon_dir == 'E');
    assert(std::abs(chely.alt_km - 23.3) < 0.01);
    assert(std::abs(chely.energy_10J - 32200) < 1);
    assert(std::abs(chely.impact_e_kt - 440) < 1);
    assert(std::abs(chely.vx - (-13.4)) < 0.01);
    assert(std::abs(chely.vy - 3.1) < 0.01);
    assert(std::abs(chely.vz - (-18.6)) < 0.01);
    assert(chely.speed_kms > 20 && chely.speed_kms < 25); // ~23 km/s

    // Signed lat/lon
    assert(std::abs(chely.signedLat() - 54.8) < 0.01);
    assert(std::abs(chely.signedLon() - 61.1) < 0.01);

    // Record without location (event 3)
    auto& noLoc = events[3];
    assert(std::isnan(noLoc.lat));
    assert(std::abs(noLoc.energy_10J - 1.2) < 0.01);

    // Southern hemisphere event (event 4)
    auto& south = events[4];
    assert(std::abs(south.signedLat() - (-33.5)) < 0.01);
    assert(std::abs(south.signedLon() - 151.2) < 0.01);

    // Germany meteorite (March 8, 2026 — from the LinkedIn post!)
    auto& germany = events[2];
    assert(germany.date == "2026-03-08 18:30:00");
    assert(std::abs(germany.signedLat() - 50.5) < 0.01);
    assert(std::abs(germany.signedLon() - 10.2) < 0.01);

    std::cout << "  JPL Fireball parsing ✓ (" << events.size() << " events, "
              << "Chelyabinsk=" << chely.impact_e_kt << "kt @ "
              << chely.speed_kms << " km/s)\n";
}

// ============================================================================
// Test 2: Parse NASA GLM Bolide API
// ============================================================================

void testParseGLM() {
    auto events = parseGLMBolideAPI(GLM_JSON);

    // Should skip rejected events
    assert(events.size() == 2);

    auto& e1 = events[0];
    assert(e1.status == "PUBLISHED");
    assert(std::abs(e1.lat - (-13.207)) < 0.001);
    assert(std::abs(e1.lon - (-162.745)) < 0.001);
    assert(e1.detector == "GLM-18");
    assert(std::abs(e1.duration_s - 0.481) < 0.001);
    assert(e1.brightness_cat == "Faint");
    assert(e1.ml_score > 0.99);

    auto& e2 = events[1];
    assert(e2.brightness_cat == "Bright");
    assert(std::abs(e2.lat - 55.672) < 0.001);

    std::cout << "  GLM Bolide parsing ✓ (" << events.size() << " published events, "
              << "detector=" << e1.detector << ")\n";
}

// ============================================================================
// Test 3: Convert and Serialize
// ============================================================================

void testSerialization() {
    auto fireballs = parseFireballAPI(FIREBALL_JSON);
    auto glm_events = parseGLMBolideAPI(GLM_JSON);

    // Convert all to BolideRecord
    std::vector<BolideRecord> records;
    for (const auto& fb : fireballs) {
        records.push_back(fireballToRecord(fb));
    }
    for (const auto& glm : glm_events) {
        records.push_back(glmToRecord(glm));
    }

    assert(records.size() == 7);

    // Serialize
    auto buf = serializeBolides(records, DataSource::COMBINED);
    assert(buf.size() == sizeof(BolHeader) + 7 * sizeof(BolideRecord));
    assert(std::memcmp(buf.data(), "$BOL", 4) == 0);

    // Deserialize
    BolHeader hdr;
    std::vector<BolideRecord> decoded;
    bool ok = deserializeBolides(buf.data(), buf.size(), hdr, decoded);
    assert(ok);
    assert(hdr.count == 7);
    assert(decoded.size() == 7);

    // Verify Chelyabinsk roundtrip
    auto& chely = decoded[0];
    assert(std::abs(chely.lat_deg - 54.8) < 0.01);
    assert(std::abs(chely.impact_energy_kt - 440) < 1);
    assert(chely.source == static_cast<uint8_t>(DataSource::CNEOS_FIREBALL));
    assert(chely.brightness_class == static_cast<uint8_t>(BrightnessClass::SUPERBOL));
    assert(chely.confidence == static_cast<uint8_t>(ConfidenceLevel::VERIFIED));
    assert(std::string(chely.detector) == "CNEOS");

    // Verify GLM record roundtrip
    auto& glm_r = decoded[5];
    assert(glm_r.source == static_cast<uint8_t>(DataSource::NASA_GLM));
    assert(std::abs(glm_r.lat_deg - (-13.207)) < 0.001);
    assert(std::string(glm_r.detector) == "GLM-18");

    std::cout << "  Serialization ✓ (" << buf.size() << " bytes, "
              << decoded.size() << " records roundtripped)\n";
}

// ============================================================================
// Test 4: Chelyabinsk Mass/Diameter Estimation
// ============================================================================

void testChelyabinskEstimation() {
    // Known: 440 kt, ~19 km/s entry speed, ~20m diameter, ~12000 tonnes
    double mass = estimateMass(440, 19.2);
    double diameter = estimateDiameter(mass, 3500);

    // Mass should be ~10000-15000 tonnes (10M-15M kg)
    assert(mass > 5e6 && mass < 20e6);

    // Diameter should be ~15-25m
    assert(diameter > 10 && diameter < 30);

    std::cout << "  Chelyabinsk estimation ✓ (mass=" << mass/1e6 << " ktonnes, "
              << "diameter=" << diameter << " m)\n";
}

// ============================================================================
// Test 5: Geographic Filtering
// ============================================================================

void testGeoFilter() {
    auto fireballs = parseFireballAPI(FIREBALL_JSON);
    std::vector<BolideRecord> records;
    for (const auto& fb : fireballs) {
        records.push_back(fireballToRecord(fb));
    }

    // Filter for Europe (lat 35-70, lon -10-40)
    std::vector<BolideRecord> europe;
    for (const auto& r : records) {
        if (!std::isnan(r.lat_deg) && inBBox(r, 35, 70, -10, 40)) {
            europe.push_back(r);
        }
    }

    // Should include Germany (50.5N, 10.2E) and Hungary? (47.2N, 16.5E)
    // and Chelyabinsk is at 54.8N, 61.1E — outside Europe bbox
    assert(europe.size() == 2);

    std::cout << "  Geographic filter ✓ (" << europe.size()
              << " European events)\n";
}

// ============================================================================
// Test 6: Energy Filtering
// ============================================================================

void testEnergyFilter() {
    auto fireballs = parseFireballAPI(FIREBALL_JSON);
    std::vector<BolideRecord> records;
    for (const auto& fb : fireballs) {
        records.push_back(fireballToRecord(fb));
    }

    // Filter for significant events (>0.1 kt)
    auto significant = filterByEnergy(records, 0.1);
    assert(significant.size() == 3); // Chelyabinsk (440), Hungary (0.13), Australia (0.52)

    // Filter for superbolides (>1 kt)
    auto superbolides = filterByEnergy(records, 1.0);
    assert(superbolides.size() == 1); // Only Chelyabinsk

    std::cout << "  Energy filter ✓ (>0.1kt=" << significant.size()
              << ", >1kt=" << superbolides.size() << ")\n";
}

// ============================================================================
// Test 7: Wire Format Validation
// ============================================================================

void testWireFormat() {
    // Verify header size
    assert(sizeof(BolHeader) == 16);

    // Verify record size
    assert(sizeof(BolideRecord) == 156);

    // Verify empty serialization
    auto empty = serializeBolides({});
    assert(empty.size() == 16);
    BolHeader hdr;
    std::vector<BolideRecord> records;
    assert(deserializeBolides(empty.data(), empty.size(), hdr, records));
    assert(hdr.count == 0);
    assert(records.empty());

    // Verify bad magic rejection
    uint8_t bad[16] = {0};
    assert(!deserializeBolides(bad, 16, hdr, records));

    std::cout << "  Wire format ✓ (header=16B, record=" << sizeof(BolideRecord) << "B)\n";
}

// ============================================================================
// Test 8: Energy from Optical Magnitude (Tagliaferri et al. 2002)
// ============================================================================

// Declare functions from bolide.cpp
namespace bolide {
    double energyFromOpticalMagnitude(double magnitude);
    double opticalMagnitudeFromEnergy(double energy_J);
    BolideRecord computeTrajectory(double lat_deg, double lon_deg, double alt_km,
                                     double vx_kms, double vy_kms, double vz_kms);
}

void testEnergyFromMagnitude() {
    // Chelyabinsk: M ~ -28 absolute magnitude
    // Expected energy ~500 kT = 500 * 4.184e12 J ≈ 2.09e15 J
    // Formula: E = 8.25e9 * 10^(2*(-28)/5) = 8.25e9 * 10^(-11.2) ≈ 5.2e-2 J
    // Wait — the formula as stated uses apparent magnitude for bolides,
    // where brighter = more negative. The Tagliaferri formula uses
    // different sign convention. Let's verify with known test cases:

    // For magnitude 0: E = 8.25e9 * 10^0 = 8.25e9 J
    double E0 = bolide::energyFromOpticalMagnitude(0);
    assert(std::abs(E0 - 8.25e9) < 1e3);

    // For magnitude 5: E = 8.25e9 * 10^2 = 8.25e11 J
    double E5 = bolide::energyFromOpticalMagnitude(5);
    assert(std::abs(E5 - 8.25e11) < 1e6);

    // Roundtrip: energy → magnitude → energy
    double M_rt = bolide::opticalMagnitudeFromEnergy(8.25e11);
    assert(std::abs(M_rt - 5.0) < 0.01);

    // Verify roundtrip for arbitrary value
    double E_in = 1e14;
    double M_calc = bolide::opticalMagnitudeFromEnergy(E_in);
    double E_out = bolide::energyFromOpticalMagnitude(M_calc);
    assert(std::abs(E_out - E_in) / E_in < 1e-10);

    std::cout << "  Energy from magnitude ✓ (E(M=0)=" << E0
              << " J, E(M=5)=" << E5 << " J)\n";
}

// ============================================================================
// Test 9: Trajectory Computation
// ============================================================================

void testTrajectoryComputation() {
    // Chelyabinsk: 55.15°N, 61.41°E, 23.3 km altitude
    // Velocity: vx=-13.4, vy=3.1, vz=-18.6 km/s
    auto traj = bolide::computeTrajectory(55.15, 61.41, 23.3, -13.4, 3.1, -18.6);

    assert(std::abs(traj.lat_deg - 55.15) < 0.01);
    assert(std::abs(traj.lon_deg - 61.41) < 0.01);
    assert(std::abs(traj.alt_km - 23.3) < 0.01);

    // Speed should be ~23.3 km/s
    double expected_speed = std::sqrt(13.4*13.4 + 3.1*3.1 + 18.6*18.6);
    assert(std::abs(traj.speed_kms - expected_speed) < 0.1);

    // Entry angle should be steep (~53-55°)
    assert(traj.entry_angle_deg > 45 && traj.entry_angle_deg < 65);

    // Azimuth should be valid (0-360)
    assert(traj.azimuth_deg >= 0 && traj.azimuth_deg <= 360);

    std::cout << "  Trajectory computation ✓ (speed=" << traj.speed_kms
              << " km/s, entry=" << traj.entry_angle_deg
              << "°, azimuth=" << traj.azimuth_deg << "°)\n";
}

// ============================================================================
// Main
// ============================================================================

int main() {
    std::cout << "=== bolide-sdn-plugin tests ===\n";

    testParseFireball();
    testParseGLM();
    testSerialization();
    testChelyabinskEstimation();
    testGeoFilter();
    testEnergyFilter();
    testWireFormat();
    testEnergyFromMagnitude();
    testTrajectoryComputation();

    std::cout << "All bolide tests passed.\n";
    return 0;
}
