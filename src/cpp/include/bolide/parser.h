#ifndef BOLIDE_PARSER_H
#define BOLIDE_PARSER_H

/**
 * JSON Parsers for NASA Bolide Data Sources
 *
 * Minimal JSON parsing — no external dependencies.
 * Handles the specific JSON formats from:
 *   1. JPL Fireball API: https://ssd-api.jpl.nasa.gov/fireball.api
 *   2. NASA Bolide Service: https://neo-bolide.ndc.nasa.gov/service/event/public
 *
 * Both return well-structured JSON. We parse with simple string scanning
 * since the formats are stable and we can't use rapidjson/nlohmann.
 */

#include "types.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>

namespace bolide {

// ============================================================================
// Utility: Minimal JSON Value Extraction
// ============================================================================

namespace json {

/// Find a key in JSON and return the value string (handles strings, numbers, null)
inline std::string findValue(const std::string& json, const std::string& key,
                              size_t startPos = 0) {
    std::string needle = "\"" + key + "\"";
    size_t pos = json.find(needle, startPos);
    if (pos == std::string::npos) return "";

    // Skip past key and colon
    pos = json.find(':', pos + needle.size());
    if (pos == std::string::npos) return "";
    pos++;

    // Skip whitespace
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t' ||
           json[pos] == '\n' || json[pos] == '\r')) pos++;

    if (pos >= json.size()) return "";

    // String value
    if (json[pos] == '"') {
        size_t end = json.find('"', pos + 1);
        if (end == std::string::npos) return "";
        return json.substr(pos + 1, end - pos - 1);
    }

    // null
    if (json.substr(pos, 4) == "null") return "";

    // Number or boolean
    size_t end = pos;
    while (end < json.size() && json[end] != ',' && json[end] != '}' &&
           json[end] != ']' && json[end] != ' ' && json[end] != '\n') end++;
    return json.substr(pos, end - pos);
}

/// Find numeric value, returns NaN if not found
inline double findDouble(const std::string& json, const std::string& key,
                          size_t startPos = 0) {
    std::string val = findValue(json, key, startPos);
    if (val.empty()) return NAN;
    return std::strtod(val.c_str(), nullptr);
}

/// Find the next array element start position (opening '{' or '[')
inline size_t findNextObject(const std::string& json, size_t pos) {
    return json.find('{', pos);
}

/// Find the matching closing brace
inline size_t findMatchingBrace(const std::string& json, size_t openPos) {
    int depth = 0;
    for (size_t i = openPos; i < json.size(); i++) {
        if (json[i] == '{') depth++;
        else if (json[i] == '}') {
            depth--;
            if (depth == 0) return i;
        }
    }
    return std::string::npos;
}

/// Parse ISO date "YYYY-MM-DD hh:mm:ss" to Unix epoch
inline double parseDateTime(const std::string& dt) {
    if (dt.size() < 19) return 0;
    struct tm tm = {};
    tm.tm_year = std::stoi(dt.substr(0, 4)) - 1900;
    tm.tm_mon  = std::stoi(dt.substr(5, 2)) - 1;
    tm.tm_mday = std::stoi(dt.substr(8, 2));
    tm.tm_hour = std::stoi(dt.substr(11, 2));
    tm.tm_min  = std::stoi(dt.substr(14, 2));
    tm.tm_sec  = std::stoi(dt.substr(17, 2));
    // timegm for UTC (POSIX)
    return static_cast<double>(timegm(&tm));
}

}  // namespace json

// ============================================================================
// Parse JPL Fireball API Response
// ============================================================================

/**
 * Parse JSON from https://ssd-api.jpl.nasa.gov/fireball.api
 *
 * Format:
 * {
 *   "count": N,
 *   "fields": ["date","lat","lat-dir","lon","lon-dir","alt","energy","impact-e",
 *              "vx","vy","vz"],
 *   "data": [
 *     ["2015-10-13 12:23:08","8.0","S","52.5","W",null,"2.3","0.082",null,null,null],
 *     ...
 *   ]
 * }
 */
inline std::vector<FireballEvent> parseFireballAPI(const std::string& jsonStr) {
    std::vector<FireballEvent> events;

    // Find "data" array
    size_t dataPos = jsonStr.find("\"data\"");
    if (dataPos == std::string::npos) return events;

    // Find opening bracket
    size_t arrStart = jsonStr.find('[', dataPos);
    if (arrStart == std::string::npos) return events;

    // Parse each row (array of arrays)
    size_t pos = arrStart + 1;
    while (pos < jsonStr.size()) {
        // Find next inner array
        size_t rowStart = jsonStr.find('[', pos);
        if (rowStart == std::string::npos) break;

        size_t rowEnd = jsonStr.find(']', rowStart);
        if (rowEnd == std::string::npos) break;

        std::string row = jsonStr.substr(rowStart + 1, rowEnd - rowStart - 1);

        // Split by comma, handling quoted strings and null
        std::vector<std::string> fields;
        size_t fieldStart = 0;
        bool inQuote = false;
        for (size_t i = 0; i <= row.size(); i++) {
            if (i == row.size() || (row[i] == ',' && !inQuote)) {
                std::string field = row.substr(fieldStart, i - fieldStart);
                // Trim whitespace and quotes
                while (!field.empty() && (field.front() == ' ' || field.front() == '"'))
                    field.erase(field.begin());
                while (!field.empty() && (field.back() == ' ' || field.back() == '"'))
                    field.pop_back();
                if (field == "null") field = "";
                fields.push_back(field);
                fieldStart = i + 1;
            }
            if (i < row.size() && row[i] == '"') inQuote = !inQuote;
        }

        // Map fields: date, lat, lat-dir, lon, lon-dir, alt, energy, impact-e, vx, vy, vz
        if (fields.size() >= 8) {
            FireballEvent fb;
            fb.date = fields[0];
            fb.epoch_s = json::parseDateTime(fb.date);

            if (!fields[1].empty()) fb.lat = std::strtod(fields[1].c_str(), nullptr);
            if (!fields[2].empty()) fb.lat_dir = fields[2][0];
            if (!fields[3].empty()) fb.lon = std::strtod(fields[3].c_str(), nullptr);
            if (!fields[4].empty()) fb.lon_dir = fields[4][0];
            if (!fields[5].empty()) fb.alt_km = std::strtod(fields[5].c_str(), nullptr);
            if (!fields[6].empty()) fb.energy_10J = std::strtod(fields[6].c_str(), nullptr);
            if (!fields[7].empty()) fb.impact_e_kt = std::strtod(fields[7].c_str(), nullptr);

            if (fields.size() >= 11) {
                if (!fields[8].empty()) fb.vx = std::strtod(fields[8].c_str(), nullptr);
                if (!fields[9].empty()) fb.vy = std::strtod(fields[9].c_str(), nullptr);
                if (!fields[10].empty()) fb.vz = std::strtod(fields[10].c_str(), nullptr);
            }
            fb.speed_kms = std::sqrt(fb.vx*fb.vx + fb.vy*fb.vy + fb.vz*fb.vz);

            events.push_back(fb);
        }

        pos = rowEnd + 1;
    }

    return events;
}

// ============================================================================
// Parse NASA Bolide Service (GLM) Response
// ============================================================================

/**
 * Parse JSON from https://neo-bolide.ndc.nasa.gov/service/event/public
 *
 * Format:
 * {
 *   "success": true,
 *   "data": [
 *     {
 *       "_id": "...",
 *       "status": "PUBLISHED",
 *       "datetime": 1772143009023,  (epoch ms)
 *       "latitude": -13.207,
 *       "longitude": -162.745,
 *       "latitudeDelta": 0.003,
 *       "longitudeDelta": 0.124,
 *       "detectedBy": "GLM-18",
 *       "duration": 0.481,
 *       "brightness": {"GLM-18": {"category": "Faint", "value": 3.72e-13}},
 *       "assessmentScore": {"score": 0.997, "assessment": "accepted"},
 *       ...
 *     }
 *   ]
 * }
 */
inline std::vector<GLMBolideEvent> parseGLMBolideAPI(const std::string& jsonStr) {
    std::vector<GLMBolideEvent> events;

    // Find "data" array
    size_t dataPos = jsonStr.find("\"data\"");
    if (dataPos == std::string::npos) return events;

    // Iterate over objects in the data array
    size_t pos = dataPos;
    while (true) {
        size_t objStart = json::findNextObject(jsonStr, pos + 1);
        if (objStart == std::string::npos) break;

        size_t objEnd = json::findMatchingBrace(jsonStr, objStart);
        if (objEnd == std::string::npos) break;

        std::string obj = jsonStr.substr(objStart, objEnd - objStart + 1);

        GLMBolideEvent evt;
        evt.id = json::findValue(obj, "_id");
        evt.status = json::findValue(obj, "status");

        // Only process published events
        if (evt.status != "PUBLISHED") {
            pos = objEnd;
            continue;
        }

        double datetime_ms = json::findDouble(obj, "datetime");
        evt.epoch_s = datetime_ms / 1000.0;

        evt.lat = json::findDouble(obj, "latitude");
        evt.lon = json::findDouble(obj, "longitude");
        evt.lat_delta = json::findDouble(obj, "latitudeDelta");
        evt.lon_delta = json::findDouble(obj, "longitudeDelta");
        evt.duration_s = json::findDouble(obj, "duration");
        evt.detector = json::findValue(obj, "detectedBy");

        // Parse brightness (nested object)
        size_t brightnessPos = obj.find("\"brightness\"");
        if (brightnessPos != std::string::npos) {
            evt.brightness_cat = json::findValue(obj, "category", brightnessPos);
            evt.brightness_W = json::findDouble(obj, "value", brightnessPos);
        }

        // Parse assessment score
        evt.ml_score = static_cast<float>(json::findDouble(obj, "score"));

        // Parse ground track
        size_t gtPos = obj.find("\"groundTrack\"");
        if (gtPos != std::string::npos) {
            evt.ground_track_cat = json::findValue(obj, "category", gtPos);
        }

        // Parse energy from attachments
        size_t attachPos = obj.find("\"attachments\"");
        if (attachPos != std::string::npos) {
            evt.energy_J = json::findDouble(obj, "energy", attachPos);
        }

        evt.is_manual = (json::findValue(obj, "isManuallyGenerated") == "true");

        events.push_back(evt);
        pos = objEnd;
    }

    return events;
}

}  // namespace bolide

#endif  // BOLIDE_PARSER_H
