#pragma once
#include "diagnostics.h"
#include "xmlParser.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <set>

namespace input_config {
inline std::string extension(const std::string &path) {
    auto ext = std::filesystem::path(path).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return std::tolower(c); });
    return ext;
}
[[noreturn]] inline void fail(const std::string &path, const std::string &message) {
    throw diagnostics::Error("INPUT.CONFIG", "Input '" + path + "': " + message);
}
inline nlohmann::json read_json(const std::string &path) {
    std::ifstream file(path);
    if (!file) fail(path, "cannot open file.");
    try {
        auto result = nlohmann::json::parse(file);
        if (!result.is_object()) fail(path, "expected a JSON object.");
        return result;
    } catch (const nlohmann::json::exception &e) {
        fail(path, std::string("invalid JSON: ") + e.what());
    }
}
inline double number(const nlohmann::json &v, const std::string &path, const std::string &field) {
    if (!v.is_number() || !std::isfinite(v.get<double>())) fail(path, field + " must be a finite number.");
    return v.get<double>();
}
inline void keys(const nlohmann::json &v, const std::set<std::string> &allowed,
                 const std::string &path) {
    if (!v.is_object()) fail(path, "expected an object.");
    for (auto it = v.begin(); it != v.end(); ++it)
        if (!allowed.count(it.key())) fail(path, "unknown field '" + it.key() + "'.");
}
// Explicit path wins. Without one, preserve the historical XML default;
// fall back to the JSON sibling when XML is absent.
inline std::string settings_path(int argc, char **argv, const std::string &stem) {
    std::string path;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--settings") {
            if (!path.empty() || i + 1 == argc) fail(stem, "--settings requires one filename and cannot be repeated.");
            path = argv[++i];
        } else if (arg == "--seed" || arg == "seed") {
            if (++i == argc) fail(stem, "--seed requires an unsigned integer.");
            std::string value = argv[i];
            if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
                fail(stem, "--seed requires an unsigned integer.");
        } else fail(stem, "unknown option '" + arg + "'.");
    }
    if (!path.empty()) return path;
    if (std::filesystem::exists(stem + ".xml")) return stem + ".xml";
    if (std::filesystem::exists(stem + ".json")) return stem + ".json";
    return stem + ".xml";
}
class Settings {
    std::string path_;
    std::map<std::string, std::string> values_;
public:
    explicit Settings(const std::string &path) : path_(path) {
        if (extension(path) == ".json") {
            auto root = read_json(path);
            for (auto it = root.begin(); it != root.end(); ++it) {
                if (it->is_string()) values_[it.key()] = it->get<std::string>();
                else if (it->is_number()) values_[it.key()] = it->dump();
                else if (it->is_boolean() && it.key() == "Spoil_Active_Pipes")
                    values_[it.key()] = it->get<bool>() ? "yes" : "no";
                else fail(path, "field '" + it.key() + "' must be a scalar string/number (Spoil_Active_Pipes also accepts boolean).");
            }
        } else if (extension(path) == ".xml") {
            auto root = XMLNode::openFileHelper(path.c_str(), "settings");
            for (int i = 0; i < root.nChildNode(); ++i) {
                auto child = root.getChildNode(i);
                std::string key = child.getName();
                if (values_.count(key)) fail(path, "duplicate field '" + key + "'.");
                values_[key] = child.getText() ? child.getText() : "";
            }
        } else fail(path, "settings extension must be .xml or .json.");
    }
    const char *text(const char *field) const {
        auto it = values_.find(field);
        if (it == values_.end() || it->second.empty()) fail(path_, "missing or empty field '" + std::string(field) + "'.");
        return it->second.c_str();
    }
    double real(const char *field) const {
        std::string value = text(field);
        try {
            std::size_t used = 0;
            double result = std::stod(value, &used);
            if (value.find_first_not_of(" \t\r\n", used) != std::string::npos || !std::isfinite(result))
                fail(path_, "field '" + std::string(field) + "' must be a finite number.");
            return result;
        } catch (const std::invalid_argument &) { fail(path_, "field '" + std::string(field) + "' must be numeric."); }
        catch (const std::out_of_range &) { fail(path_, "field '" + std::string(field) + "' is out of range."); }
    }
    int integer(const char *field) const {
        double value = real(field);
        if (value != std::floor(value) || value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
            fail(path_, "field '" + std::string(field) + "' must be an integer in range.");
        return static_cast<int>(value);
    }
};
}
