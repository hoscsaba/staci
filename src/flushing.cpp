#include "flushing.h"
#include "HydrantOutlet.h"
#include "Staci.h"
#include "StaciException.h"
#include "epanet_document.h"
#include "flushing_plan.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <set>
#include <sstream>
#include <stdexcept>

namespace flushing {
namespace {
namespace fs = std::filesystem;
using json = nlohmann::json;
constexpr double head_tolerance = 1e-6, mass_tolerance = 1e-6;
struct Options {
  fs::path inp, hydrants, output, config;
  std::vector<std::string> hydrant_nodes;
  bool write_network_files = false;
  double area = 0, k = 0, velocity = 0, min_pressure = 0;
};
struct Hydrant {
  std::string asset, node;
};
using Records = std::map<std::string, std::vector<std::vector<std::string>>>;
struct Input {
  Records records;
  std::set<std::string> junctions, nodes, pipes, pumps, links;
  std::map<std::string, bool> pump_status;
  std::vector<std::string> warnings;
  double density = 1000;
};
void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}
std::string upper(std::string s) {
  for (char &c : s)
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}
double number(const std::string &s) {
  size_t end = 0;
  const double value = std::stod(s, &end);
  require(end == s.size() && std::isfinite(value),
          "Invalid finite number: " + s);
  return value;
}
std::string csv(const std::string &s) {
  std::string out = "\"";
  for (char c : s) {
    if (c == '"')
      out += '"';
    out += c;
  }
  return out + '"';
}
std::string one_decimal(double value) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(1) << value;
  return out.str();
}
std::string two_decimals(double value) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(2) << value;
  return out.str();
}
std::ofstream output(const fs::path &p) {
  std::ofstream out(p);
  require(bool(out), "Cannot write " + p.string());
  out.exceptions(std::ios::failbit | std::ios::badbit);
  out << std::setprecision(12);
  return out;
}
void write_json(const fs::path &p, const json &value) {
  auto out = output(p);
  out << value.dump(2) << '\n';
}
std::string fingerprint(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  require(bool(in), "Cannot read " + p.string());
  uint64_t hash = 14695981039346656037ULL;
  char c;
  while (in.get(c)) {
    hash ^= static_cast<unsigned char>(c);
    hash *= 1099511628211ULL;
  }
  std::ostringstream out;
  out << std::hex << std::setw(16) << std::setfill('0') << hash;
  return out.str();
}
Options parse(int argc, char **argv) {
  std::map<std::string, std::string> args;
  const std::set<std::string> allowed = {"--inp",
                                         "--hydrants",
                                         "--output-dir",
                                         "--hydrant-area-m2",
                                         "--loss-coefficient",
                                         "--velocity-threshold-mps",
                                         "--min-pressure-head-m",
                                         "--snapshot",
                                         "--config"};
  for (int i = 1; i < argc; ++i) {
    std::string key = argv[i];
    require(allowed.count(key) && i + 1 < argc,
            "Unknown option or missing value: " + key);
    require(args.emplace(key, argv[++i]).second, "Repeated option: " + key);
  }
  Options o;
  if (args.count("--config")) {
    const auto config_path = fs::absolute(args.at("--config"));
    o.config = config_path;
    std::ifstream file(config_path);
    require(bool(file), "Cannot read config: " + config_path.string());
    json config;
    file >> config;
    require(config.is_object(), "Config must be a JSON object");
    const std::map<std::string, std::string> keys = {
        {"hydrant_area_m2", "--hydrant-area-m2"},
        {"total_loss_coefficient", "--loss-coefficient"},
        {"velocity_threshold_mps", "--velocity-threshold-mps"},
        {"output_dir", "--output-dir"},
        {"min_pressure_head_m", "--min-pressure-head-m"}};
    for (auto it = config.begin(); it != config.end(); ++it) {
      if (it.key() == "hydrant_node_ids") {
        require(it.value().is_array() && !it.value().empty(),
                "hydrant_node_ids must be a nonempty array");
        for (const auto &node : it.value()) {
          require(node.is_string() && !node.get<std::string>().empty(),
                  "Hydrant node IDs must be nonempty strings");
          o.hydrant_nodes.push_back(node.get<std::string>());
        }
        continue;
      }
      if (it.key() == "write_network_files") {
        require(it.value().is_boolean(), "write_network_files must be boolean");
        o.write_network_files = it.value().get<bool>();
        continue;
      }
      require(keys.count(it.key()), "Unknown config key: " + it.key());
      const auto &option = keys.at(it.key());
      require(!args.count(option),
              "Config setting also supplied on command line: " + option);
      if (it.key() == "output_dir") {
        require(it.value().is_string(), "output_dir must be a string");
        const auto path = it.value().get<std::string>();
        require(!path.empty(), "output_dir must not be empty");
        args[option] = (config_path.parent_path() / fs::path(path))
                           .lexically_normal()
                           .string();
      } else {
        require(it.value().is_number(),
                "Config setting must be numeric: " + it.key());
        args[option] = it.value().dump();
      }
    }
    for (const auto &entry : keys) {
      if (entry.first != "min_pressure_head_m")
        require(config.contains(entry.first),
                "Required config key: " + entry.first);
    }
  }
  for (auto key : {"--inp", "--output-dir", "--hydrant-area-m2",
                   "--loss-coefficient", "--velocity-threshold-mps"})
    require(args.count(key), std::string("Required option: ") + key);
  require(!args.count("--snapshot") || args["--snapshot"] == "initial",
          "Only --snapshot initial is supported");
  require(o.hydrant_nodes.empty() != (args.count("--hydrants") == 0),
          "Provide hydrant_node_ids in config OR --hydrants, but not both");
  o.inp = fs::absolute(args["--inp"]);
  if (args.count("--hydrants"))
    o.hydrants = fs::absolute(args["--hydrants"]);
  o.output = fs::absolute(args["--output-dir"]);
  o.area = number(args["--hydrant-area-m2"]);
  o.k = number(args["--loss-coefficient"]);
  o.velocity = number(args["--velocity-threshold-mps"]);
  if (args.count("--min-pressure-head-m"))
    o.min_pressure = number(args["--min-pressure-head-m"]);
  require(o.area > 0 && o.k > 0 && o.velocity > 0 && o.min_pressure >= 0,
          "Area, total K and velocity threshold must be positive; minimum "
          "pressure must be nonnegative");
  require(upper(o.inp.extension().string()) == ".INP",
          "--inp must be an EPANET .inp file");
  require(fs::is_regular_file(o.inp) &&
              (o.hydrants.empty() || fs::is_regular_file(o.hydrants)),
          "Input file not found");
  require(!fs::exists(o.output) ||
              (fs::is_directory(o.output) && fs::is_empty(o.output)),
          "Output directory must be new or empty (existing results are never "
          "overwritten)");
  return o;
}
Input inspect(const fs::path &path) {
  Input in;
  auto document = EpanetDocument::read(path.string());
  for (const auto &line : document.lines()) {
    if (line.section_header || line.content.empty())
      continue;
    std::istringstream words(line.content);
    std::vector<std::string> row;
    std::string word;
    while (words >> word)
      row.push_back(word);
    if (!row.empty())
      in.records[upper(line.section)].push_back(row);
  }
  auto &r = in.records;
  for (auto section : {"VALVES", "EMITTERS"})
    require(r[section].empty(),
            std::string("Unsupported hydraulic section [") + section + "]");
  std::set<std::string> ids;
  for (auto section : {"JUNCTIONS", "RESERVOIRS", "TANKS", "PIPES", "PUMPS"}) {
    size_t minimum = std::string(section) == "TANKS"   ? 6
                     : std::string(section) == "PIPES" ? 6
                     : std::string(section) == "PUMPS" ? 5
                                                       : 2;
    for (const auto &row : r[section]) {
      require(row.size() >= minimum,
              std::string("Incomplete [") + section + "] row: " + row[0]);
      require(ids.insert(row[0]).second, "Duplicate node/link ID: " + row[0]);
      if (std::string(section) == "PIPES" || std::string(section) == "PUMPS")
        in.links.insert(row[0]);
      else
        in.nodes.insert(row[0]);
    }
  }
  for (const auto &row : r["JUNCTIONS"]) {
    in.junctions.insert(row[0]);
    number(row[1]);
    if (row.size() > 2)
      number(row[2]);
  }
  for (const auto &row : r["PIPES"]) {
    in.pipes.insert(row[0]);
    require(in.nodes.count(row[1]) && in.nodes.count(row[2]),
            "Unknown pipe endpoint: " + row[0]);
    require(number(row[3]) > 0 && number(row[4]) > 0 && number(row[5]) > 0,
            "Invalid pipe geometry/roughness: " + row[0]);
    require(row.size() <= 8 && (row.size() < 7 || number(row[6]) == 0) &&
                (row.size() < 8 || upper(row[7]) == "OPEN"),
            "Nonzero minor loss or non-open pipe is unsupported: " + row[0]);
  }
  for (const auto &row : r["PUMPS"]) {
    in.pumps.insert(row[0]);
    in.pump_status[row[0]] = true;
    require(in.nodes.count(row[1]) && in.nodes.count(row[2]),
            "Unknown pump endpoint: " + row[0]);
    require(row.size() == 5 &&
                (upper(row[3]) == "POWER" || upper(row[3]) == "HEAD"),
            "Unsupported pump settings: " + row[0]);
    if (upper(row[3]) == "POWER")
      require(number(row[4]) > 0, "Pump power must be positive");
    else {
      size_t points = 0;
      for (const auto &curve : r["CURVES"])
        if (curve[0] == row[4])
          ++points;
      require(points >= 3,
              "Pump head curve needs at least three points: " + row[0]);
    }
  }
  for (const auto &row : r["STATUS"]) {
    require(row.size() == 2 && in.pumps.count(row[0]) &&
                (upper(row[1]) == "OPEN" || upper(row[1]) == "CLOSED"),
            "Only OPEN/CLOSED pump [STATUS] records supported");
    in.pump_status[row[0]] = upper(row[1]) == "OPEN";
  }
  for (const auto &row : r["OPTIONS"]) {
    const auto key = upper(row[0]);
    if (key == "HEADLOSS")
      require(row.size() == 2 &&
                  (upper(row[1]) == "H-W" || upper(row[1]) == "D-W"),
              "Only H-W and D-W headloss supported");
    if (key == "UNITS") {
      const std::set<std::string> units = {"LPS", "LPM", "MLD", "CMH",  "CMD",
                                           "CFS", "GPM", "MGD", "IMGD", "AFD"};
      require(row.size() == 2 && units.count(upper(row[1])),
              "Unsupported INP flow units");
    }
    if (key == "SPECIFIC") {
      require(row.size() == 3 && upper(row[1]) == "GRAVITY",
              "Invalid specific gravity option");
      in.density = 1000 * number(row[2]);
    }
    if (key == "VISCOSITY")
      require(row.size() == 2 && number(row[1]) == 1,
              "Non-default viscosity is not represented by STACI importer");
    if (key == "DEMAND" && row.size() > 1 && upper(row[1]) == "MODEL")
      require(row.size() == 3 && upper(row[2]) == "DDA",
              "Pressure-dependent customer demands are not supported");
  }
  require(in.density > 0 && !in.junctions.empty() && !in.pipes.empty(),
          "Empty network or invalid density");
  require(!r["RESERVOIRS"].empty() || !r["TANKS"].empty(),
          "At least one fixed-head supply is required");
  for (auto section : {"PATTERNS", "CURVES"})
    for (const auto &row : r[section]) {
      require(row.size() >= (std::string(section) == "CURVES" ? 3u : 2u),
              "Incomplete pattern/curve");
      for (size_t i = 1; i < row.size(); ++i)
        number(row[i]);
    }
  for (const auto &row : r["DEMANDS"]) {
    require(row.size() >= 2 && in.junctions.count(row[0]),
            "Invalid demand node");
    number(row[1]);
  }
  for (const auto &row : r["RESERVOIRS"])
    number(row[1]);
  for (const auto &row : r["TANKS"])
    for (size_t i = 1; i < std::min(size_t(7), row.size()); ++i)
      number(row[i]);
  std::set<std::string> patterns;
  for (const auto &row : r["PATTERNS"])
    patterns.insert(row[0]);
  auto require_pattern = [&](const std::string &id) {
    require(patterns.count(id), "Unknown pattern: " + id);
  };
  for (const auto &row : r["JUNCTIONS"])
    if (row.size() > 3)
      require_pattern(row[3]);
  for (const auto &row : r["DEMANDS"])
    if (row.size() > 2)
      require_pattern(row[2]);
  for (const auto &row : r["RESERVOIRS"])
    if (row.size() > 2)
      require_pattern(row[2]);
  for (const auto &row : r["OPTIONS"])
    if (upper(row[0]) == "PATTERN") {
      require(row.size() == 2, "Invalid default pattern option");
      require_pattern(row[1]);
    }
  in.warnings.push_back(
      "Initial snapshot: first demand/reservoir pattern multiplier; tank "
      "initial levels fixed; no elapsed-time simulation.");
  if (!r["RULES"].empty() || !r["CONTROLS"].empty())
    in.warnings.push_back("RULES/CONTROLS are frozen, not executed; all "
                          "initial pump states remain fixed during flushing.");
  in.warnings.push_back(
      "Customer demands remain pressure-independent; coverage does not prove "
      "sediment removal or flushing duration.");
  return in;
}
std::vector<Hydrant> load_hydrants(const fs::path &p, const Input &in) {
  std::vector<Hydrant> result;
  std::ifstream file(p);
  require(bool(file), "Cannot read hydrants");
  if (upper(p.extension().string()) == ".JSON") {
    json j;
    file >> j;
    require(j.contains("hydrants") && j["hydrants"].is_array(),
            "JSON requires a hydrants array");
    for (const auto &h : j["hydrants"]) {
      require(h.value("status", "") == "matched" && h.contains("node_id") &&
                  h["node_id"].is_string(),
              "Unmatched or invalid JSON hydrant");
      std::string node = h["node_id"].get<std::string>();
      result.push_back({h.value("dxf_handle", node), node});
    }
  } else {
    std::string line;
    while (std::getline(file, line)) {
      auto start = line.find_first_not_of(" \t\r");
      if (start == std::string::npos || line[start] == '#')
        continue;
      auto node =
          line.substr(start, line.find_last_not_of(" \t\r") - start + 1);
      result.push_back({node, node});
    }
  }
  std::set<std::string> assets;
  for (const auto &h : result) {
    require(!h.asset.empty() && assets.insert(h.asset).second,
            "Duplicate/empty hydrant asset ID: " + h.asset);
    require(in.junctions.count(h.node),
            "Hydrant is not an existing junction: " + h.node);
  }
  require(!result.empty(), "Empty hydrant list");
  return result;
}
// Staci owns the imported network. Keep only the added outlets here; they
// outlive the solver and its non-owning compatibility views.
struct Network {
  std::vector<std::unique_ptr<HydrantOutlet>> added_outlets;
  Staci system;
  explicit Network(const std::string &path) : system(path) {}
  Network(const Network &) = delete;
  Network &operator=(const Network &) = delete;
};
struct State {
  std::vector<double> heads, flows;
  explicit State(Staci &s) {
    for (auto *n : s.cspok)
      heads.push_back(n->Get_p());
    for (auto *e : s.agelemek)
      flows.push_back(e->Get_mp());
  }
  void restore(Staci &s) const {
    for (size_t i = 0; i < heads.size(); ++i)
      s.cspok[i]->Set_p(heads[i]);
    for (size_t i = 0; i < flows.size(); ++i)
      s.agelemek[i]->Set_mp(flows[i]);
  }
};
struct Check {
  bool finite = true;
  double min_head = std::numeric_limits<double>::infinity(), balance = 0,
         residual = 0;
  std::string min_node;
};
Check check(Staci &s, const Input &in) {
  Check c;
  for (auto *n : s.cspok) {
    c.finite &= std::isfinite(n->Get_p());
    if (in.junctions.count(n->Get_nev()) && n->Get_p() < c.min_head) {
      c.min_head = n->Get_p();
      c.min_node = n->Get_nev();
    }
    double balance = -n->Get_fogy();
    for (int i : n->ag_be)
      balance += s.agelemek[i]->Get_mp();
    for (int i : n->ag_ki)
      balance -= s.agelemek[i]->Get_mp();
    c.finite &= std::isfinite(balance);
    c.balance = std::max(c.balance, std::abs(balance));
  }
  for (auto *e : s.agelemek) {
    c.finite &= std::isfinite(e->Get_mp()) && std::isfinite(e->Get_v());
    auto *a = s.cspok[e->Get_Cspe_Index()];
    std::vector<double> x = {a->Get_p(), 0, a->Get_h(), 0};
    if (e->Get_Csp_db() == 2) {
      auto *b = s.cspok[e->Get_Cspv_Index()];
      x[1] = b->Get_p();
      x[3] = b->Get_h();
    }
    double r = e->f(x);
    c.finite &= std::isfinite(r);
    c.residual = std::max(c.residual, std::abs(r));
  }
  return c;
}
std::string status(bool converged, const Check &c, double minimum) {
  if (!converged)
    return "not_converged";
  if (!c.finite)
    return "nonfinite";
  if (c.balance > 10 * mass_tolerance || c.residual > 10 * head_tolerance)
    return "residual_failed";
  if (c.min_head < -head_tolerance)
    return "negative_pressure";
  if (c.min_head < minimum)
    return "below_min_pressure";
  return "ok";
}
void verify_supply(Staci &s, const Input &in) {
  std::map<std::string, size_t> index;
  for (size_t i = 0; i < s.cspok.size(); ++i)
    index[s.cspok[i]->Get_nev()] = i;
  std::vector<std::vector<size_t>> graph(s.cspok.size());
  std::vector<bool> visited(s.cspok.size(), false);
  std::vector<size_t> queue;
  for (auto *e : s.agelemek) {
    if (!e->Is_enabled())
      continue;
    if (e->Get_Csp_db() == 1) {
      auto i = index.at(e->Get_Cspe_Nev());
      if (!visited[i]) {
        visited[i] = true;
        queue.push_back(i);
      }
    } else {
      auto a = index.at(e->Get_Cspe_Nev()), b = index.at(e->Get_Cspv_Nev());
      graph[a].push_back(b);
      graph[b].push_back(a);
    }
  }
  for (size_t k = 0; k < queue.size(); ++k)
    for (auto i : graph[queue[k]])
      if (!visited[i]) {
        visited[i] = true;
        queue.push_back(i);
      }
  for (size_t i = 0; i < visited.size(); ++i)
    require(visited[i],
            "Node has no enabled supply path: " + s.cspok[i]->Get_nev());
}
// EPANET emitter coefficient uses flow units per sqrt(metre or psi).
void export_network(const fs::path &path, const Options &o, const Input &in,
                    const Hydrant &h, const std::string &state,
                    const std::map<std::string, Csomopont *> &nodes,
                    const std::set<std::string> &flushed) {
  std::string units = "LPS";
  for (const auto &row : in.records.at("OPTIONS"))
    if (upper(row[0]) == "UNITS")
      units = upper(row[1]);
  const std::map<std::string, double> flow_factor = {
      {"LPS", 1000},
      {"LPM", 60000},
      {"MLD", 86.4},
      {"CMH", 3600},
      {"CMD", 86400},
      {"CFS", 1 / 0.028316846592},
      {"GPM", 1 / 0.0000630901964},
      {"MGD", 1 / 0.0438126363888889},
      {"IMGD", 1 / 0.0526167824074074},
      {"AFD", 1 / 0.0142764101851852}};
  const bool us = units == "CFS" || units == "GPM" || units == "MGD" ||
                  units == "IMGD" || units == "AFD";
  const double length_factor = us ? 1 / 0.3048 : 1;
  const double head_per_pressure =
      us ? 0.3048 / (0.4333 * in.density / 1000) : 1;
  auto out = output(path);
  out << "; STACI flushing snapshot: hydrant " << h.node << ", status " << state
      << "\n; Tags describe STACI results, not recomputed EPANET results.\n";
  const std::set<std::string> replaced = {
      "JUNCTIONS", "RESERVOIRS", "DEMANDS",  "PATTERNS", "CONTROLS", "RULES",
      "TIMES",     "OPTIONS",    "EMITTERS", "TAGS",     "END"};
  out << "\n[JUNCTIONS]\n";
  for (const auto &id : in.junctions) {
    auto *n = nodes.at(id);
    out << id << ' ' << n->Get_h() * length_factor << ' '
        << n->Get_fogy() / in.density * flow_factor.at(units) << '\n';
  }
  out << "\n[RESERVOIRS]\n";
  for (const auto &row : in.records.at("RESERVOIRS")) {
    double head = number(row[1]);
    if (row.size() > 2) {
      for (const auto &pattern : in.records.at("PATTERNS"))
        if (pattern[0] == row[2]) {
          head *= number(pattern[1]);
          break;
        }
    }
    out << row[0] << ' ' << head << '\n';
  }
  const auto document = EpanetDocument::read(o.inp.string());
  for (const auto &line : document.lines()) {
    if (replaced.count(upper(line.section)))
      continue;
    out << line.raw << '\n';
  }
  out << "\n[OPTIONS]\n";
  for (const auto &row : in.records.at("OPTIONS")) {
    const auto key = upper(row[0]);
    if (key == "UNITS" || key == "PATTERN" || key == "DEMAND" ||
        key == "EMITTER" || key == "PRESSURE" || key == "HYDRAULICS" ||
        key == "QUALITY")
      continue;
    for (const auto &word : row)
      out << word << ' ';
    out << '\n';
  }
  out << "UNITS " << units
      << "\nDEMAND MULTIPLIER 1\nDEMAND MODEL DDA\nEMITTER EXPONENT "
         "0.5\nQUALITY "
         "NONE\n"
      << "\n[TIMES]\nDURATION 0\n\n[EMITTERS]\n"
      << h.node << ' '
      << o.area * std::sqrt(2 * 9.81 * head_per_pressure / o.k) *
             flow_factor.at(units)
      << '\n';
  out << "\n[TAGS]\n";
  for (const auto &row : in.records.at("TAGS")) {
    if (row.size() >= 2 && upper(row[0]) == "LINK" && in.pipes.count(row[1]))
      continue;
    for (const auto &word : row)
      out << word << ' ';
    out << '\n';
  }
  for (const auto &id : in.pipes)
    out << "LINK " << id << ' '
        << (state != "ok"       ? "INVALID_SCENARIO"
            : flushed.count(id) ? "FLUSHED"
                                : "BELOW_THRESHOLD")
        << '\n';
  out << "\n[END]\n";
}
int execute(const Options &o) {
  auto in = inspect(o.inp);
  in.records["TAGS"];
  std::vector<Hydrant> hydrants;
  if (!o.hydrants.empty())
    hydrants = load_hydrants(o.hydrants, in);
  else {
    std::set<std::string> seen;
    for (const auto &node : o.hydrant_nodes) {
      require(in.junctions.count(node),
              "Hydrant is not an existing junction: " + node);
      require(seen.insert(node).second, "Duplicate hydrant node ID: " + node);
      hydrants.push_back({node, node});
    }
  }
  if (o.write_network_files)
    for (const auto &h : hydrants)
      require(h.node.find_first_of("/\\") == std::string::npos,
              "Hydrant node ID contains a path separator: " + h.node);
  fs::create_directories(o.output / "work");
  auto work_input = o.output / "work" / "network.inp";
  fs::copy_file(o.inp, work_input);
  if (!o.config.empty())
    fs::copy_file(o.config, o.output / "config.json");
  if (!o.hydrants.empty())
    fs::copy_file(o.hydrants,
                  o.output / ("hydrants" + o.hydrants.extension().string()));
  json manifest = {
      {"schema_version", 1},
      {"program", "staci_flush"},
      {"status", "running"},
      {"input_path", o.inp.string()},
      {"hydrants_path", o.hydrants.string()},
      {"input_fnv1a64", fingerprint(o.inp)},
      {"hydrants_fnv1a64", o.hydrants.empty() ? "" : fingerprint(o.hydrants)},
      {"hydrant_node_ids", o.hydrant_nodes},
      {"write_network_files", o.write_network_files},
      {"snapshot", "initial"},
      {"hydrant_area_m2", o.area},
      {"total_loss_coefficient", o.k},
      {"loss_convention", "h=K*(Q/A)^2/(2*g); K includes outlet kinetic head"},
      {"outlet_height_above_junction_m", 0},
      {"velocity_threshold_mps", o.velocity},
      {"comparison", "abs(v)>threshold"},
      {"min_pressure_head_m", o.min_pressure},
      {"pressure_head_tolerance_m", head_tolerance},
      {"mass_tolerance_kgs", mass_tolerance},
      {"initial_pump_states", in.pump_status},
      {"warnings", in.warnings}};
  write_json(o.output / "run.json", manifest);
  for (const auto &w : in.warnings)
    std::cerr << "Warning: " << w << '\n';
  // If a later exception escapes, do not leave a misleading running manifest.
  struct FailureManifest {
    fs::path path;
    json &value;
    ~FailureManifest() {
      if (value["status"] == "running") {
        value["status"] = "fatal_error";
        try {
          write_json(path, value);
        } catch (...) {
        }
      }
    }
  } failure_manifest{o.output / "run.json", manifest};
  Network network(work_input.string());
  auto &s = network.system;
  s.Set_debug_level(0);
  s.set_do_save_file(false);
  std::map<std::string, Csomopont *> nodes;
  std::map<std::string, Agelem *> edges;
  for (auto *n : s.cspok)
    nodes.emplace(n->Get_nev(), n);
  for (auto *e : s.agelemek)
    edges.emplace(e->Get_nev(), e);
  require(nodes.size() == in.nodes.size(), "Importer changed node count");
  for (const auto &id : in.nodes)
    require(nodes.count(id), "Importer omitted node: " + id);
  for (const auto &id : in.links)
    require(edges.count(id), "Importer omitted link: " + id);
  require(edges.size() == in.links.size() + in.records["RESERVOIRS"].size() +
                              in.records["TANKS"].size(),
          "Unexpected imported link count");
  for (const auto &p : in.pump_status)
    edges.at(p.first)->Set_enabled(p.second);
  verify_supply(s, in);
  std::map<std::string, HydrantOutlet *> outlets;
  std::set<std::string> outlet_nodes;
  for (const auto &h : hydrants)
    outlet_nodes.insert(h.node);
  for (const auto &node : outlet_nodes) {
    std::string id = "__FLUSH_" + node;
    while (nodes.count(id) || edges.count(id))
      id += '_';
    auto outlet =
        std::make_unique<HydrantOutlet>(id, node, in.density, o.area, o.k);
    outlets[node] = outlet.get();
    s.agelemek.push_back(outlet.get());
    network.added_outlets.push_back(std::move(outlet));
  }
  // The solver uses RMS residuals. Scale tolerances so large networks meet
  // the same per-equation accuracy as small fixtures (RMS * sqrt(N) = L2).
  s.set_solver_tolerances(head_tolerance / std::sqrt(double(s.agelemek.size())),
                          mass_tolerance / std::sqrt(double(s.cspok.size())));
  s.build_system();
  s.ini();
  for (const auto &h : outlets)
    h.second->Set_mp(0);
  bool solved = s.solve_system();
  Check initial = check(s, in);
  // Baseline must be physically valid, even if the optional service target is
  // higher.
  auto baseline_status = status(solved, initial, 0);
  if (baseline_status != "ok") {
    manifest["status"] = "baseline_" + baseline_status;
    manifest["baseline_min_pressure_head_m"] = initial.min_head;
    manifest["baseline_max_mass_residual_kgs"] = initial.balance;
    manifest["baseline_max_edge_residual"] = initial.residual;
    write_json(o.output / "run.json", manifest);
    throw std::runtime_error("Baseline failed: " + baseline_status);
  }
  State baseline(s);
  std::map<std::string, double> base_v, max_v;
  std::map<std::string, std::vector<std::string>> coverage;
  std::map<std::string, std::string> best;
  auto baseline_out = output(o.output / "baseline_pipes.csv");
  baseline_out << "pipe_id,node_from,node_to,length_m,diameter_m,flow_m3s,"
                  "velocity_mps,above_threshold\n";
  for (const auto &id : in.pipes) {
    auto *e = edges.at(id);
    base_v[id] = e->Get_v();
    max_v[id] = -1;
    baseline_out << csv(id) << ',' << csv(e->Get_Cspe_Nev()) << ','
                 << csv(e->Get_Cspv_Nev()) << ',' << e->Get_dprop("length")
                 << ',' << e->Get_dprop("diameter") << ',' << e->Get_Q() << ','
                 << e->Get_v() << ',' << (std::abs(e->Get_v()) > o.velocity)
                 << '\n';
  }
  auto summary = output(o.output / "scenarios.csv");
  summary << "hydrant_id,node_id,status,converged,flow_m3s,flow_lps,hydrant_"
             "pressure_head_m,min_pressure_head_m,min_pressure_node,max_mass_"
             "residual_kgs,max_edge_residual,qualifying_pipe_count,qualifying_"
             "length_m,qualifying_volume_m3\n";
  auto details = output(o.output / "scenario_pipes.csv");
  auto above = output(o.output / "pipes_above_threshold.csv");
  const std::string header =
      "hydrant_id,node_id,pipe_id,flow_m3s,velocity_mps,absolute_velocity_mps,"
      "baseline_velocity_mps,above_threshold,newly_above_threshold\n";
  details << header;
  above << header;
  auto text_summary = output(o.output / "summary.txt");
  std::ofstream network_index;
  if (o.write_network_files) {
    fs::create_directories(o.output / "networks");
    network_index = output(o.output / "networks.csv");
    network_index << "hydrant_id,node_id,status,network_file\n";
  }
  std::vector<PlanCandidate> candidates;
  std::map<std::string, double> pipe_volumes;
  for (const auto &id : in.pipes) {
    auto *e = edges.at(id);
    const double d = e->Get_dprop("diameter");
    pipe_volumes[id] = std::acos(-1.0) / 4 * d * d * e->Get_dprop("length");
  }
  std::set<std::string> storage_nodes;
  for (const auto &id : in.nodes)
    if (!in.junctions.count(id))
      storage_nodes.insert(id);
  auto travel_out = output(o.output / "pipe_travel_times.csv");
  travel_out << "hydrant_id,node_id,pipe_id,status,travel_time_s\n";
  int failures = 0;
  size_t count = 0;
  for (const auto &h : hydrants) {
    baseline.restore(s);
    for (const auto &outlet : outlets) {
      outlet.second->Set_enabled(false);
      outlet.second->Set_mp(0);
    }
    auto *outlet = outlets.at(h.node);
    outlet->Set_enabled(true);
    outlet->Set_mp(in.density * outlet->discharge(nodes.at(h.node)->Get_p()));
    solved = s.solve_system();
    auto c = check(s, in);
    std::string state = status(solved, c, o.min_pressure);
    if (state == "ok" && outlet->Get_Q() < 0)
      state = "reverse_hydrant_flow";
    size_t qualifying = 0;
    double length = 0, volume = 0;
    std::set<std::string> flushed;
    if (state == "ok")
      for (const auto &id : in.pipes) {
        auto *e = edges.at(id);
        double v = e->Get_v();
        bool hit = std::abs(v) > o.velocity;
        std::ostringstream row;
        row << std::setprecision(12) << csv(h.asset) << ',' << csv(h.node)
            << ',' << csv(id) << ',' << e->Get_Q() << ',' << v << ','
            << std::abs(v) << ',' << base_v.at(id) << ',' << hit << ','
            << (hit && std::abs(base_v.at(id)) <= o.velocity) << '\n';
        details << row.str();
        if (hit) {
          above << row.str();
          ++qualifying;
          length += e->Get_dprop("length");
          const double diameter = e->Get_dprop("diameter");
          volume += std::acos(-1.0) / 4 * diameter * diameter *
                    e->Get_dprop("length");
          flushed.insert(id);
          coverage[id].push_back(h.asset);
        }
        if (std::abs(v) > max_v[id] ||
            (std::abs(v) == max_v[id] && h.asset < best[id])) {
          max_v[id] = std::abs(v);
          best[id] = h.asset;
        }
      }
    else
      ++failures;
    if (state == "ok") {
      std::vector<TravelArc> arcs;
      for (const auto &id : in.links) {
        auto *e = edges.at(id);
        if (!e->Is_enabled() || std::abs(e->Get_Q()) <= 1e-12)
          continue;
        std::string from = e->Get_Cspe_Nev(), to = e->Get_Cspv_Nev();
        if (e->Get_Q() < 0)
          std::swap(from, to);
        const double seconds =
            in.pipes.count(id) ? e->Get_dprop("length") / std::abs(e->Get_v())
                               : 0;
        arcs.push_back({id, from, to, seconds});
      }
      auto timing = opening_time(arcs, flushed, h.node, storage_nodes);
      for (const auto &p : timing.pipes) {
        travel_out << csv(h.asset) << ',' << csv(h.node) << ',' << csv(p.pipe)
                   << ',' << p.status << ',';
        if (p.status == "ok")
          travel_out << p.seconds;
        travel_out << '\n';
      }
      candidates.push_back(
          {h.asset, h.node, flushed, std::move(timing), outlet->Get_Q()});
    }
    summary << csv(h.asset) << ',' << csv(h.node) << ',' << state << ','
            << solved << ',';
    if (solved && c.finite)
      summary << outlet->Get_Q() << ',' << 1000 * outlet->Get_Q() << ','
              << nodes.at(h.node)->Get_p() << ',' << c.min_head << ','
              << csv(c.min_node) << ',' << c.balance << ',' << c.residual;
    else
      summary << ",,,,,,";
    summary << ',' << qualifying << ',' << length << ',';
    if (state == "ok")
      summary << two_decimals(volume);
    summary << '\n';
    ++count;
    if (o.write_network_files) {
      const std::string name =
          "networks/" + o.inp.stem().string() + "_hydrant_" + h.node + ".inp";
      export_network(o.output / name, o, in, h, state, nodes, flushed);
      network_index << csv(h.asset) << ',' << csv(h.node) << ',' << state << ','
                    << csv(name) << '\n';
    }
    const std::string volume_text = state == "ok"
                                        ? two_decimals(volume) + " m3"
                                        : "unavailable (invalid scenario)";
    text_summary << "Hydrant " << h.asset << " (" << h.node << "): " << state
                 << ", " << qualifying
                 << " pipes; qualifying pipe volume: " << volume_text << '\n';
    std::cout << "Hydrant " << count << '/' << hydrants.size() << ' ' << h.asset
              << " (" << h.node << "): " << state << ", " << qualifying
              << " pipes; qualifying pipe volume: " << volume_text << "\n"
              << std::flush;
  }
  auto aggregate = output(o.output / "pipe_coverage.csv");
  aggregate << "pipe_id,baseline_above_threshold,hydrant_count,hydrant_ids_"
               "json,max_valid_velocity_mps,best_hydrant_id\n";
  for (const auto &id : in.pipes) {
    auto &list = coverage[id];
    std::sort(list.begin(), list.end());
    aggregate << csv(id) << ',' << (std::abs(base_v.at(id)) > o.velocity) << ','
              << list.size() << ',' << csv(json(list).dump()) << ',';
    if (max_v[id] >= 0)
      aggregate << max_v[id];
    aggregate << ',' << csv(best[id]) << '\n';
  }
  double unique_volume = 0;
  for (const auto &id : in.pipes)
    if (!coverage[id].empty()) {
      auto *e = edges.at(id);
      const double d = e->Get_dprop("diameter");
      unique_volume += std::acos(-1.0) / 4 * d * d * e->Get_dprop("length");
    }
  text_summary << "Unique pipe volume covered by valid scenarios: "
               << two_decimals(unique_volume) << " m3\n";
  std::cout << "Unique pipe volume covered by valid scenarios: "
            << two_decimals(unique_volume) << " m3\n";
  double total_network_volume = 0;
  for (const auto &entry : pipe_volumes)
    total_network_volume += entry.second;
  auto plan = output(o.output / "flushing_plan.csv");
  auto plan_text = output(o.output / "flushing_plan.txt");
  plan << "rank,hydrant_id,node_id,qualifying_volume_m3,additional_volume_m3,"
          "cumulative_volume_m3,cumulative_volume_percent,additional_pipe_ids_"
          "json,"
          "opening_time_min,timing_status,hydrant_flow_m3s,"
          "volume_over_flow_time_min,volume_over_flow_status,critical_pipe_id,"
          "redundant\n";
  plan_text
      << "Greedy volume ranking; each step credits only previously uncovered "
         "qualifying pipes.\n"
      << "Opening time: longest advective pipe travel time to the active "
         "hydrant.\n"
      << "All qualifying pipes are timed, including previously covered pipes.\n"
      << "Volume/flow time: total qualifying pipe volume for this hydrant "
         "divided by its discharge.\n"
      << "These are transport estimates, not sediment-cleaning completion "
         "guarantees.\n"
      << "Branches may carry contamination to consumers; storage terminates "
         "routes.\n\n";
  plan_text << "Overall network pipe volume: "
            << two_decimals(total_network_volume) << " m3\n\n";
  const auto ranked = rank_hydrants(candidates, pipe_volumes);
  size_t undetermined = 0, rank = 0;
  for (const auto &step : ranked) {
    const auto &candidate = candidates[step.candidate];
    const auto &timing = candidate.opening;
    const double percent =
        total_network_volume > 0
            ? 100 * step.cumulative_volume / total_network_volume
            : 0;
    const double volume_flow_seconds =
        step.total_volume == 0
            ? 0
            : (candidate.hydrant_flow_m3s > 0
                   ? step.total_volume / candidate.hydrant_flow_m3s
                   : std::numeric_limits<double>::quiet_NaN());
    const bool volume_flow_valid = std::isfinite(volume_flow_seconds);
    const std::string volume_flow_status =
        step.total_volume == 0
            ? "no_qualifying_pipes"
            : (volume_flow_valid ? "volume_over_flow_estimate"
                                 : "undetermined");
    const bool timed = timing.status != "undetermined";
    if (!timed)
      ++undetermined;
    plan << ++rank << ',' << csv(candidate.hydrant) << ','
         << csv(candidate.node) << ',' << two_decimals(step.total_volume) << ','
         << two_decimals(step.additional_volume) << ','
         << two_decimals(step.cumulative_volume) << ',' << two_decimals(percent)
         << ',' << csv(json(step.added_pipes).dump()) << ',';
    if (timed)
      plan << one_decimal(timing.seconds / 60);
    plan << ',' << timing.status << ',' << candidate.hydrant_flow_m3s << ',';
    if (volume_flow_valid)
      plan << one_decimal(volume_flow_seconds / 60);
    plan << ',' << volume_flow_status << ',' << csv(timing.critical_pipe) << ','
         << (step.additional_volume == 0) << '\n';
    std::ostringstream line;
    line << rank << ". " << candidate.hydrant;
    if (candidate.hydrant != candidate.node)
      line << " (node " << candidate.node << ")";
    line << ": additional " << two_decimals(step.additional_volume)
         << " m3; cumulative " << two_decimals(step.cumulative_volume)
         << " m3 (" << two_decimals(percent) << "%); opening ";
    if (timed)
      line << one_decimal(timing.seconds / 60) << " min";
    else
      line << "UNDETERMINED: see pipe_travel_times.csv";
    line << "; volume/flow time ";
    if (volume_flow_valid)
      line << one_decimal(volume_flow_seconds / 60) << " min";
    else
      line << "UNDETERMINED: nonpositive hydrant discharge";
    if (step.additional_volume == 0)
      line << "; redundant coverage";
    plan_text << line.str() << '\n';
    std::cout << line.str() << '\n';
  }
  plan_text << "\nHydraulically invalid scenarios excluded: " << failures
            << '\n';
  manifest["plan"] = {
      {"method", "greedy_additional_pipe_volume"},
      {"tie_break", "hydrant_id_lexicographic"},
      {"ranked_hydrants", candidates.size()},
      {"undetermined_opening_times", undetermined},
      {"timing", "longest_directed_advective_path_all_qualifying_pipes"},
      {"minimum_transport_flow_m3s", 1e-12},
      {"additional_timing",
       "total_qualifying_pipe_volume_divided_by_hydrant_discharge"}};
  manifest["unique_qualifying_volume_m3"] =
      std::round(unique_volume * 100) / 100;
  manifest["total_network_pipe_volume_m3"] =
      std::round(total_network_volume * 100) / 100;
  manifest["cumulative_volume_percent"] =
      total_network_volume > 0
          ? std::round(10000 * unique_volume / total_network_volume) / 100
          : 0;
  manifest["status"] = failures ? "partial_failure" : "complete";
  manifest["scenario_count"] = hydrants.size();
  manifest["failed_scenarios"] = failures;
  manifest["pipe_count"] = in.pipes.size();
  manifest["baseline_min_pressure_head_m"] = initial.min_head;
  write_json(o.output / "run.json", manifest);
  std::cout << "Flushing complete: " << hydrants.size() << " scenarios, "
            << failures << " failed. Results: " << o.output << '\n';
  return failures ? 2 : 0;
}
} // namespace
int run(int argc, char **argv) {
  if (argc == 2 && std::string(argv[1]) == "--help") {
    std::cout << "Usage: staci_flush --inp network.inp --config "
                 "flushing_config.json\n"
                 "Config: hydrant_node_ids array, optional write_network_files "
                 "boolean.\n"
                 "Alternatively use --hydrants nodes.txt|hydrants.json.\n"
                 "Legacy usage without --config:\n"
                 "  --hydrant-area-m2 A --loss-coefficient K "
                 "--velocity-threshold-mps V\n"
                 "  --output-dir NEW_DIRECTORY [--min-pressure-head-m H] "
                 "[--snapshot initial]\n"
                 "K is TOTAL resistance: h=K*(Q/A)^2/(2g), including outlet "
                 "kinetic head.\n"
                 "Each scenario opens one hydrant; network valve/pump states "
                 "remain fixed.\n"
                 "Exit codes: 0 complete; 1 input/baseline error; 2 partial "
                 "scenario failure.\n";
    return 0;
  }
  try {
    return execute(parse(argc, argv));
  } catch (StaciException &e) {
    std::cerr << "staci_flush: " << e.getDescription() << '\n';
  } catch (const std::exception &e) {
    std::cerr << "staci_flush: " << e.what() << '\n';
  }
  return 1;
}
} // namespace flushing
