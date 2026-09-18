#ifndef STACI_FLUSHING_PLAN_H
#define STACI_FLUSHING_PLAN_H
#include <map>
#include <set>
#include <string>
#include <vector>
namespace flushing {
struct TravelArc {
  std::string id, from, to;
  double seconds;
};
struct PipeTravel {
  std::string pipe, status;
  double seconds = 0;
};
struct OpeningTime {
  std::string status = "no_qualifying_pipes", critical_pipe;
  double seconds = 0;
  std::vector<PipeTravel> pipes;
};
// Longest advective route, starting at each qualifying pipe's upstream end.
OpeningTime opening_time(const std::vector<TravelArc> &arcs,
                         const std::set<std::string> &qualifying,
                         const std::string &hydrant,
                         const std::set<std::string> &storage_nodes);
struct PlanCandidate {
  std::string hydrant, node;
  std::set<std::string> pipes;
  OpeningTime opening;
  double hydrant_flow_m3s = 0;
};
struct PlanStep {
  size_t candidate;
  double total_volume, additional_volume, cumulative_volume;
  std::vector<std::string> added_pipes;
};
std::vector<PlanStep>
rank_hydrants(const std::vector<PlanCandidate> &candidates,
              const std::map<std::string, double> &volumes);
} // namespace flushing
#endif
