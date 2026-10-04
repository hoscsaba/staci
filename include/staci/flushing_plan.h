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
  bool timing_allowed = true; // Pipes require abs(v) > threshold; pumps may connect.
};
struct PipeTravel {
  std::string pipe, status;
  double seconds = 0;
  std::vector<TravelArc> route;
};
struct OpeningTime {
  std::string status = "no_qualifying_pipes", critical_pipe;
  double seconds = 0;
  std::vector<PipeTravel> pipes;
};
// Longest allowed advective route; disconnected qualifying pipes are ignored.
OpeningTime opening_time(const std::vector<TravelArc> &arcs,
                         const std::set<std::string> &qualifying,
                         const std::string &hydrant,
                         const std::set<std::string> &storage_nodes);
// Multiple outlets: include routes continuing through an open junction to
// another open outlet, because withdrawal need not consume all incoming flow.
OpeningTime opening_time(const std::vector<TravelArc> &arcs,
                         const std::set<std::string> &qualifying,
                         const std::set<std::string> &hydrants,
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
