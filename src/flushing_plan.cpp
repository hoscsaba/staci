#include "flushing_plan.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
namespace flushing {
OpeningTime opening_time(const std::vector<TravelArc> &arcs,
                         const std::set<std::string> &qualifying,
                         const std::string &hydrant,
                         const std::set<std::string> &storage_nodes) {
  OpeningTime result;
  if (qualifying.empty())
    return result;
  std::map<std::string, std::vector<const TravelArc *>> outgoing, incoming;
  std::map<std::string, const TravelArc *> by_id;
  for (const auto &a : arcs) {
    if (!std::isfinite(a.seconds) || a.seconds < 0)
      throw std::invalid_argument("Invalid transport time for " + a.id);
    by_id[a.id] = &a;
    // A path terminates at storage or at the operating hydrant.
    if (storage_nodes.count(a.from) || a.from == hydrant)
      continue;
    outgoing[a.from].push_back(&a);
    incoming[a.to].push_back(&a);
  }
  std::set<std::string> reachable{hydrant};
  std::vector<std::string> queue{hydrant};
  for (size_t i = 0; i < queue.size(); ++i)
    for (const auto *a : incoming[queue[i]])
      if (reachable.insert(a->from).second)
        queue.push_back(a->from);
  std::map<std::string, int> color;
  std::map<std::string, double> memo;
  const double infinity = std::numeric_limits<double>::infinity();
  std::function<double(const std::string &)> longest =
      [&](const std::string &node) {
        if (node == hydrant)
          return 0.0;
        if (color[node] == 1)
          return infinity; // reachable recirculation: no finite upper bound
        if (color[node] == 2)
          return memo.at(node);
        color[node] = 1;
        double value = 0;
        for (const auto *a : outgoing[node])
          if (reachable.count(a->to))
            value = std::max(value, a->seconds + longest(a->to));
        color[node] = 2;
        return memo[node] = value;
      };
  bool unavailable = false;
  result.status = "advective_estimate";
  for (const auto &id : qualifying) {
    PipeTravel p{id, "ok", 0};
    const auto found = by_id.find(id);
    if (found == by_id.end() || !reachable.count(found->second->to)) {
      p.status = "no_path_to_hydrant";
    } else {
      p.seconds = found->second->seconds + longest(found->second->to);
      if (!std::isfinite(p.seconds))
        p.status = "flow_cycle";
    }
    if (p.status != "ok")
      unavailable = true;
    else if (result.critical_pipe.empty() || p.seconds > result.seconds) {
      result.seconds = p.seconds;
      result.critical_pipe = id;
    }
    result.pipes.push_back(p);
  }
  if (unavailable) {
    result.status = "undetermined";
    result.seconds = std::numeric_limits<double>::quiet_NaN();
    result.critical_pipe.clear();
  }
  return result;
}
std::vector<PlanStep>
rank_hydrants(const std::vector<PlanCandidate> &candidates,
              const std::map<std::string, double> &volumes) {
  std::set<std::string> covered;
  std::vector<bool> used(candidates.size(), false);
  std::vector<PlanStep> result;
  double cumulative = 0;
  while (result.size() < candidates.size()) {
    size_t best = candidates.size();
    double best_gain = -1;
    for (size_t i = 0; i < candidates.size(); ++i)
      if (!used[i]) {
        double gain = 0;
        for (const auto &id : candidates[i].pipes)
          if (!covered.count(id))
            gain += volumes.at(id);
        if (best == candidates.size() || gain > best_gain ||
            (gain == best_gain &&
             candidates[i].hydrant < candidates[best].hydrant)) {
          best = i;
          best_gain = gain;
        }
      }
    PlanStep step{best, 0, best_gain, cumulative + best_gain, {}};
    for (const auto &id : candidates[best].pipes) {
      step.total_volume += volumes.at(id);
      if (covered.insert(id).second)
        step.added_pipes.push_back(id);
    }
    cumulative = step.cumulative_volume;
    used[best] = true;
    result.push_back(step);
  }
  return result;
}
} // namespace flushing
