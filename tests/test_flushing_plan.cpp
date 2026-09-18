#include "flushing_plan.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace flushing;
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("Flushing plan regression failed");
}
int main() {
  const std::map<std::string, double> volumes{
      {"p1", 10}, {"p2", 8}, {"p3", 7}, {"p4", 9}};
  std::vector<PlanCandidate> candidates{{"B", "B", {"p1", "p3"}, {}},
                                        {"D", "D", {"p1"}, {}},
                                        {"C", "C", {"p4"}, {}},
                                        {"A", "A", {"p1", "p2"}, {}}};
  const auto ranked = rank_hydrants(candidates, volumes);
  check(ranked.size() == 4);
  check(candidates[ranked[0].candidate].hydrant == "A" &&
        ranked[0].additional_volume == 18);
  check(candidates[ranked[1].candidate].hydrant == "C" &&
        ranked[1].additional_volume == 9);
  check(candidates[ranked[2].candidate].hydrant == "B" &&
        ranked[2].additional_volume == 7);
  check(candidates[ranked[3].candidate].hydrant == "D" &&
        ranked[3].additional_volume == 0);
  check(ranked[3].cumulative_volume == 34);
  check(rank_hydrants({}, volumes).empty());
  auto ties =
      rank_hydrants({{"Z", "Z", {"p1"}, {}}, {"A", "A", {"p1"}, {}}}, volumes);
  check(ties[0].candidate == 1);
  std::vector<TravelArc> arcs{{"first", "R", "A", 5},  {"p", "A", "B", 10},
                              {"short", "B", "H", 20}, {"slow1", "B", "C", 30},
                              {"slow2", "C", "H", 40}, {"away", "B", "X", 2},
                              {"cycle1", "U", "V", 1}, {"cycle2", "V", "U", 1}};
  auto timing = opening_time(arcs, {"first", "p"}, "H", {"R"});
  check(timing.status == "advective_estimate" && timing.seconds == 85 &&
        timing.critical_pipe == "first");
  check(opening_time(arcs, {}, "H", {}).seconds == 0);
  // Downstream connecting pipes are counted even when they do not qualify.
  check(opening_time(arcs, {"p"}, "H", {}).seconds == 80);
  auto blocked = opening_time(arcs, {"p"}, "H", {"B"});
  check(blocked.status == "undetermined" && !std::isfinite(blocked.seconds));
  auto unreachable = opening_time(arcs, {"p", "away"}, "H", {});
  check(unreachable.status == "undetermined" &&
        unreachable.pipes[0].status == "no_path_to_hydrant");
  arcs.push_back({"return", "C", "B", 1});
  auto cyclic = opening_time(arcs, {"p"}, "H", {});
  check(cyclic.status == "undetermined" &&
        cyclic.pipes[0].status == "flow_cycle");
  // Pump transit can be approximated as zero; the surrounding pipe times
  // remain.
  check(opening_time({{"p", "R", "B", 5}, {"pump", "B", "H", 0}}, {"p"}, "H",
                     {"R"})
            .seconds == 5);
  std::cout << "Greedy ranking and transport graph checks passed\n";
}
