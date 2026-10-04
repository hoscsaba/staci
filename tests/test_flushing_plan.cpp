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
  // Allowed downstream connecting pipes are included.
  check(opening_time(arcs, {"p"}, "H", {}).seconds == 80);
  auto blocked = opening_time(arcs, {"p"}, "H", {"B"});
  check(blocked.status == "no_connected_qualifying_pipes" && blocked.seconds == 0);
  auto unreachable = opening_time(arcs, {"p", "away"}, "H", {});
  check(unreachable.status == "advective_estimate" && unreachable.seconds == 80 &&
        unreachable.pipes[0].status == "ignored_no_qualifying_path");
  check(timing.pipes.front().route.size() == 4);
  check(timing.pipes.front().route[2].id == "slow1");
  // A below-threshold connector cuts off the remote qualifying pipe.
  auto slow = opening_time({{"fast", "R", "A", 100},
                            {"slow", "A", "H", 100000, false}}, {"fast"}, "H", {"R"});
  check(slow.seconds == 0 && slow.status == "no_connected_qualifying_pipes" &&
        slow.pipes[0].status == "ignored_no_qualifying_path" && slow.pipes[0].route.empty());
  auto filtered = arcs;
  filtered[3].timing_allowed = false;
  auto alternate = opening_time(filtered, {"first", "p", "away"}, "H", {"R"});
  check(alternate.seconds == 35 && alternate.status == "advective_estimate");
  check(alternate.pipes[1].route.back().id == "short");
  // Two independently fed outlets; keep the longest reachable route.
  const std::set<std::string> both{"H1", "H2"};
  auto multi = opening_time({{"a", "R", "H1", 10},
                             {"b", "R", "H2", 20}}, {"a", "b"}, both, {"R"});
  check(multi.seconds == 20 && multi.critical_pipe == "b");
  // Withdrawal at H1 need not consume all incoming water: some reaches H2.
  multi = opening_time({{"a", "R", "H1", 10}, {"b", "H1", "H2", 20}},
                       {"a", "b"}, both, {"R"});
  check(multi.seconds == 30 && multi.pipes[0].route.back().to == "H2");
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
