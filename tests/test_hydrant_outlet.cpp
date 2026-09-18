#include "HydrantOutlet.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("Hydrant outlet regression failed");
}
int main() {
  HydrantOutlet h("h", "j", 1000, .002, 2);
  h.Set_enabled(true);
  const double head = 20, m = 1000 * .002 * std::sqrt(2 * 9.81 * head / 2);
  h.Set_mp(m);
  std::vector<double> x = {head, 0, 0, 0};
  check(std::abs(h.f(x)) < 1e-12);
  check(std::abs(h.Get_Q() - h.discharge(head)) < 1e-12);
  const auto jac = h.df(x);
  const double delta = 1e-5;
  h.Set_mp(m + delta);
  double plus = h.f(x);
  h.Set_mp(m - delta);
  double minus = h.f(x);
  check(std::abs((plus - minus) / (2 * delta) - jac[2]) < 1e-8);
  h.Set_mp(m);
  x[0] += delta;
  plus = h.f(x);
  x[0] -= 2 * delta;
  minus = h.f(x);
  check(std::abs((plus - minus) / (2 * delta) - jac[0]) < 1e-8);
  HydrantOutlet bigger("b", "j", 1000, .004, 2), lossy("l", "j", 1000, .002, 8);
  check(std::abs(bigger.discharge(head) / h.discharge(head) - 2) < 1e-12);
  check(std::abs(lossy.discharge(head) / h.discharge(head) - .5) < 1e-12);
  h.Set_enabled(false);
  h.Set_mp(3);
  check(h.f(x) == 3 && h.df(x)[2] == 1);
  h.Set_enabled(true);
  h.Set_mp(0);
  x[0] = -1;
  check(h.f(x) == 0 && h.discharge(-1) == 0);
  x[0] = 0;
  check(h.f(x) == 0 && h.df(x)[2] == 1);
  x[0] = 1e-10;
  h.Set_mp(1000 * h.discharge(x[0]));
  check(std::abs(h.f(x)) < 1e-20);
  bool threw = false;
  try {
    HydrantOutlet bad("bad", "j", 1000, 0, 2);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  check(threw);
  std::cout << "Hydrant law, Jacobian, scaling and active-set checks passed\n";
}
