// Exercise the real optimizer objective without adding a public test-only CLI.
#define main calibration_application_main
#include "../src/staci_calibrate.cpp"
#undef main

static void check(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

static void check_pool_balance() {
    for (std::size_t i = 0; i < wds.size(); ++i) {
        for (auto *edge : wds[i]->agelemek) {
            if (edge->GetType() != "Vegakna") continue;
            auto *node = wds[i]->cspok.at(edge->Get_Cspe_Index());
            check(std::abs(node->Get_p() + node->Get_h() -
                  edge->Get_dprop("bottom_level") - edge->Get_dprop("water_level")) < 0.01,
                  "Solved hydraulic head disagrees with the tank boundary");
            if (i == 0) continue;
            Agelem *previous = nullptr;
            // Locate the preceding tank by ID, independently of node ordering.
            for (auto *candidate : wds[i - 1]->agelemek)
                if (candidate->Get_nev() == edge->Get_nev()) previous = candidate;
            check(previous != nullptr, "Preceding tank is missing");
            const double expected = previous->Get_dprop("water_level") +
                previous->Get_Q() * dt * 3600.0 / edge->Get_Aref();
            check(std::abs(expected - edge->Get_dprop("water_level")) < 1e-9,
                  "Tank storage balance does not use the current candidate's preceding period");
        }
    }
}

int main(int argc, char **argv) {
    try {
        check(run_application(argc, argv) == 0, "Calibration failed");
        pagmo::vector_double a, b;
        for (std::size_t i = 0; i < pipe_name.size(); ++i)
            if (pipe_is_active[i]) {
                a.push_back(pipe_origD[i] * 1.05);
                b.push_back(pipe_origD[i] * 0.95);
            }
        best_obj = -1.0; // Avoid best-result output during the probe.
        const double first = Objective(a);
        check_pool_balance();
        const double repeated = Objective(a);
        check_pool_balance();
        Objective(b);
        const double after_other = Objective(a);
        check_pool_balance();
        check(std::isfinite(first) && std::abs(first - repeated) < 1e-7 &&
              std::abs(first - after_other) < 1e-7,
              "Objective depends on previous candidate evaluations");
        std::cout << "PASS: repeated and interleaved candidates preserve objective and tank balance\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
