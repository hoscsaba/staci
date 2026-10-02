#ifndef STACI_DIAGNOSTICS_H
#define STACI_DIAGNOSTICS_H
#include <functional>
#include <stdexcept>
#include <string>

namespace diagnostics {
enum ExitCode { success = 0, calculation_error = 1, input_error = 2, partial_failure = 3 };
class Error : public std::runtime_error {
public:
    Error(std::string code, std::string message, int exit_code = input_error)
        : std::runtime_error(std::move(message)), code(std::move(code)), exit_code(exit_code) {}
    std::string code;
    int exit_code;
};
// A single CLI boundary for all four applications. Common options are removed
// before the original application receives argc/argv.
int run(const char *program, int argc, char **argv,
        const std::function<int(int, char **)> &application);
// Common CLI overrides apply to every solver constructed during this run.
void apply_solver_overrides(double &head_m, double &mass_kgs, int &iterations);
void warning(const std::string &code, const std::string &message);
void error(const std::string &code, const std::string &message);
[[noreturn]] void fail_legacy(const char *source, int line);
// Numerical failures in discarded optimizer candidates are warnings.
class CandidateScope {
public:
    CandidateScope();
    ~CandidateScope();
    CandidateScope(const CandidateScope &) = delete;
    CandidateScope &operator=(const CandidateScope &) = delete;
};
}
#endif
