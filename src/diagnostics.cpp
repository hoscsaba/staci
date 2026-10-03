#include "diagnostics.h"
#include <cmath>
#include <limits>
#include <optional>
#include "epanet_reader.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <ctime>
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <regex>
#include <sstream>
#include <streambuf>
#include <vector>
#ifdef _WIN32
// Keep Windows min/max macros from colliding with the C++ standard library.
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <process.h>
#else
#include <fcntl.h>
#include <unistd.h>
#include <sys/file.h>
#endif
namespace diagnostics {
namespace {
using json = nlohmann::json;
thread_local unsigned candidate_depth = 0;
std::string timestamp() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t seconds = std::chrono::system_clock::to_time_t(now);
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &seconds);
#else
    gmtime_r(&seconds, &utc);
#endif
    std::ostringstream out;
    out << std::put_time(&utc, "%Y-%m-%dT%H:%M:%S") << '.' << std::setfill('0')
        << std::setw(3) << std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count() % 1000 << 'Z';
    return out.str();
}
class Session;
Session *active = nullptr;
class Capture : public std::streambuf {
public:
    Capture(Session &session, std::streambuf *original, bool echo) : session_(session), original_(original), echo_(echo) {}
    int sync() override;
    void finish();
protected:
    int_type overflow(int_type c) override;
    std::streamsize xsputn(const char *s, std::streamsize n) override;
private:
    Session &session_;
    std::streambuf *original_;
    std::string line_;
    bool echo_;
    std::mutex capture_mutex_;
};
class Session {
public:
    std::string program, run_id, recent, last_error;
    unsigned errors = 0, warnings = 0;
#ifdef _WIN32
    HANDLE file = INVALID_HANDLE_VALUE;
#else
    int file = -1;
#endif
    std::streambuf *out = std::cout.rdbuf(), *err = std::cerr.rdbuf();
    Capture capture_out{*this, out, true}, capture_err{*this, err, false};
    std::mutex mutex, recent_mutex;
    bool attached = false, write_failed = false;
    unsigned sequence = 0;
    std::optional<double> head_override, mass_override;
    std::optional<int> iteration_override;
    Session(std::string name, const std::filesystem::path &path) : program(std::move(name)) {
        const auto ticks = std::chrono::high_resolution_clock::now().time_since_epoch().count();
#ifdef _WIN32
        const auto pid = _getpid();
#else
        const auto pid = getpid();
#endif
        run_id = program + "-" + std::to_string(pid) + "-" + std::to_string(ticks);
#ifdef _WIN32
        file = CreateFileW(path.wstring().c_str(), GENERIC_READ | FILE_APPEND_DATA,
                           FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr, OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
        const bool opened = file != INVALID_HANDLE_VALUE;
#else
        file = ::open(path.c_str(), O_CREAT | O_APPEND | O_WRONLY, 0666);
        const bool opened = file >= 0;
#endif
        if (!opened) throw std::runtime_error("Cannot open diagnostics file '" + path.string() + "' for append. Check its parent directory and write permissions.");
    }
    void attach() {
        active = this;
        std::cout.rdbuf(&capture_out); std::cerr.rdbuf(&capture_err); attached = true;
    }
    void detach() {
        if (!attached) return;
        capture_out.finish(); capture_err.finish();
        std::cout.rdbuf(out); std::cerr.rdbuf(err); attached = false; active = nullptr;
    }
    ~Session() {
        detach();
#ifdef _WIN32
        if (file != INVALID_HANDLE_VALUE) CloseHandle(file);
#else
        if (file >= 0) ::close(file);
#endif
    }
    std::string recent_context() {
        std::lock_guard<std::mutex> guard(recent_mutex);
        return recent;
    }
    bool append(const std::string &line) {
#ifdef _WIN32
        OVERLAPPED lock{};
        if (!LockFileEx(file, LOCKFILE_EXCLUSIVE_LOCK, 0, MAXDWORD, MAXDWORD, &lock)) return false;
        DWORD written = 0;
        const bool ok = WriteFile(file, line.data(), static_cast<DWORD>(line.size()), &written, nullptr) && written == line.size();
        UnlockFileEx(file, 0, MAXDWORD, MAXDWORD, &lock);
        return ok;
#else
        int locked;
        do { locked = flock(file, LOCK_EX); } while (locked < 0 && errno == EINTR);
        if (locked < 0) return false;
        std::size_t offset = 0;
        while (offset < line.size()) {
            const auto written = ::write(file, line.data() + offset, line.size() - offset);
            if (written < 0 && errno == EINTR) continue;
            if (written <= 0) break;
            offset += static_cast<std::size_t>(written);
        }
        flock(file, LOCK_UN);
        return offset == line.size();
#endif
    }
    void emit(std::string severity, const std::string &code, const std::string &message,
              const std::string &event = "diagnostic", int exit_code = -1, bool console = false) {
        if (severity == "error" && candidate_depth) severity = "warning";
        std::lock_guard<std::mutex> guard(mutex);
        if (severity == "error") { ++errors; last_error = message; }
        if (severity == "warning") ++warnings;
        json record = {{"schema_version",1}, {"timestamp",timestamp()}, {"program",program},
            {"run_id",run_id}, {"sequence",++sequence}, {"event",event},
            {"severity",severity}, {"code",code}, {"message",message}};
        if (exit_code >= 0) {
            record["exit_code"] = exit_code; record["error_count"] = errors; record["warning_count"] = warnings;
        }
        // Preserve native EPANET section/element/line information in machine fields.
        static const std::regex location(R"(\[([A-Z_ ]+)\](?: element '([^']*)')? line ([0-9]+))");
        std::smatch match;
        if (std::regex_search(message, match, location)) {
            record["section"] = match[1].str(); record["element"] = match[2].str();
            record["line"] = std::stoul(match[3].str());
        }
        static const std::regex network(R"(Network '([^']*)')");
        if (std::regex_search(message, match, network)) record["network"] = match[1].str();
        const std::string line = record.dump(-1, ' ', false, json::error_handler_t::replace) + '\n';
        if (!append(line)) write_failed = true;
        if (console) {
            const std::string text = (severity == "warning" ? "WARNING" : "ERROR") +
                std::string(" [") + program + "][" + code + "]: " + message + '\n';
            err->sputn(text.data(), static_cast<std::streamsize>(text.size())); err->pubsync();
        }
    }
    void observe_noexcept(const std::string &line, bool echo) noexcept {
        try { observe(line, echo); }
        catch (...) { write_failed = true; }
    }
    void observe(const std::string &line, bool echo) {
        {
            std::lock_guard<std::mutex> guard(recent_mutex);
            recent += line + '\n';
            if (recent.size() > 8192) recent.erase(0, recent.size() - 8192);
        }
        static const std::regex diagnostic(R"((^\s*(ERROR|WARNING|Warning|Error)\b)|(\b(ERROR|WARNING)\s*(!|:|\[)))");
        if (std::regex_search(line, diagnostic)) {
            const bool warn = line.find("WARNING") != std::string::npos || line.find("Warning") != std::string::npos;
            static const std::regex tagged(R"((ERROR|WARNING) \[([A-Z0-9_ .-]+)\]\[([A-Z0-9_ .-]+)\])");
            std::smatch match;
            std::string code = warn ? "LEGACY_WARNING" : "LEGACY_ERROR";
            if (std::regex_search(line, match, tagged)) code = match[2].str() + "." + match[3].str();
            emit(warn ? "warning" : "error", code, line, "diagnostic", -1, echo);
        } else if (!line.empty() && line.find_first_not_of(" \t") > 0 &&
                   (line.find("[VALVES]") != std::string::npos || line.find("[OPTIONS]") != std::string::npos ||
                    line.find("[TOPOLOGY]") != std::string::npos || line.find("[EMITTERS]") != std::string::npos)) {
            emit("error", "EPANET.COMPATIBILITY_DETAIL", line);
        }
    }
};
Capture::int_type Capture::overflow(int_type c) {
    if (traits_type::eq_int_type(c, traits_type::eof())) return traits_type::not_eof(c);
    std::lock_guard<std::mutex> guard(capture_mutex_);
    const char value = traits_type::to_char_type(c);
    original_->sputc(value);
    if (value == '\n') { session_.observe_noexcept(line_, echo_); line_.clear(); }
    else line_ += value;
    return c;
}
std::streamsize Capture::xsputn(const char *s, std::streamsize n) {
    std::lock_guard<std::mutex> guard(capture_mutex_);
    original_->sputn(s,n);
    for (std::streamsize i=0; i<n; ++i)
        if (s[i] == '\n') { session_.observe_noexcept(line_, echo_); line_.clear(); }
        else line_ += s[i];
    return n;
}
int Capture::sync() { return original_->pubsync(); }
void Capture::finish() { if (!line_.empty()) { session_.observe_noexcept(line_, echo_); line_.clear(); } }
}
void warning(const std::string &code, const std::string &message) {
    if (active) active->emit("warning",code,message,"diagnostic",-1,true);
    else std::cerr << "WARNING [" << code << "]: " << message << '\n';
}
void error(const std::string &code, const std::string &message) {
    if (active) active->emit("error",code,message,"diagnostic",-1,true);
    else std::cerr << "ERROR [" << code << "]: " << message << '\n';
}
[[noreturn]] void fail_legacy(const char *source, int line) {
    const std::string context = active ? active->recent_context() : "See preceding application diagnostics.";
    throw Error("INPUT_OR_CONFIGURATION", std::string(source) + ":" + std::to_string(line) +
                ": Application stopped. Context:\n" + context);
}
CandidateScope::CandidateScope() { ++candidate_depth; }
CandidateScope::~CandidateScope() { --candidate_depth; }
void apply_solver_overrides(double &head_m, double &mass_kgs, int &iterations) {
    if (!active) return;
    if (active->head_override) head_m = *active->head_override;
    if (active->mass_override) mass_kgs = *active->mass_override;
    if (active->iteration_override) iterations = *active->iteration_override;
}

int run(const char *program, int argc, char **argv, const std::function<int(int,char **)> &application) {
    std::filesystem::path path = "staci-diagnostics.jsonl";
    if (const char *env = std::getenv("STACI_DIAGNOSTICS_FILE")) if (*env) path = env;
    std::vector<char *> arguments{argv[0]};
    std::string argument_error;
    std::optional<double> head_override, mass_override;
    std::optional<int> iteration_override;
    bool show_solver_help = false;
    for (int i=1;i<argc;++i) {
        const std::string argument = argv[i];
        if (argument == "--help" || argument == "-h") show_solver_help = true;
        const auto equals = argument.find('=');
        const std::string key = argument.substr(0, equals);
        if (key == "--head-tolerance-m" || key == "--mass-tolerance-kg-s" || key == "--max-iterations") {
            std::string value;
            if (equals != std::string::npos) value = argument.substr(equals + 1);
            else if (i + 1 < argc && std::string(argv[i+1]).find("--") != 0) value = argv[++i];
            try {
                std::size_t consumed = 0;
                const double number = std::stod(value, &consumed);
                if (consumed != value.size() || !std::isfinite(number) || number <= 0)
                    throw std::invalid_argument("value");
                if (key == "--max-iterations") {
                    if (number != std::floor(number) || number > std::numeric_limits<int>::max())
                        throw std::invalid_argument("integer");
                    iteration_override = static_cast<int>(number);
                } else if (key == "--head-tolerance-m") head_override = number;
                else mass_override = number;
            } catch (const std::exception &) {
                argument_error = key + " requires a positive finite " +
                    (key == "--max-iterations" ? std::string("integer") : std::string("number")) +
                    "; received '" + value + "'.";
            }
        } else if (argument == "--diagnostics-file") {
            if (i+1 == argc || std::string(argv[i+1]).empty()) argument_error = "--diagnostics-file requires a nonempty path.";
            else path = argv[++i];
        } else arguments.push_back(argv[i]);
    }
    arguments.push_back(nullptr);
    try {
        Session session(program,path);
        session.emit("info","RUN_START","Application started.","run_start");
        session.head_override = head_override; session.mass_override = mass_override;
        session.iteration_override = iteration_override;
        session.attach();
        int result = success;
        try {
            if (!argument_error.empty()) throw Error("CLI_ARGUMENT",argument_error);
            result = application(static_cast<int>(arguments.size()-1), arguments.data());
        } catch (const EpanetCompatibilityError &e) {
            error("EPANET.COMPATIBILITY",e.what()); result = input_error;
        } catch (const Error &e) {
            error(e.code,e.what()); result = e.exit_code;
        } catch (const std::invalid_argument &e) {
            error("INVALID_ARGUMENT",e.what()); result = input_error;
        } catch (const std::exception &e) {
            error("INPUT_OR_CALCULATION",e.what()); result = calculation_error;
        } catch (const std::string &e) {
            error("INPUT_OR_CALCULATION",e); result = calculation_error;
        } catch (const char *e) {
            error("INPUT_OR_CALCULATION",e ? e : "Unknown C-string exception."); result = calculation_error;
        } catch (...) {
            error("UNEXPECTED_EXCEPTION","Unknown exception escaped the application."); result = calculation_error;
        }
        if (show_solver_help && result == success)
            std::cout << "\nCommon hydraulic solver options (all four applications):\n"
                "  --head-tolerance-m VALUE    RMS head residual limit [m]; INP default 0.0001 (0.1 mm).\n"
                "  --mass-tolerance-kg-s VALUE RMS continuity residual limit [kg/s]; INP default 1e-8.\n"
                "  --max-iterations INTEGER   Override the input's maximum hydraulic iterations.\n"
                "Command-line values override network/application solver settings.\n";
        session.detach();
        if (result != success && session.errors == 0) session.emit("error","RUN_FAILED","Application returned a failure code; inspect preceding warnings and result files.","diagnostic",-1,true);
        if (session.write_failed) {
            result = calculation_error;
            session.emit("error","DIAGNOSTICS_WRITE","Cannot write the diagnostics file (disk full or write error).","diagnostic",-1,true);
        }
        session.emit("info","RUN_END","Application finished.","run_end",result);
        return session.write_failed ? calculation_error : result;
    } catch (const std::exception &e) {
        std::cerr << "ERROR [" << program << "][DIAGNOSTICS_OPEN]: " << e.what() << '\n';
        return calculation_error;
    }
}
}
