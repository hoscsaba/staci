#include "epanet_steady_quality.h"

#ifdef STACI_USE_EIGEN_SPARSE
#include <Eigen/SparseLU>
#else
#include "umfpack.h"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>

namespace {

constexpr double kFlowEpsilon = 1.0e-12;
constexpr double kVolumeEpsilon = 1.0e-12;

struct DirectedLink {
    std::size_t original = 0;
    std::size_t upstream = 0;
    std::size_t downstream = 0;
    double flow = 0.0;
    double direction = 1.0;
    double travel_time = 0.0;
    double transfer = 1.0;
};

#ifndef STACI_USE_EIGEN_SPARSE
struct SymbolicDeleter {
    void operator()(void *pointer) const noexcept {
        if (pointer != nullptr)
            umfpack_di_free_symbolic(&pointer);
    }
};

struct NumericDeleter {
    void operator()(void *pointer) const noexcept {
        if (pointer != nullptr)
            umfpack_di_free_numeric(&pointer);
    }
};

class SparseFactorization {
public:
    explicit SparseFactorization(const std::vector<std::map<int, double> > &rows) {
        const int n = static_cast<int>(rows.size());
        columns_.resize(rows.size());
        for (int row = 0; row < n; ++row)
            for (const auto &entry : rows[static_cast<std::size_t>(row)])
                if (std::abs(entry.second) > 0.0)
                    columns_[static_cast<std::size_t>(entry.first)].push_back(
                        std::make_pair(row, entry.second));

        ap_.reserve(static_cast<std::size_t>(n) + 1);
        ap_.push_back(0);
        for (const auto &column : columns_) {
            for (const auto &entry : column) {
                ai_.push_back(entry.first);
                ax_.push_back(entry.second);
            }
            ap_.push_back(static_cast<int>(ai_.size()));
        }

        std::array<double, UMFPACK_CONTROL> control;
        std::array<double, UMFPACK_INFO> info;
        umfpack_di_defaults(control.data());
        void *symbolic = nullptr;
        int status = umfpack_di_symbolic(
            n, n, ap_.data(), ai_.data(), ax_.data(), &symbolic,
            control.data(), info.data());
        if (status != UMFPACK_OK || symbolic == nullptr)
            throw std::runtime_error("Steady-quality UMFPACK symbolic factorization failed.");
        symbolic_.reset(symbolic);

        void *numeric = nullptr;
        status = umfpack_di_numeric(
            ap_.data(), ai_.data(), ax_.data(), symbolic_.get(), &numeric,
            control.data(), info.data());
        if (status != UMFPACK_OK || numeric == nullptr)
            throw std::runtime_error("Steady-quality matrix is singular or numerically invalid.");
        numeric_.reset(numeric);
    }

    std::vector<double> solve(const std::vector<double> &rhs) const {
        const int n = static_cast<int>(columns_.size());
        if (static_cast<int>(rhs.size()) != n)
            throw std::invalid_argument("Steady-quality right-hand side has invalid size.");
        std::vector<double> solution(rhs.size(), 0.0);
        std::array<double, UMFPACK_CONTROL> control;
        std::array<double, UMFPACK_INFO> info;
        umfpack_di_defaults(control.data());
        const int status = umfpack_di_solve(
            UMFPACK_A, ap_.data(), ai_.data(), ax_.data(), solution.data(),
            rhs.data(), numeric_.get(), control.data(), info.data());
        if (status != UMFPACK_OK)
            throw std::runtime_error("Steady-quality UMFPACK solve failed.");
        return solution;
    }

private:
    std::vector<std::vector<std::pair<int, double> > > columns_;
    std::vector<int> ap_;
    std::vector<int> ai_;
    std::vector<double> ax_;
    std::unique_ptr<void, SymbolicDeleter> symbolic_;
    std::unique_ptr<void, NumericDeleter> numeric_;
};
#else
class SparseFactorization {
public:
    explicit SparseFactorization(const std::vector<std::map<int, double> > &rows) {
        const int n = static_cast<int>(rows.size());
        std::vector<Eigen::Triplet<double> > entries;
        for (int row = 0; row < n; ++row)
            for (const auto &entry : rows[static_cast<std::size_t>(row)])
                if (entry.second != 0.0)
                    entries.emplace_back(row, entry.first, entry.second);
        matrix_.resize(n, n);
        matrix_.setFromTriplets(entries.begin(), entries.end());
        matrix_.makeCompressed();
        solver_.analyzePattern(matrix_);
        solver_.factorize(matrix_);
        if (solver_.info() != Eigen::Success)
            throw std::runtime_error(
                "Steady-quality matrix is singular or numerically invalid.");
    }

    std::vector<double> solve(const std::vector<double> &rhs) const {
        if (rhs.size() != static_cast<std::size_t>(matrix_.rows()))
            throw std::invalid_argument(
                "Steady-quality right-hand side has invalid size.");
        const Eigen::Map<const Eigen::VectorXd> right_hand_side(
            rhs.data(), static_cast<Eigen::Index>(rhs.size()));
        const Eigen::VectorXd solution = solver_.solve(right_hand_side);
        if (solver_.info() != Eigen::Success || !solution.allFinite())
            throw std::runtime_error("Steady-quality Eigen solve failed.");
        return std::vector<double>(solution.data(),
                                   solution.data() + solution.size());
    }

private:
    Eigen::SparseMatrix<double> matrix_;
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver_;
};
#endif

double exponential_average(double value) {
    if (std::abs(value) < 1.0e-7)
        return 1.0 + value / 2.0 + value * value / 6.0;
    return std::expm1(value) / value;
}

double exponential_average_derivative(double value) {
    if (std::abs(value) < 1.0e-5)
        return 0.5 + value / 3.0 + value * value / 8.0;
    return (std::exp(value) * (value - 1.0) + 1.0) / (value * value);
}

void validate_sensitivity_size(const std::vector<double> &values,
                               std::size_t expected,
                               const char *name) {
    if (!values.empty() && values.size() != expected)
        throw std::invalid_argument(std::string("Steady-quality sensitivity vector has invalid size: ") + name);
}

} // namespace

EpanetSteadyQualityModel::EpanetSteadyQualityModel(
    std::vector<EpanetSteadyQualityNode> nodes,
    std::vector<EpanetSteadyQualityLink> links)
    : nodes_(std::move(nodes)), links_(std::move(links)) {
    if (nodes_.empty())
        throw std::invalid_argument("Steady-quality model requires at least one node.");
    for (const EpanetSteadyQualityNode &node : nodes_) {
        if (!std::isfinite(node.fixed_concentration_kgm3) ||
            node.fixed_concentration_kgm3 < 0.0 ||
            !std::isfinite(node.external_inflow_m3s) ||
            node.external_inflow_m3s < 0.0 ||
            !std::isfinite(node.external_concentration_kgm3) ||
            node.external_concentration_kgm3 < 0.0 ||
            !std::isfinite(node.mass_source_kgs))
            throw std::invalid_argument("Steady-quality node data must be finite and non-negative where applicable.");
    }
    for (const EpanetSteadyQualityLink &link : links_) {
        if (link.from_node >= nodes_.size() || link.to_node >= nodes_.size() ||
            !std::isfinite(link.flow_m3s) || !std::isfinite(link.volume_m3) ||
            link.volume_m3 < 0.0 ||
            !std::isfinite(link.reaction_coefficient_per_s))
            throw std::invalid_argument("Steady-quality link data are invalid.");
    }
}

EpanetSteadyQualityResult EpanetSteadyQualityModel::solve(
    bool solve_age, bool solve_chemical,
    const EpanetSteadyQualitySensitivityInput *sensitivity) const {
    if (!solve_age && !solve_chemical)
        throw std::invalid_argument("Steady-quality solve requires AGE and/or CHEMICAL mode.");
    if (sensitivity != nullptr) {
        validate_sensitivity_size(sensitivity->link_flow_derivative_m3s,
                                  links_.size(), "link flow");
        validate_sensitivity_size(sensitivity->link_volume_derivative_m3,
                                  links_.size(), "link volume");
        validate_sensitivity_size(sensitivity->link_reaction_derivative_per_s,
                                  links_.size(), "link reaction");
        validate_sensitivity_size(sensitivity->fixed_concentration_derivative_kgm3,
                                  nodes_.size(), "fixed concentration");
        validate_sensitivity_size(sensitivity->external_inflow_derivative_m3s,
                                  nodes_.size(), "external inflow");
        validate_sensitivity_size(sensitivity->external_concentration_derivative_kgm3,
                                  nodes_.size(), "external concentration");
        validate_sensitivity_size(sensitivity->mass_source_derivative_kgs,
                                  nodes_.size(), "mass source");
    }

    std::vector<DirectedLink> directed;
    directed.reserve(links_.size());
    std::vector<std::vector<std::size_t> > incoming(nodes_.size());
    EpanetSteadyQualityResult result;
    result.link_travel_time_s.assign(links_.size(), 0.0);
    result.link_transfer_factor.assign(links_.size(), 1.0);
    for (std::size_t index = 0; index < links_.size(); ++index) {
        const EpanetSteadyQualityLink &link = links_[index];
        if (std::abs(link.flow_m3s) <= kFlowEpsilon)
            continue;
        DirectedLink item;
        item.original = index;
        item.direction = link.flow_m3s >= 0.0 ? 1.0 : -1.0;
        item.upstream = link.flow_m3s >= 0.0 ? link.from_node : link.to_node;
        item.downstream = link.flow_m3s >= 0.0 ? link.to_node : link.from_node;
        item.flow = std::abs(link.flow_m3s);
        item.travel_time = link.volume_m3 > kVolumeEpsilon
            ? link.volume_m3 / item.flow : 0.0;
        item.transfer = std::exp(
            link.reaction_coefficient_per_s * item.travel_time);
        result.link_travel_time_s[index] = item.travel_time;
        result.link_transfer_factor[index] = item.transfer;
        incoming[item.downstream].push_back(directed.size());
        directed.push_back(item);
    }

    const auto value_or_zero = [](const std::vector<double> &values,
                                  std::size_t index) {
        return values.empty() ? 0.0 : values[index];
    };

    const auto link_derivatives = [&](const DirectedLink &item,
                                      double &dq, double &dtau,
                                      double &dk, double &dbeta) {
        dq = dtau = dk = dbeta = 0.0;
        if (sensitivity == nullptr)
            return;
        const std::size_t index = item.original;
        const double signed_dq = value_or_zero(
            sensitivity->link_flow_derivative_m3s, index);
        dq = item.direction * signed_dq;
        const double dvolume = value_or_zero(
            sensitivity->link_volume_derivative_m3, index);
        dk = value_or_zero(
            sensitivity->link_reaction_derivative_per_s, index);
        if (links_[index].volume_m3 > kVolumeEpsilon)
            dtau = dvolume / item.flow - item.travel_time * dq / item.flow;
        const double dx = dk * item.travel_time +
            links_[index].reaction_coefficient_per_s * dtau;
        dbeta = item.transfer * dx;
    };

    if (solve_age) {
        std::vector<std::map<int, double> > matrix(nodes_.size());
        std::vector<double> rhs(nodes_.size(), 0.0);
        for (std::size_t node = 0; node < nodes_.size(); ++node) {
            if (nodes_[node].fixed_age) {
                matrix[node][static_cast<int>(node)] = 1.0;
                continue;
            }
            double total = nodes_[node].external_inflow_m3s;
            rhs[node] += nodes_[node].external_inflow_m3s * 0.0;
            for (std::size_t directed_index : incoming[node]) {
                const DirectedLink &item = directed[directed_index];
                total += item.flow;
                matrix[node][static_cast<int>(item.upstream)] -= item.flow;
                rhs[node] += item.flow * item.travel_time;
            }
            if (total <= kFlowEpsilon)
                throw std::runtime_error("Steady water age has an unfixed node without inflow: " + nodes_[node].id);
            matrix[node][static_cast<int>(node)] += total;
        }
        SparseFactorization factorization(matrix);
        result.node_age_s = factorization.solve(rhs);
        result.link_average_age_s.assign(links_.size(), 0.0);
        for (const DirectedLink &item : directed)
            result.link_average_age_s[item.original] =
                result.node_age_s[item.upstream] + 0.5 * item.travel_time;

        if (sensitivity != nullptr) {
            std::vector<double> derivative_rhs(nodes_.size(), 0.0);
            for (std::size_t node = 0; node < nodes_.size(); ++node) {
                if (nodes_[node].fixed_age)
                    continue;
                const double dexternal = value_or_zero(
                    sensitivity->external_inflow_derivative_m3s, node);
                derivative_rhs[node] += dexternal *
                    (0.0 - result.node_age_s[node]);
                for (std::size_t directed_index : incoming[node]) {
                    const DirectedLink &item = directed[directed_index];
                    double dq, dtau, dk, dbeta;
                    link_derivatives(item, dq, dtau, dk, dbeta);
                    derivative_rhs[node] += dq *
                        (result.node_age_s[item.upstream] + item.travel_time -
                         result.node_age_s[node]);
                    derivative_rhs[node] += item.flow * dtau;
                }
            }
            result.node_age_sensitivity = factorization.solve(derivative_rhs);
            result.link_average_age_sensitivity.assign(links_.size(), 0.0);
            for (const DirectedLink &item : directed) {
                double dq, dtau, dk, dbeta;
                link_derivatives(item, dq, dtau, dk, dbeta);
                result.link_average_age_sensitivity[item.original] =
                    result.node_age_sensitivity[item.upstream] + 0.5 * dtau;
            }
        }
    }

    if (solve_chemical) {
        std::vector<std::map<int, double> > matrix(nodes_.size());
        std::vector<double> rhs(nodes_.size(), 0.0);
        for (std::size_t node = 0; node < nodes_.size(); ++node) {
            if (nodes_[node].fixed_chemical) {
                matrix[node][static_cast<int>(node)] = 1.0;
                rhs[node] = nodes_[node].fixed_concentration_kgm3;
                continue;
            }
            double total = nodes_[node].external_inflow_m3s;
            rhs[node] = nodes_[node].mass_source_kgs +
                nodes_[node].external_inflow_m3s *
                nodes_[node].external_concentration_kgm3;
            for (std::size_t directed_index : incoming[node]) {
                const DirectedLink &item = directed[directed_index];
                total += item.flow;
                matrix[node][static_cast<int>(item.upstream)] -=
                    item.flow * item.transfer;
            }
            if (total <= kFlowEpsilon)
                throw std::runtime_error("Steady chemical quality has an unfixed node without inflow: " + nodes_[node].id);
            matrix[node][static_cast<int>(node)] += total;
        }
        SparseFactorization factorization(matrix);
        result.node_concentration_kgm3 = factorization.solve(rhs);
        result.link_average_concentration_kgm3.assign(links_.size(), 0.0);
        for (const DirectedLink &item : directed) {
            const double x = links_[item.original].reaction_coefficient_per_s *
                item.travel_time;
            result.link_average_concentration_kgm3[item.original] =
                result.node_concentration_kgm3[item.upstream] *
                exponential_average(x);
        }

        if (sensitivity != nullptr) {
            std::vector<double> derivative_rhs(nodes_.size(), 0.0);
            for (std::size_t node = 0; node < nodes_.size(); ++node) {
                if (nodes_[node].fixed_chemical) {
                    derivative_rhs[node] = value_or_zero(
                        sensitivity->fixed_concentration_derivative_kgm3, node);
                    continue;
                }
                const double dexternal = value_or_zero(
                    sensitivity->external_inflow_derivative_m3s, node);
                const double dexternal_concentration = value_or_zero(
                    sensitivity->external_concentration_derivative_kgm3, node);
                derivative_rhs[node] += value_or_zero(
                    sensitivity->mass_source_derivative_kgs, node);
                derivative_rhs[node] += dexternal *
                    (nodes_[node].external_concentration_kgm3 -
                     result.node_concentration_kgm3[node]);
                derivative_rhs[node] += nodes_[node].external_inflow_m3s *
                    dexternal_concentration;
                for (std::size_t directed_index : incoming[node]) {
                    const DirectedLink &item = directed[directed_index];
                    double dq, dtau, dk, dbeta;
                    link_derivatives(item, dq, dtau, dk, dbeta);
                    derivative_rhs[node] += dq *
                        (item.transfer *
                             result.node_concentration_kgm3[item.upstream] -
                         result.node_concentration_kgm3[node]);
                    derivative_rhs[node] += item.flow * dbeta *
                        result.node_concentration_kgm3[item.upstream];
                }
            }
            result.node_concentration_sensitivity =
                factorization.solve(derivative_rhs);
            result.link_average_concentration_sensitivity.assign(
                links_.size(), 0.0);
            for (const DirectedLink &item : directed) {
                double dq, dtau, dk, dbeta;
                link_derivatives(item, dq, dtau, dk, dbeta);
                const double k = links_[item.original].reaction_coefficient_per_s;
                const double x = k * item.travel_time;
                const double dx = dk * item.travel_time + k * dtau;
                result.link_average_concentration_sensitivity[item.original] =
                    result.node_concentration_sensitivity[item.upstream] *
                        exponential_average(x) +
                    result.node_concentration_kgm3[item.upstream] *
                        exponential_average_derivative(x) * dx;
            }
        }
    }

    return result;
}
