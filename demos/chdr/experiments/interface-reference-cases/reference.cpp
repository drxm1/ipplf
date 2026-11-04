/** @file reference.cpp
 * @brief Evaluate independently implemented planar-interface references without a simulation.
 * @details Inputs use normalized solver units. The continuum and aligned-lattice formulas
 * are documented in references.org and implemented in InterfaceOracles.h. The fixed
 * component is E_y for TE and H_y for TM; no measurement is read or changed here.
 */
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "tests/apps/planar-interface/InterfaceOracles.h"

#if __has_include("ChdrSourceIdentity.h")
#include "ChdrSourceIdentity.h"
#endif
#ifndef CHDR_SOURCE_REVISION
#define CHDR_SOURCE_REVISION "unrecorded"
#endif
#ifndef CHDR_SOURCE_DIRTY
#define CHDR_SOURCE_DIRTY "unrecorded"
#endif

namespace {
    namespace plane = chdr::apps::planarInterface;

    /// Explicit reference inputs; no application preset or numerical setup is inferred.
    struct Request {
        std::string polarization_m, frequencies_m;
        double epsilonIncident_m = 0, epsilonTransmitted_m = 0, dt_m = 0;
        double kz_m = std::numeric_limits<double>::quiet_NaN();
        std::array<double, 3> spacing_m{};
        bool continuumOnly_m = false;
    };

    /// Applicability belongs to the reference, independently of measurement validity.
    struct Answer {
        plane::LatticeCoefficients coefficients_m{};
        std::string incidentRoot_m = "not_evaluated", transmittedRoot_m = "not_evaluated";
        std::string reason_m;
        bool valid_m         = false;
        double reflectance_m = 0, transmittance_m = 0, angleDegrees_m = 0;
        double transmissionWeight_m = 0;
    };

    /// Parse a finite scalar, rejecting trailing material instead of accepting a prefix.
    double number(const std::string& text) {
        std::istringstream input(text);
        double value;
        std::string extra;
        if (!(input >> value) || (input >> extra) || !std::isfinite(value))
            throw std::invalid_argument("Expected a finite number: " + text);
        return value;
    }

    /// Three spacings are needed because the reference uses the recorded Cartesian mesh.
    std::array<double, 3> spacing(const std::string& text) {
        std::istringstream input(text);
        std::array<double, 3> result{};
        std::string item;
        unsigned count = 0;
        while (std::getline(input, item, ',')) {
            if (count == 3)
                throw std::invalid_argument("--spacing expects hx,hy,hz");
            result[count++] = number(item);
        }
        if (count != 3 || text.back() == ',')
            throw std::invalid_argument("--spacing expects hx,hy,hz");
        return result;
    }

    /// Named options keep the incident/transmitted material ordering explicit.
    void assign(Request& result, const std::string& key, const std::string& value) {
        if (key == "--polarization")
            result.polarization_m = value;
        else if (key == "--epsilon-incident")
            result.epsilonIncident_m = number(value);
        else if (key == "--epsilon-transmitted")
            result.epsilonTransmitted_m = number(value);
        else if (key == "--spacing")
            result.spacing_m = spacing(value);
        else if (key == "--dt")
            result.dt_m = number(value);
        else if (key == "--kz")
            result.kz_m = number(value);
        else if (key == "--frequencies")
            result.frequencies_m = value;
        else
            throw std::invalid_argument("Unknown option: " + key);
    }

    /// Reference input validation does not impose a simulation stability or accuracy claim.
    void validate(const Request& request) {
        if (request.polarization_m != "TE" && request.polarization_m != "TM")
            throw std::invalid_argument("--polarization must be TE or TM");
        if (!(request.epsilonIncident_m > 0 && request.epsilonTransmitted_m > 0 && request.dt_m > 0)
            || !std::isfinite(request.kz_m))
            throw std::invalid_argument("Positive permittivities/dt and finite kz are required");
        for (double value : request.spacing_m)
            if (!(value > 0))
                throw std::invalid_argument("All three spacings must be positive");
        if (request.frequencies_m.empty())
            throw std::invalid_argument("--frequencies FILE is required; use - for stdin");
    }

    Request parse(int argc, char** argv) {
        Request result;
        for (int index = 1; index < argc; ++index) {
            const std::string key(argv[index]);
            if (key == "--continuum-only")
                result.continuumOnly_m = true;
            else {
                if (++index == argc)
                    throw std::invalid_argument("Missing value for " + key);
                assign(result, key, argv[index]);
            }
        }
        validate(result);
        return result;
    }

    /// The file contains one positive normalized angular frequency per nonempty line.
    std::vector<double> readFrequencies(std::istream& input) {
        std::vector<double> result;
        std::string line;
        while (std::getline(input, line)) {
            if (line.find_first_not_of(" \t\r") == std::string::npos)
                continue;
            const double omega = number(line);
            if (!(omega > 0))
                throw std::invalid_argument("Frequencies must be positive");
            result.push_back(omega);
        }
        if (input.bad() || result.empty())
            throw std::invalid_argument("Cannot read a nonempty frequency vector");
        return result;
    }

    std::vector<double> frequencies(const Request& request) {
        if (request.frequencies_m == "-")
            return readFrequencies(std::cin);
        std::ifstream input(request.frequencies_m);
        if (!input)
            throw std::runtime_error("Cannot open " + request.frequencies_m);
        return readFrequencies(input);
    }

    plane::Polarization polarization(const Request& request) {
        return request.polarization_m == "TE" ? plane::Polarization::TE : plane::Polarization::TM;
    }

    bool finite(plane::Complex value) {
        return std::isfinite(value.real()) && std::isfinite(value.imag());
    }

    /// A zero normal root is labelled separately from an ordinary propagating mode.
    std::string continuumRoot(plane::Complex normal) {
        if (!finite(normal))
            return "nonfinite";
        if (normal.imag() > 0)
            return "evanescent";
        return normal.real() > 0 ? "propagating" : "cutoff";
    }

    /// Classify the strict interior branch; latticeNormal also labels its endpoints real.
    std::string latticeRoot(const plane::LatticeNormal& normal, double h) {
        if (!finite(normal.kx))
            return "nonfinite";
        if (normal.root == plane::LatticeRoot::Evanescent)
            return "evanescent";
        if (normal.root == plane::LatticeRoot::Sawtooth)
            return "sawtooth";
        return normal.kx.real() > 0 && normal.kx.real() * h < std::numbers::pi ? "propagating"
                                                                               : "cutoff";
    }

    /** @brief Continuum normal-flux fractions in the measured component basis.
     * @details From the Poynting vector [ellingson2020ElectromagneticsVol2, Eq. (3.98)],
     * TE has \f$T=\operatorname{Re}k_2/k_1\,|t_E|^2\f$, while fixed H_y TM has
     * an additional \f$\varepsilon_1/\varepsilon_2\f$ factor; see references.org.
     * An evanescent transmitted wave in a lossless half-space has zero normal flux.
     */
    void continuumPower(Answer& answer, const Request& request) {
        const auto& c               = answer.coefficients_m;
        const double weight         = request.polarization_m == "TM"
                                          ? request.epsilonIncident_m / request.epsilonTransmitted_m
                                          : 1;
        answer.reflectance_m        = std::norm(c.r);
        answer.transmissionWeight_m = c.k2.real() == 0 ? 0 : weight * c.k2.real() / c.k1.real();
        answer.transmittance_m      = answer.transmissionWeight_m * std::norm(c.t);
        answer.angleDegrees_m =
            std::atan2(std::abs(request.kz_m), c.k1.real()) * 180 / std::numbers::pi;
        if (!std::isfinite(answer.reflectance_m) || !std::isfinite(answer.transmittance_m)) {
            answer.valid_m  = false;
            answer.reason_m = "nonfinite_power";
        }
    }

    Answer continuum(const Request& request, double omega) {
        Answer result;
        result.coefficients_m =
            plane::fresnelCoefficients(polarization(request), request.epsilonIncident_m,
                                       request.epsilonTransmitted_m, 0, request.kz_m, omega);
        result.incidentRoot_m    = continuumRoot(result.coefficients_m.k1);
        result.transmittedRoot_m = continuumRoot(result.coefficients_m.k2);
        if (result.incidentRoot_m != "propagating")
            result.reason_m = "incident_not_propagating";
        else if (result.transmittedRoot_m == "cutoff")
            result.reason_m = "transmitted_cutoff";
        else if (!finite(result.coefficients_m.r) || !finite(result.coefficients_m.t))
            result.reason_m = "nonfinite_coefficients";
        result.valid_m = result.reason_m.empty();
        if (result.valid_m)
            continuumPower(result, request);
        return result;
    }

    /// The aligned oracle is intentionally disabled for a continuum-only request.
    std::string latticeDomain(const Request& request, double omega) {
        if (request.continuumOnly_m)
            return "not_requested";
        if (!(omega * request.dt_m < std::numbers::pi))
            return "temporal_nyquist";
        if (!(std::abs(request.kz_m) * request.spacing_m[2] < std::numbers::pi))
            return "transverse_nyquist";
        return "";
    }

    void latticeRoots(Answer& result, const Request& request, double omega) {
        const auto first  = plane::latticeNormal(omega, 0, request.kz_m, request.epsilonIncident_m,
                                                 request.spacing_m, request.dt_m);
        const auto second = plane::latticeNormal(
            omega, 0, request.kz_m, request.epsilonTransmitted_m, request.spacing_m, request.dt_m);
        result.incidentRoot_m    = latticeRoot(first, request.spacing_m[0]);
        result.transmittedRoot_m = latticeRoot(second, request.spacing_m[0]);
        if (result.incidentRoot_m != "propagating")
            result.reason_m = "incident_not_propagating";
        else if (result.transmittedRoot_m == "cutoff")
            result.reason_m = "transmitted_cutoff";
        else if (result.transmittedRoot_m == "sawtooth")
            result.reason_m = "transmitted_sawtooth";
    }

    Answer lattice(const Request& request, double omega) {
        Answer result;
        result.reason_m = latticeDomain(request, omega);
        if (!result.reason_m.empty())
            return result;
        latticeRoots(result, request, omega);
        if (!result.reason_m.empty())
            return result;
        result.coefficients_m = plane::latticeCoefficients(
            polarization(request), request.epsilonIncident_m, request.epsilonTransmitted_m,
            request.spacing_m, request.dt_m, 0, request.kz_m, omega);
        if (!finite(result.coefficients_m.r) || !finite(result.coefficients_m.t))
            result.reason_m = "nonfinite_coefficients";
        result.valid_m = result.reason_m.empty();
        return result;
    }

    void numberColumn(double value, bool valid) {
        std::cout << ',';
        if (valid)
            std::cout << value;
    }

    void complexColumns(plane::Complex value, bool valid) {
        numberColumn(value.real(), valid);
        numberColumn(value.imag(), valid);
    }

    void answerColumns(const Answer& answer, bool power) {
        std::cout << ',' << answer.valid_m << ',' << answer.reason_m << ',' << answer.incidentRoot_m
                  << ',' << answer.transmittedRoot_m;
        const auto& c = answer.coefficients_m;
        for (auto value : {c.k1, c.k2, c.r, c.t})
            complexColumns(value, answer.valid_m);
        if (!power)
            return;
        for (double value : {answer.reflectance_m, answer.transmittance_m, answer.angleDegrees_m,
                             answer.transmissionWeight_m})
            numberColumn(value, answer.valid_m);
    }

    void header() {
        std::cout
            << "omega_code,continuum_valid,continuum_reason,continuum_incident_root,"
               "continuum_transmitted_root,continuum_k1_real,continuum_k1_imag,continuum_k2_real,"
               "continuum_k2_imag,continuum_r_real,continuum_r_imag,continuum_t_real,continuum_t_"
               "imag,"
               "continuum_R,continuum_T,continuum_angle_degrees,continuum_transmission_weight,"
               "lattice_valid,lattice_reason,lattice_incident_root,"
               "lattice_transmitted_root,lattice_k1_real,lattice_k1_imag,lattice_k2_real,"
               "lattice_k2_imag,lattice_r_real,lattice_r_imag,lattice_t_real,lattice_t_imag\n";
    }

    void write(const Request& request) {
        const auto bins = frequencies(request);
        std::cout << std::setprecision(17);
        header();
        for (double omega : bins) {
            std::cout << omega;
            answerColumns(continuum(request, omega), true);
            answerColumns(lattice(request, omega), false);
            std::cout << '\n';
        }
        if (!std::cout)
            throw std::runtime_error("Failed to write reference CSV");
    }

    void identity() {
        std::cout << "{\"program\":\"interface-reference\",\"reference_schema_version\":1,"
                     "\"source_revision\":\""
                  << CHDR_SOURCE_REVISION << "\",\"source_dirty\":\"" << CHDR_SOURCE_DIRTY
                  << "\"}\n";
    }
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--identity")
            identity();
        else if (argc == 2 && std::string(argv[1]) == "--help")
            std::cout
                << "interface-reference --polarization TE|TM --epsilon-incident E1 "
                   "--epsilon-transmitted E2 --spacing hx,hy,hz --dt DT --kz KZ "
                   "--frequencies FILE [--continuum-only]\nFILE: one normalized angular frequency "
                   "per line; - reads stdin. --identity reports build provenance.\n";
        else
            write(parse(argc, argv));
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "interface-reference: " << error.what() << '\n';
        return 2;
    }
}
