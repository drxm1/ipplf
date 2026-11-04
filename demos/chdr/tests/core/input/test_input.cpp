/** @file test_input.cpp
 * @brief Check strict common input and one independently evaluated SI normalization.
 * @details readRunInput() reads temporary YAML files. It must throw std::invalid_argument for
 * 8 malformed run or observation entries and for 7 invalid unit declarations, and
 * booleanParameter() must throw for 3 non-boolean values. applicationInput() must reject 3
 * incomplete or mixed command-line selectors, and a missing file and unclosed YAML must throw.
 * The valid input checks partial schedule overrides. For \f$L_0 = 1\f$ mm and
 * \f$\rho_0 = 8.8541878188\times10^{-9}\ \mathrm{C\,m^{-3}}\f$, the SI factors must equal fixed
 * values (EXPECT_DOUBLE_EQ, or absolute tolerances such as \f$10^{-26}\f$ s on
 * \f$t_0 = 3.3356409519815207\times10^{-12}\f$ s). An SI duration of 8 steps of 0.25 must give
 * 8 steps; 8.5 steps and the largest double must throw.
 *
 * Reference: refusal; closed-form SI factors. CTest entry: chdr.input.unit.serial (one rank,
 * timeout 30 s).
 */
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <unistd.h>

#include "core/input/Application.h"

namespace {
    using namespace chdr::input;

    /// Exercise the actual file reader and remove its temporary input on either outcome.
    RunInput parseYaml(const std::string& text) {
        auto path = (std::filesystem::temp_directory_path() / "chdr-input-XXXXXX").string();
        const int descriptor = mkstemp(path.data());
        if (descriptor < 0)
            throw std::runtime_error("Cannot create temporary YAML input");
        close(descriptor);
        try {
            std::ofstream output(path);
            output << text;
            output.close();
            auto result = readRunInput(path);
            std::filesystem::remove(path);
            return result;
        } catch (...) {
            std::filesystem::remove(path);
            throw;
        }
    }

    /// Supply only the required common keys so each malformed suffix isolates one rule.
    std::string normalizedInput() {
        return "schema_version: 1\ncase: coarse\nunits: {system: normalized}\n";
    }

    /// Check partial schedule overrides and application-specific validation through the public API.
    void checkValidInput() {
        auto input = parseYaml(normalizedInput()
            + "parameters: {smoke: true}\n"
              "run: {steps: 7}\nobservations: {fields: {enabled: false, first: 1}}\n");
        EXPECT_EQ(resolvedSteps(input, 0.25, 10), 7);
        EXPECT_TRUE(booleanParameter(input, "smoke", false));
        EXPECT_FALSE(booleanParameter(input, "omitted", false));
        EXPECT_NO_THROW(validateParameters(input, {"smoke"}));
        EXPECT_THROW(validateParameters(input, {}), std::invalid_argument);
        const auto schedule = observationSchedule(input, "fields", {4, 0, 20, true});
        EXPECT_EQ(schedule.first_m, 1);
        EXPECT_EQ(schedule.every_m, 4);
        EXPECT_EQ(schedule.last_m, 20);
        EXPECT_TRUE(schedule.includeFinal_m);
        EXPECT_FALSE(observationEnabled(input, "fields"));
        EXPECT_NO_THROW(validateObservationNames(input, {"fields"}));
        EXPECT_THROW(validateObservationNames(input, {"energy"}), std::invalid_argument);
    }

    /// Catch silent key typos, coercion of unsuitable types, and invalid disabled schedules.
    void checkInvalidInput() {
        for (const std::string suffix :
             {"rnu: {steps: 1}\n", "run: {steps: 1, duration: 1}\n", "run: {steps: -1}\n",
              "run: {steps: 1.5}\n", "run: {steps: 2147483648}\n", "run: {steps: [1, 2]}\n",
              "observations: {fields: {ever: 2}}\n",
              "observations: {fields: {enabled: false, every: 0}}\n"})
            EXPECT_THROW(parseYaml(normalizedInput() + suffix), std::invalid_argument) << suffix;
        for (const std::string value : {"yes", "2", "[true, false]"}) {
            auto input = parseYaml(normalizedInput() + "parameters: {smoke: " + value + "}\n");
            EXPECT_THROW(booleanParameter(input, "smoke", false), std::invalid_argument);
        }
    }

    /// Reject incomplete or mixed selectors before accessing a file or the legacy parser.
    void checkInputSelection() {
        const auto legacy = [](const std::vector<std::string>&) {
            RunInput input;
            input.caseName_m = "legacy";
            return input;
        };
        EXPECT_EQ(applicationInput({"coarse"}, legacy).caseName_m, "legacy");
        for (const std::vector<std::string>& args : {std::vector<std::string>{"--input"},
                                                     {"--input", "file.yaml", "--steps", "1"},
                                                     {"coarse", "--input", "file.yaml"}})
            EXPECT_THROW(applicationInput(args, legacy), std::invalid_argument);
    }

    /// Exercise missing files and parser syntax failures without launching an application.
    void checkInputFiles() {
        auto path = (std::filesystem::temp_directory_path() / "chdr-missing-XXXXXX").string();
        const int descriptor = mkstemp(path.data());
        ASSERT_GE(descriptor, 0);
        close(descriptor);
        std::filesystem::remove(path);
        EXPECT_THROW(readRunInput(path), std::exception);
        EXPECT_THROW(parseYaml(normalizedInput() + "parameters: [unclosed\n"), std::exception);
    }

    /// Valid input with partial overrides, 8 malformed entries, 3 non-boolean parameter values,
    /// 3 invalid selector lists, a missing file and unclosed YAML.
    TEST(RunInput, StrictSchemaAndOverrides) {
        checkValidInput();
        checkInvalidInput();
        checkInputSelection();
        checkInputFiles();
    }

    /** @brief A millimetre scale with one volt per metre for unit code electric field.
     * @details The fixed values below evaluate the normalization derived in UnitSystem,
     * \f$t_0 = L_0/c\f$, \f$E_0 = \rho_0L_0/\varepsilon_0\f$, \f$B_0 = E_0/c\f$,
     * \f$\phi_0 = E_0L_0\f$, \f$A_0 = B_0L_0\f$ and \f$U_0 = \varepsilon_0E_0^2L_0^3\f$
     * [fitzpatrick2014ClassicalElectromagnetism, Eqs. (1.1), (1.3), (1.11)-(1.16), (1.109)],
     * with the CODATA 2022 values of \f$c\f$ and \f$\varepsilon_0\f$
     * [mohr2025Codata2022, Table XXXIII, p. 025002-45], not the implementation's factor methods.
     * They detect missing powers of length and confusion between scalar/vector potentials.
     */
    void checkSIFactors(const UnitSystem& units) {
        EXPECT_DOUBLE_EQ(units.lengthFactor(), 0.001);
        EXPECT_NEAR(units.timeFactor(), 3.3356409519815207e-12, 1e-26);
        EXPECT_DOUBLE_EQ(units.electricFactor(), 1.0);
        EXPECT_NEAR(units.magneticFactor(), 3.3356409519815204e-9, 1e-23);
        EXPECT_DOUBLE_EQ(units.scalarPotentialFactor(), 0.001);
        EXPECT_NEAR(units.vectorPotentialFactor(), 3.3356409519815207e-12, 1e-26);
        EXPECT_NEAR(units.energyFactor(), 8.8541878188e-21, 3e-35);
        EXPECT_DOUBLE_EQ(units.velocityFactor(), 299792458.0);
        EXPECT_DOUBLE_EQ(units.frequencyFactor(), 299792458000.0);
        EXPECT_DOUBLE_EQ(UnitSystem{}.energyFactor(), 1.0);
    }

    /// Check explicit scales and conversion overflow before any application allocates fields.
    void checkInvalidSI() {
        const std::string prefix = "schema_version: 1\ncase: coarse\nunits: ";
        for (const std::string units : {"{system: si}", "{system: si, length_m: 0.001}",
                                        "{system: si, length_m: -1, charge_density_C_m3: 1}",
                                        "{system: si, length_m: 1e300, charge_density_C_m3: 1e300}",
                                        "{system: si, length_m: wrong, charge_density_C_m3: 1}",
                                        "{system: cgs}", "{system: normalized, length_m: 1}"})
            EXPECT_THROW(parseYaml(prefix + units + "\n"), std::invalid_argument) << units;
    }

    /// SI factors for 1 mm and \f$\rho_0 = 8.8541878188\times10^{-9}\ \mathrm{C\,m^{-3}}\f$, the
    /// duration-to-step conversion with step 0.25, and 7 invalid unit declarations.
    TEST(RunInput, SIUnitsAndDuration) {
        auto input = parseYaml(
            "schema_version: 1\ncase: coarse\n"
            "units: {system: si, length_m: 0.001, charge_density_C_m3: 8.8541878188e-9}\n"
            "run: {duration: 6.6712819039630414e-12}\n");
        checkSIFactors(input.units_m);
        EXPECT_EQ(resolvedSteps(input, 0.25, 20), 8);
        input.duration_m = 8 * 0.25 * input.units_m.timeFactor();
        EXPECT_EQ(resolvedSteps(input, 0.25, 20), 8);
        input.duration_m = 8.5 * 0.25 * input.units_m.timeFactor();
        EXPECT_THROW(resolvedSteps(input, 0.25, 20), std::invalid_argument);
        input.duration_m = std::numeric_limits<double>::max();
        EXPECT_THROW(resolvedSteps(input, 0.25, 20), std::invalid_argument);
        checkInvalidSI();
    }
}  // namespace
