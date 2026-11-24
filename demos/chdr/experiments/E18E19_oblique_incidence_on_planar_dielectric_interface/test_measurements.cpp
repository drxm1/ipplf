/** @file test_measurements.cpp
 * @brief Check E18/E19 reductions on prescribed fields, without evolving potentials.
 */
#include "Ippl.h"

#include <Kokkos_MathematicalConstants.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "core/solver/BoxSetup.h"
#include "experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/E18E19Measurements.h"

namespace {
    using namespace chdr;
    using namespace chdr::experiments::e18e19;
    constexpr double KX          = 2 * Kokkos::numbers::pi;
    constexpr double Q           = 4 * Kokkos::numbers::pi;
    constexpr double INTERFACE_Z = 0.5;

    /// Keep every numerical check active in release builds.
    void expectNear(double actual, double expected, const std::string& name) {
        if (!std::isfinite(actual) || std::abs(actual - expected) > 2e-12)
            throw std::runtime_error(name + ": actual " + std::to_string(actual) + ", expected "
                                     + std::to_string(expected));
    }

    /// Supply an independently stated two-material array, with its face between k=3 and k=4.
    struct DielectricMap {
        KOKKOS_INLINE_FUNCTION double operator()(material::CellIndex cell, material::Position,
                                                 bool& marked) const {
            marked = cell[2] >= 4;
            return marked ? 4.0 : 1.0;
        }
    };

    /// Prescribe owned vector values while retaining NaNs in every ghost cell.
    template <class Field, class Sample>
    void plant(Field& field, Sample sample) {
        field            = std::numeric_limits<double>::quiet_NaN();
        const auto view  = field.getView();
        const auto first = field.getLayout().getLocalNDIndex().first();
        const int ghost  = field.getNghost();
        Kokkos::parallel_for(
            "E18/E19 prescribed measurement fixture", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                for (unsigned component = 0; component < 3; ++component)
                    view(i, j, k)[component] = sample(i + first[0] - ghost, j + first[1] - ghost,
                                                      k + first[2] - ghost, component);
            });
    }

    /** @brief Check native positions and inclusion of the magnetic sample on the interface.
     * @details One nonzero sample gives its own negative complex exponential exactly.
     * The adjacent dielectric cell has a much larger value and must be excluded.
     * Positions below are stated from the 8 by 2 by 8 grid independently of the
     * measurement's LevelContract lookup. Both rank contributions are reduced first.
     */
    template <class Field>
    void checkSamplePositions(Field& field, LevelContract contract, Quantity quantity) {
        plant(
            field, KOKKOS_LAMBDA(int i, int j, int k, unsigned component) {
                return i == 0 && j == 0 && component == 1 ? (k == 3   ? 1.0
                                                             : k == 4 ? 99.0
                                                                      : 0.0)
                                                          : 0.0;
            });
        std::uint64_t invalid = 0;
        auto sum = measurementDetail::localFourierSum(field, contract, quantity, KX, 1.7,
                                                      INTERFACE_Z, invalid);
        MPI_Allreduce(MPI_IN_PLACE, sum.data(), 2, MPI_DOUBLE, MPI_SUM, field.getLayout().comm);
        const bool magnetic = quantity == Quantity::Magnetic;
        const double x      = magnetic ? 1.0 / 8 : 1.0 / 16;
        const double z      = magnetic ? 0.5 : 7.0 / 16;
        expectNear(sum[0], std::cos(KX * x + 1.7 * z), "native Fourier real part");
        expectNear(sum[1], -std::sin(KX * x + 1.7 * z), "native Fourier imaginary part");
        expectNear(static_cast<double>(invalid), 0, "finite owned impulse");
    }

    /// A full tangential period cancels globally even when each MPI slab has nonzero magnitude.
    template <class Field>
    void checkGlobalCancellation(Field& field, LevelContract contract) {
        plant(field, KOKKOS_LAMBDA(int, int, int, unsigned) { return 1.0; });
        std::array<double, 12> wavenumbers{};
        const auto magnitudes = measureFourierMagnitude(field, contract, Quantity::Electric, KX,
                                                        wavenumbers, INTERFACE_Z, false);
        for (double magnitude : magnitudes)
            expectNear(magnitude, 0, "global tangential cancellation");
    }

    /// Integer periods give magnitude N/2 for the selected traveling wave and zero for its partner.
    template <class Field>
    void checkDirection(Field& field, LevelContract contract, bool reflected) {
        const int sign = reflected ? -1 : 1;
        plant(
            field, KOKKOS_LAMBDA(int i, int, int k, unsigned component) {
                if (k >= 4)
                    return 99.0;
                const double phase = KX * (i + 0.5) / 8 + sign * Q * (k + 0.5) / 8;
                return component == 1 ? Kokkos::cos(phase) : 0.0;
            });
        std::array<double, 12> wavenumbers{};
        wavenumbers.fill(Q);
        for (bool selectReflected : {false, true}) {
            const auto magnitudes = measureFourierMagnitude(
                field, contract, Quantity::Electric, KX, wavenumbers, INTERFACE_Z, selectReflected);
            for (double magnitude : magnitudes)
                expectNear(magnitude, reflected == selectReflected ? 32 : 0, "wave direction");
        }
    }

    /** @brief Check both energy quadratures against closed cell counts.
     * @details The prescribed E=(1,2,3), B=(4,5,6) give squared sums 14 and 77.
     * Each side contains 64 cells. Its 16 positive-z face samples use harmonic
     * permittivity 1.6, on the material face and on the periodic second face.
     * The common volume factor is half a cell volume, as defined in E18/E19.meas02.
     */
    template <class Field>
    void checkEnergy(Field& electric, Field& magnetic, const material::MaterialOwner& owner) {
        plant(electric, KOKKOS_LAMBDA(int, int, int, unsigned c) { return double(c + 1); });
        plant(magnetic, KOKKOS_LAMBDA(int, int, int, unsigned c) { return double(c + 4); });
        const auto te           = measureSideEnergy(electric, magnetic, owner, false, INTERFACE_Z);
        const auto tm           = measureSideEnergy(electric, magnetic, owner, true, INTERFACE_Z);
        constexpr double FACTOR = 0.5 / (8 * 8 * 8);
        expectNear(te[0], FACTOR * 64 * (14 + 77), "TE vacuum energy");
        expectNear(te[1], FACTOR * 64 * (4 * 14 + 77), "TE dielectric energy");
        expectNear(tm[0], te[0] + FACTOR * 16 * 9 * (1.6 - 1), "TM material face energy");
        expectNear(tm[1], te[1] + FACTOR * 16 * 9 * (1.6 - 4), "TM periodic face energy");
    }

    /// Distinguish the maximum of the two outer layers from the whole-domain maximum.
    template <class Field>
    void checkPeak(Field& electric) {
        plant(
            electric, KOKKOS_LAMBDA(int, int, int k, unsigned component) {
                if (component != 0)
                    return 0.0;
                return k == 0 ? -2.0 : k == 7 ? -3.0 : k == 4 ? 7.0 : 1.0;
            });
        expectNear(measureElectricPeak(electric, false), 7, "whole-domain electric peak");
        expectNear(measureElectricPeak(electric, true), 3, "second-face electric peak");
    }

    /// Every rank must observe failure when an owned field value is nonfinite.
    template <class Function>
    void expectNonfiniteFailure(Function measure) {
        bool failed = false;
        try {
            measure();
        } catch (const std::runtime_error&) {
            failed = true;
        }
        if (!failed)
            throw std::runtime_error("Nonfinite owned data was accepted");
    }

    /// Inject invalid data on one x slab to check collective rejection on the other rank.
    template <class Field>
    void checkNonfinite(Field& electric, Field& magnetic, const material::MaterialOwner& owner,
                        LevelContract contract) {
        plant(
            electric, KOKKOS_LAMBDA(int i, int, int, unsigned) {
                return i == 0 ? std::numeric_limits<double>::quiet_NaN() : 1.0;
            });
        std::array<double, 12> wavenumbers{};
        expectNonfiniteFailure([&] {
            measureElectricPeak(electric, true);
        });
        expectNonfiniteFailure([&] {
            measureSideEnergy(electric, magnetic, owner, true, 0.5);
        });
        expectNonfiniteFailure([&] {
            measureFourierMagnitude(electric, contract, Quantity::Electric, KX, wavenumbers,
                                    INTERFACE_Z, false);
        });
    }

    /// The x decomposition deliberately gives both ranks complex contributions from vacuum.
    void checkMeasurements() {
        auto fields = solver::makeBoxFieldContainer({8, 2, 8}, {1, 0.25, 1}, 0);
        initializeFieldStorage(fields);
        const material::MaterialOwner owner(fields.getMesh(), fields.getFL(), DielectricMap{});
        LevelContract contract;
        contract.spatialDiscretization_m = SpatialDiscretization::Staggered;
        checkSamplePositions(fields.getE(), contract, Quantity::Electric);
        checkSamplePositions(fields.getB(), contract, Quantity::Magnetic);
        checkGlobalCancellation(fields.getE(), contract);
        checkDirection(fields.getE(), contract, false);
        checkDirection(fields.getE(), contract, true);
        checkEnergy(fields.getE(), fields.getB(), owner);
        checkPeak(fields.getE());
        checkNonfinite(fields.getE(), fields.getB(), owner, contract);
    }
}  // namespace

int main(int argc, char** argv) {
    ippl::initialize(argc, argv);
    int failed = 0;
    try {
        checkMeasurements();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        failed = 1;
    }
    MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_MAX, *ippl::Comm);
    if (!failed && ippl::Comm->rank() == 0)
        std::cout << "E18/E19 prescribed-field measurements passed\n";
    ippl::finalize();
    return failed;
}
