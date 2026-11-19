/** @file E03NativeLines.h
 * @brief Samples of the potentials and fields along three grid lines, at their own positions.
 */
#pragma once

#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include "core/input/RunInput.h"
#include "core/output/RunRecord.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03Measurements.h"

namespace chdr::experiments::e03 {
    namespace nativeLineDetail {
        /// Stored per component: the three coordinates, the numerical and the continuum value.
        inline constexpr int VALUES_PER_COMPONENT = 5;

        /// The name of the quantity in the `quantity` column of `native-lines.csv`.
        inline const char* quantityName(Quantity quantity) {
            switch (quantity) {
                case Quantity::Electric:
                    return "E";
                case Quantity::Magnetic:
                    return "B";
                case Quantity::Phi:
                    return "phi";
                default:
                    return "A";
            }
        }

        /** @brief Store the position, the numerical and the continuum value of one cell.
         * @details Nothing is interpolated. Every component uses its position offset from the
         * LevelContract [ryu2016potentialfdtd, Appendix A, Eqs. (19)--(24)], and the continuum
         * value is taken at that position and at the comparison time of the quantity. The
         * four-potential field holds \f$\phi\f$ in component 0 and \f$\mathbf{A}\f$ in
         * components 1 to 3; \f$\mathbf{E}\f$ and \f$\mathbf{B}\f$ start at component 0.
         */
        template <class View, class Line, DipolePulseReference Reference>
        KOKKOS_INLINE_FUNCTION void storeSample(View view, Line line, int m,
                                                Kokkos::Array<int, 3> index,
                                                Kokkos::Array<double, 3> centre, double h,
                                                LevelContract contract, const Reference& ref,
                                                Quantity quantity, double time) {
            const unsigned count = quantity == Quantity::Phi ? 1 : 3;
            const unsigned first = quantity == Quantity::VectorPotential ? 1 : 0;
            for (unsigned c = 0; c < count; ++c) {
                const auto offset = contract.positionOffset(quantity, c);
                const int column  = VALUES_PER_COMPONENT * c;
                for (unsigned d = 0; d < 3; ++d)
                    line(m, column + d) = centre[d] + offset[d] * h;
                line(m, column + 3) = view(index[0], index[1], index[2])[c + first];
                line(m, column + 4) = e03MeasurementDetail::referenceValue(
                    ref, quantity, line(m, column), line(m, column + 1), line(m, column + 2), time,
                    c);
            }
        }

        /** @brief Sample one quantity in the cells of one grid line through the lower central cell.
         * @details The line runs along `axis` through the cell with the index \f$N/2-1\f$ in
         * the two other directions. A rank stores the cells that it owns. The entries of the
         * other cells stay zero, so that a sum over the ranks assembles the line.
         */
        template <class Field, class Line, DipolePulseReference Reference>
        void sampleLine(Field& field, Line line, unsigned axis, LevelContract contract,
                        Reference ref, Quantity quantity, double time) {
            const auto view  = field.getView();
            const auto local = field.getLayout().getLocalNDIndex();
            const auto first = local.first(), last = local.last();
            const auto origin = field.get_mesh().getOrigin();
            const double h    = field.get_mesh().getMeshSpacing()[0];
            const int n = line.extent_int(0), ghost = field.getNghost();
            using Exec = typename Field::view_type::execution_space;
            Kokkos::parallel_for(
                "E03 native line", Kokkos::RangePolicy<Exec>(0, n), KOKKOS_LAMBDA(int m) {
                    Kokkos::Array<int, 3> cell{n / 2 - 1, n / 2 - 1, n / 2 - 1};
                    cell[axis] = m;
                    Kokkos::Array<double, 3> centre{0.0, 0.0, 0.0};
                    for (unsigned d = 0; d < 3; ++d) {
                        if (cell[d] < first[d] || cell[d] > last[d])
                            return;
                        centre[d] = origin[d] + (cell[d] + 0.5) * h;
                        cell[d] += ghost - first[d];
                    }
                    storeSample(view, line, m, cell, centre, h, contract, ref, quantity, time);
                });
        }

        /// Assemble the line on rank zero; the result is empty on the other ranks.
        template <class Line>
        std::vector<double> gatherLine(Line line, MPI_Comm comm) {
            Kokkos::fence();
            const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, line);
            int rank;
            MPI_Comm_rank(comm, &rank);
            std::vector<double> samples(rank == 0 ? host.size() : 0);
            MPI_Reduce(host.data(), samples.data(), static_cast<int>(host.size()), MPI_DOUBLE,
                       MPI_SUM, 0, comm);
            return samples;
        }

        /// Append one row per sample, with the coordinates relative to the pulse centre.
        inline void writeLineRows(std::ofstream& file, Quantity quantity, unsigned axis,
                                  long long step, double time, const std::vector<double>& samples) {
            static constexpr std::array<const char*, 3> NAMES{"x", "y", "z"};
            const unsigned count  = quantity == Quantity::Phi ? 1 : 3;
            const std::size_t row = 3 * VALUES_PER_COMPONENT;
            for (std::size_t m = 0; m < samples.size() / row; ++m) {
                for (unsigned c = 0; c < count; ++c) {
                    const double* s = &samples[m * row + c * VALUES_PER_COMPONENT];
                    file << quantityName(quantity) << ',' << (count == 1 ? "-" : NAMES[c]) << ','
                         << NAMES[axis] << ',' << step << ',' << time << ',' << s[0] - E03BoxCentre
                         << ',' << s[1] - E03BoxCentre << ',' << s[2] - E03BoxCentre << ',' << s[3]
                         << ',' << s[4] << '\n';
                }
            }
        }
    }  // namespace nativeLineDetail

    /** @brief Optional samples of one run along the three grid lines through the lower central
     * cell, written to `native-lines.csv`.
     * @details The parameter `write_native_lines` enables the output. Each row holds one
     * component of one quantity at its own position, the numerical value and the continuum
     * value at that position and at the comparison time of the quantity. The samples enter no
     * measurement. They show the pointwise difference to the continuum without the
     * interpolation of `plane-z0.csv`.
     * @pre The run directory has been prepared collectively before construction.
     */
    class SphericalNativeLines {
    public:
        SphericalNativeLines(const input::ChdrYamlRunInput& request, MPI_Comm comm)
            : path_m(std::filesystem::path(request.outputDirectory_m) / "native-lines.csv")
            , enabled_m(input::booleanParameter(request, "write_native_lines", false)) {
            int rank;
            MPI_Comm_rank(comm, &rank);
            if (enabled_m && rank == 0)
                std::ofstream(path_m)
                    << "quantity,component,line_axis,step,time,x,y,z,numerical,continuum\n";
        }

        /// True if the run writes `native-lines.csv`.
        bool writes() const { return enabled_m; }

        /// Name the file in the run record when the run writes it.
        void record(output::RunRecordDescription& description) const {
            if (enabled_m)
                description.effective_m["native_lines/file"] = "native-lines.csv";
        }

        /** @brief Append the three lines of one quantity of this state.
         * @param field The four-potential field for \f$\mathbf{A}\f$ and \f$\phi\f$, otherwise
         * the reconstructed field of the quantity.
         * @param time Comparison time of the quantity in this state (comparisonTime).
         */
        template <class Field, DipolePulseReference Reference>
        void write(Field& field, LevelContract contract, const Reference& ref, Quantity quantity,
                   long long step, double time) const {
            const int n  = field.getLayout().getDomain()[0].length();
            using Device = typename Field::view_type::device_type;
            for (unsigned axis = 0; axis < 3; ++axis) {
                Kokkos::View<double**, Kokkos::LayoutRight, Device> line(
                    "E03 native line", n, 3 * nativeLineDetail::VALUES_PER_COMPONENT);
                nativeLineDetail::sampleLine(field, line, axis, contract, ref, quantity, time);
                append(nativeLineDetail::gatherLine(line, field.getLayout().comm), quantity, axis,
                       step, time);
            }
        }

    private:
        /// Rank zero holds the assembled samples; the other ranks have nothing to write.
        void append(const std::vector<double>& samples, Quantity quantity, unsigned axis,
                    long long step, double time) const {
            if (samples.empty())
                return;
            std::ofstream file;
            file.exceptions(std::ios::failbit | std::ios::badbit);
            file.open(path_m, std::ios::app);
            file << std::setprecision(std::numeric_limits<double>::max_digits10);
            nativeLineDetail::writeLineRows(file, quantity, axis, step, time, samples);
        }

        std::filesystem::path path_m;
        bool enabled_m;
    };
}  // namespace chdr::experiments::e03
