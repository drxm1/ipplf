/** @file VacuumOutput.h
 * @brief Synchronous owned-cell VTI/PVTI output for the vacuum demonstration.
 * @details Format sources (consulted 2026-10-02): VTK XML File Formats,
 * https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html ;
 * VTK Field Data as Time Meta-Data,
 * https://docs.vtk.org/en/latest/design_documents/IOXMLTimeInFieldData.html ;
 * ParaView User Guide, Handling Temporal File Series,
 * https://docs.paraview.org/en/latest/UsersGuide/dataIngestion.html#handling-temporal-file-series .
 * Logical x-fastest traversal follows VTK 9.5.2 vtkStructuredData::ComputeCellId
 * and GetLinearIndex (vtkStructuredData.h, lines 280-283 and 329-333):
 * https://github.com/Kitware/VTK/blob/v9.5.2/Common/DataModel/vtkStructuredData.h .
 * Caller-supplied timestamps are serialized unchanged; this writer derives no
 * physical array times from potential time or dt. No GPU validation is claimed.
 */
#ifndef CHDR_VACUUM_OUTPUT_H
#define CHDR_VACUUM_OUTPUT_H

#include <Kokkos_Core.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <mpi.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chdr::vacuum {

    /** @brief Zero-based global cell counts and global lower-point coordinates. */
    struct OutputGeometry {
        std::array<int, 3> cells_m;
        std::array<double, 3> origin_m;
        std::array<double, 3> spacing_m;
    };

    /** @brief Completed solver state: step is separate from output-frame numbering. */
    struct FrameRecord {
        int step_m;              ///< Completed solver step, not frame ordinal.
        double potentialTime_m;  ///< Dataset/series time supplied by the caller.
        double electricTime_m;   ///< Actual raw E sample time supplied by the caller.
        double magneticTime_m;   ///< Actual raw B sample time supplied by the caller.
        double dt_m;             ///< Actual positive solver timestep.
    };

    namespace output_detail {
        using Extent = std::array<int, 6>;  ///< Global point ranges, inclusive at both ends.
        static_assert(sizeof(Extent) == 6 * sizeof(int),
                      "MPI extent gather requires packed arrays");

        /** @brief Locale-independent numbers with double round-trip precision. */
        inline void configureStream(std::ostream& out) {
            out.imbue(std::locale::classic());
            out << std::setprecision(std::numeric_limits<double>::max_digits10);
        }

        /** @brief Format a small coordinate/extent array without assuming memory layout. */
        template <class T, std::size_t N>
        std::string arrayText(const std::array<T, N>& values) {
            std::ostringstream out;
            configureStream(out);
            for (std::size_t i = 0; i < N; ++i)
                out << (i ? " " : "") << values[i];
            return out.str();
        }

        /** @brief Global point extent for zero-based cell indices. */
        inline Extent wholeExtent(const OutputGeometry& geometry) {
            return {0, geometry.cells_m[0], 0, geometry.cells_m[1], 0, geometry.cells_m[2]};
        }

        /** @brief Controlled basename: step is nonnegative and rank is communicator-local. */
        inline std::string frameName(int step, int rank = -1) {
            const auto stem = "wave-" + std::to_string(step);
            return rank < 0 ? stem + ".pvti" : stem + "-rank-" + std::to_string(rank) + ".vti";
        }

        /** @brief One numerical metadata tuple; names/type strings are internal constants. */
        template <class T>
        void scalar(std::ostream& out, const char* name, const char* type, T value) {
            out << "<DataArray type=\"" << type << "\" Name=\"" << name
                << "\" NumberOfTuples=\"1\" format=\"ascii\">" << value << "</DataArray>\n";
        }

        /** @brief Explicit array times; normalization code 1 means dimensionless c=1, no SI scales.
         */
        inline void metadata(std::ostream& out, const FrameRecord& frame) {
            out << "<FieldData>\n";
            scalar(out, "TimeValue", "Float64", frame.potentialTime_m);
            scalar(out, "potential_time", "Float64", frame.potentialTime_m);
            scalar(out, "E_time", "Float64", frame.electricTime_m);
            scalar(out, "B_time", "Float64", frame.magneticTime_m);
            scalar(out, "dt", "Float64", frame.dt_m);
            scalar(out, "step", "Int64", frame.step_m);
            scalar(out, "normalization", "Int32", 1);
            out << "</FieldData>\n";
        }

        /** @brief Check geometry without assigning any physical unit conversion. */
        inline void validateGeometry(const OutputGeometry& geometry) {
            for (int d = 0; d < 3; ++d) {
                if (geometry.cells_m[d] < 1 || !std::isfinite(geometry.origin_m[d])
                    || !std::isfinite(geometry.spacing_m[d]) || geometry.spacing_m[d] <= 0)
                    throw std::runtime_error("Invalid output geometry");
            }
        }

        /** @brief Convert inclusive IPPL owned cell ranges to VTK point extents. */
        template <class Field>
        Extent fieldExtent(const Field& field) {
            const auto local = field.getLayout().getLocalNDIndex();
            Extent result{};
            for (int d = 0; d < 3; ++d) {
                result[2 * d]     = local[d].first();
                result[2 * d + 1] = local[d].last() + 1;
            }
            return result;
        }

        /** @brief Reject an incompatible global domain, local extent or halo allocation. */
        template <class Field>
        void validateField(const Field& field, const OutputGeometry& geometry,
                           const Extent& extent) {
            const auto domain = field.getLayout().getDomain();
            const int ghost   = field.getNghost();
            if (ghost < 0)
                throw std::runtime_error("Negative output halo width");
            for (int d = 0; d < 3; ++d) {
                const int first = extent[2 * d], end = extent[2 * d + 1];
                if (domain[d].first() != 0
                    || !std::cmp_equal(domain[d].length(), geometry.cells_m[d]) || first < 0
                    || end > geometry.cells_m[d] || end <= first
                    || field.getView().extent(d)
                           != static_cast<std::size_t>(end - first + 2 * ghost))
                    throw std::runtime_error("Output field geometry/allocation mismatch");
            }
        }

        /** @brief Validate matching E/B ownership before reading either host snapshot. */
        template <class EField, class BField>
        Extent checkedExtent(const EField& electric, const BField& magnetic,
                             const OutputGeometry& geometry) {
            const auto local = fieldExtent(electric);
            validateField(electric, geometry, local);
            validateField(magnetic, geometry, local);
            if (local != fieldExtent(magnetic) || electric.getNghost() != magnetic.getNghost())
                throw std::runtime_error("Output E/B ownership differs");
            return local;
        }

        /** @brief Write only owned vector tuples, x fastest; reject nonfinite components. */
        template <class View>
        void vectors(std::ostream& out, const char* name, const View& view, const Extent& e,
                     int ghost) {
            out << "<DataArray type=\"Float64\" Name=\"" << name
                << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (int k = 0; k < e[5] - e[4]; ++k)
                for (int j = 0; j < e[3] - e[2]; ++j)
                    for (int i = 0; i < e[1] - e[0]; ++i) {
                        const auto value = view(i + ghost, j + ghost, k + ghost);
                        for (int c = 0; c < 3; ++c) {
                            if (!std::isfinite(value[c]))
                                throw std::runtime_error("Nonfinite output field");
                            out << value[c] << (c == 2 ? '\n' : ' ');
                        }
                    }
            out << "</DataArray>\n";
        }

        /** @brief Flush/close a temporary file before atomically publishing its final name. */
        template <class Write>
        void atomicFile(const std::filesystem::path& path, Write write) {
            const auto temporary = path.string() + ".tmp";
            try {
                std::ofstream out;
                out.exceptions(std::ios::failbit | std::ios::badbit);
                configureStream(out);
                out.open(temporary, std::ios::out | std::ios::trunc);
                write(out);
                out.flush();
                out.close();
                std::filesystem::rename(temporary, path);
            } catch (const std::exception& error) {
                throw std::runtime_error("Output file " + path.string() + ": " + error.what());
            }
        }

        /** @brief Serialize one rank with global origin and its global-index point extent. */
        template <class EView, class BView>
        void piece(std::ostream& out, const OutputGeometry& geometry, const Extent& e,
                   const EView& electric, const BView& magnetic, int ghost,
                   const FrameRecord& frame) {
            out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"ImageData\" version=\"1.0\" "
                   "byte_order=\"LittleEndian\">\n"
                << "<ImageData WholeExtent=\"" << arrayText(e) << "\" Origin=\""
                << arrayText(geometry.origin_m) << "\" Spacing=\"" << arrayText(geometry.spacing_m)
                << "\">\n";
            metadata(out, frame);
            out << "<Piece Extent=\"" << arrayText(e) << "\"><CellData Vectors=\"E_raw\">\n";
            vectors(out, "E_raw", electric, e, ghost);
            vectors(out, "B_raw", magnetic, e, ghost);
            out << "</CellData><PointData/></Piece>\n</ImageData>\n</VTKFile>\n";
        }

        /** @brief Publish a parallel wrapper, including the serial one-piece case. */
        inline void wrapper(std::ostream& out, const OutputGeometry& geometry,
                            const std::vector<Extent>& extents, const FrameRecord& frame) {
            out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PImageData\" version=\"1.0\" "
                   "byte_order=\"LittleEndian\">\n"
                << "<PImageData WholeExtent=\"" << arrayText(wholeExtent(geometry))
                << "\" GhostLevel=\"0\" Origin=\"" << arrayText(geometry.origin_m)
                << "\" Spacing=\"" << arrayText(geometry.spacing_m) << "\">\n";
            metadata(out, frame);
            out << "<PPointData/><PCellData Vectors=\"E_raw\">\n"
                << "<PDataArray type=\"Float64\" Name=\"E_raw\" NumberOfComponents=\"3\"/>\n"
                << "<PDataArray type=\"Float64\" Name=\"B_raw\" "
                   "NumberOfComponents=\"3\"/>\n</PCellData>\n";
            for (std::size_t rank = 0; rank < extents.size(); ++rank)
                out << "<Piece Extent=\"" << arrayText(extents[rank]) << "\" Source=\""
                    << frameName(frame.step_m, static_cast<int>(rank)) << "\"/>\n";
            out << "</PImageData>\n</VTKFile>\n";
        }

        /** @brief Explicit potential-time list; no inference from filenames/frame numbers. */
        inline void series(std::ostream& out, const std::vector<FrameRecord>& frames) {
            out << "{\"file-series-version\":\"1.0\",\"files\":[\n";
            for (std::size_t i = 0; i < frames.size(); ++i) {
                out << (i ? ",\n" : "") << "{\"name\":\"" << frameName(frames[i].step_m)
                    << "\",\"time\":" << frames[i].potentialTime_m << "}";
            }
            out << "\n]}\n";
        }
    }  // namespace output_detail

    /** @brief Small collective output owner; constructor and calls require every rank.
     * @details Directory must be new or empty on a shared filesystem. Call with
     * identical geometry, explicit frame metadata and cadence on every rank. E/B must share
     * owned extents and halo width. The caller must pass the fields' actual mesh
     * origin/spacing: this generic adapter checks index domains, not mesh objects.
     * Host mirrors exist only during output; fields
     * stay alive and unchanged until return. Normal I/O exceptions are reduced
     * before subsequent collectives; any thrown exception is terminal for this
     * writer. Fatal MPI/process/allocation failure recovery is outside this demo.
     * Each piece closes before collective success, then rank zero publishes PVTI.
     * Only successful PVTI frames enter writeSeries(); temporary/orphan pieces
     * after a failed frame are never advertised. No filesystem cleanup is hidden.
     */
    class VacuumOutput {
    public:
        /** @brief Collectively prepare a fresh output directory; borrow the communicator. */
        VacuumOutput(std::filesystem::path directory, OutputGeometry geometry,
                     MPI_Comm communicator)
            : directory_m(std::move(directory))
            , geometry_m(geometry)
            , communicator_m(communicator) {
            MPI_Comm_rank(communicator_m, &rank_m);
            MPI_Comm_size(communicator_m, &size_m);
            collectively([&] {
                output_detail::validateGeometry(geometry_m);
                extents_m.resize(size_m);
                if (rank_m == 0)
                    prepareDirectory();
            });
        }

        /** @brief Write a completed raw E/B state with caller-supplied measured timestamps. */
        template <class EField, class BField>
        void writeFrame(const EField& electric, const BField& magnetic, const FrameRecord& frame) {
            output_detail::Extent local{};
            collectively([&] {
                validateFrame(frame);
                local = writeLocalPiece(electric, magnetic, frame);
            });
            MPI_Allgather(local.data(), 6, MPI_INT, extents_m.data(), 6, MPI_INT, communicator_m);
            collectively([&] {
                if (rank_m == 0)
                    publishWrapper(frame);
            });
            collectively([&] {
                frames_m.push_back(frame);
            });
        }

        /** @brief Atomically publish the accepted frames; collective, safe after each frame. */
        void writeSeries() const {
            collectively([&] {
                if (rank_m == 0)
                    output_detail::atomicFile(directory_m / "wave.pvti.series",
                                              [&](std::ostream& out) {
                                                  output_detail::series(out, frames_m);
                                              });
            });
        }

    private:
        /** @brief Turn caught rank-local errors into the same failed call on every rank. */
        template <class Action>
        void collectively(Action action) const {
            std::string message;
            try {
                action();
            } catch (const std::exception& error) {
                message = error.what();
            }
            const int failed = !message.empty();
            int anyFailed    = 0;
            MPI_Allreduce(&failed, &anyFailed, 1, MPI_INT, MPI_MAX, communicator_m);
            if (!message.empty())
                std::cerr << "Output rank " << rank_m << ": " << message << '\n';
            if (anyFailed)
                throw std::runtime_error("Collective vacuum output failed in "
                                         + directory_m.string());
        }

        /** @brief Avoid accidentally presenting stale output from an earlier run. */
        void prepareDirectory() const {
            if (std::filesystem::exists(directory_m) && !std::filesystem::is_empty(directory_m))
                throw std::runtime_error("Output directory is not empty: " + directory_m.string());
            std::filesystem::create_directories(directory_m);
        }

        /** @brief Require finite times and advancing state labels; impose no E/B offset formula. */
        void validateFrame(const FrameRecord& frame) const {
            if (frame.step_m < 0 || !std::isfinite(frame.potentialTime_m)
                || !std::isfinite(frame.dt_m) || frame.dt_m <= 0
                || !std::isfinite(frame.electricTime_m) || !std::isfinite(frame.magneticTime_m))
                throw std::runtime_error("Invalid output step/time/dt");
            if (!frames_m.empty()
                && (frame.step_m <= frames_m.back().step_m
                    || frame.potentialTime_m <= frames_m.back().potentialTime_m))
                throw std::runtime_error("Output frames must advance in step and time");
        }

        /** @brief Mirror E/B by logical index and synchronously finish this rank's file. */
        template <class EField, class BField>
        output_detail::Extent writeLocalPiece(const EField& electric, const BField& magnetic,
                                              const FrameRecord& frame) const {
            const auto local = output_detail::checkedExtent(electric, magnetic, geometry_m);
            Kokkos::fence();
            const auto e =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), electric.getView());
            const auto b =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), magnetic.getView());
            output_detail::atomicFile(directory_m / output_detail::frameName(frame.step_m, rank_m),
                                      [&](std::ostream& out) {
                                          output_detail::piece(out, geometry_m, local, e, b,
                                                               electric.getNghost(), frame);
                                      });
            return local;
        }

        /** @brief Rank-zero wrapper publication follows collective piece success. */
        void publishWrapper(const FrameRecord& frame) const {
            output_detail::atomicFile(directory_m / output_detail::frameName(frame.step_m),
                                      [&](std::ostream& out) {
                                          output_detail::wrapper(out, geometry_m, extents_m, frame);
                                      });
        }

        std::filesystem::path directory_m;
        OutputGeometry geometry_m;
        MPI_Comm communicator_m;
        int rank_m = 0;
        int size_m = 1;
        std::vector<output_detail::Extent> extents_m;
        std::vector<FrameRecord> frames_m;
    };
}  // namespace chdr::vacuum
#endif
