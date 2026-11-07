/** @file FieldOutput.h
 * @brief Synchronous owned-cell VTI/PVTI field and material output for ChDR demonstrations.
 * @details Format sources (consulted 2026-10-02): VTK XML File Formats,
 * https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html ;
 * VTK Field Data as Time Meta-Data,
 * https://docs.vtk.org/en/latest/design_documents/IOXMLTimeInFieldData.html ;
 * ParaView User Guide, Handling Temporal File Series,
 * https://docs.paraview.org/en/latest/UsersGuide/dataIngestion.html#handling-temporal-file-series .
 * Logical x-fastest traversal follows VTK 9.5.2 vtkStructuredData::ComputeCellId
 * and GetLinearIndex (vtkStructuredData.h, lines 280-283 and 329-333):
 * https://github.com/Kitware/VTK/blob/v9.5.2/Common/DataModel/vtkStructuredData.h .
 * LevelContract supplies each quantity's time and component positions. The caller
 * selects the arrays, series clock and output units without redefining that contract.
 * No GPU validation is claimed.
 */
#ifndef CHDR_FIELD_OUTPUT_H
#define CHDR_FIELD_OUTPUT_H

#include <Kokkos_Core.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <mpi.h>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/base/FieldContainer.h"
#include "core/base/LevelContract.h"

namespace chdr::output {

    /** @brief Zero-based global cell counts and global lower-point coordinates. */
    struct OutputGeometry {
        std::array<int, 3> cells_m;
        std::array<double, 3> origin_m;
        std::array<double, 3> spacing_m;
    };

    /** @brief Completed solver state: step is separate from output-frame numbering. */
    struct FrameRecord {
        int step_m;           ///< Completed solver step, not frame ordinal.
        double seriesTime_m;  ///< Dataset/series clock; array times follow LevelContract.
        double dt_m;          ///< Actual positive solver timestep, in solver units.
    };

    /** @brief Bind a native vector field to its identity and output conversion.
     * @details The reference is borrowed only during writeFrame(). Quantity identifies
     * the sampling contract; it does not reconstruct or modify the supplied values.
     * Scalar quantities require a different storage interface and are rejected here.
     */
    struct VectorFieldSnapshot {
        std::string name_m;
        const EMFieldsType& field_m;
        Quantity quantity_m;
        double valueScale_m;
    };

    /** @brief Multiplicative conversions applied only to serialized output.
     * @details Factors convert the supplied solver quantities to the selected output units.
     * Geometry and times use these factors; each vector description supplies its value factor.
     * Component offsets remain dimensionless fractions of mesh spacing. No solver array
     * is changed. si_m labels SI output; default factors preserve normalized output.
     */
    struct OutputScales {
        double lengthFactor_m = 1, timeFactor_m = 1;
        bool si_m = false;

        /// Reject unusable scales before writing any output file.
        void validate() const {
            for (double scale : {lengthFactor_m, timeFactor_m})
                if (!std::isfinite(scale) || scale <= 0)
                    throw std::invalid_argument("Output scale must be finite and positive");
        }
    };

    /// Obtain output factors from a unit convention without making the writer own input parsing.
    template <class UnitSystem>
    OutputScales fieldOutputScales(const UnitSystem& units) {
        return {units.lengthFactor(), units.timeFactor(), units.isSI()};
    }

    namespace output_detail {
        using Extent = std::array<int, 6>;  ///< Global point ranges, inclusive at both ends.
        static_assert(sizeof(Extent) == 6 * sizeof(int),
                      "MPI extent gather requires packed arrays");

        /// Evaluated metadata owns no solver arrays; offsets are in mesh-spacing units.
        struct VectorMetadata {
            std::string name_m;
            Quantity quantity_m;
            double sampleTime_m, valueScale_m;
            /** @brief Flattened 3-by-3 table of native component positions.
             * @details Entry [3*c+d] is the displacement along coordinate d of field
             * component c; both indices follow x, y, z order. Displacements are relative
             * to the logical cell centre, in units of mesh spacing along d.
             * LevelContract::positionOffset supplies each row. This ordering writes
             * three VTK tuples: the (x,y,z) offset of field-x, then field-y, then field-z.
             */
            std::array<double, 9> componentOffsets_m;
        };

        /// One resolved frame, with converted times and all per-vector sampling metadata.
        struct FrameMetadata {
            FrameRecord frame_m;
            LevelContract contract_m;
            bool phiIsZero_m;
            std::vector<VectorMetadata> vectors_m;
        };

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
        inline std::string frameName(const std::string& basename, int step, int rank = -1) {
            const auto stem = basename + "-" + std::to_string(step);
            return rank < 0 ? stem + ".pvti" : stem + "-rank-" + std::to_string(rank) + ".vti";
        }

        /** @brief One numerical metadata tuple; names/type strings are internal constants. */
        template <class T>
        void scalar(std::ostream& out, const std::string& name, const char* type, T value) {
            out << "<DataArray type=\"" << type << "\" Name=\"" << name
                << "\" NumberOfTuples=\"1\" format=\"ascii\">" << value << "</DataArray>\n";
        }

        /** @brief Write a vector's contract-derived time, offsets and conversion factor.
         * @details Quantity codes follow Quantity: E=0, B=1, raw B=2, phi=3, A=4,
         * rho=5, J=6. Only vector quantities are accepted by this writer.
         */
        inline void vectorMetadata(std::ostream& out, const VectorMetadata& vector) {
            scalar(out, vector.name_m + "_time", "Float64", vector.sampleTime_m);
            scalar(out, vector.name_m + "_scale_from_solver", "Float64", vector.valueScale_m);
            scalar(out, vector.name_m + "_quantity", "Int32", static_cast<int>(vector.quantity_m));
            out << "<DataArray type=\"Float64\" Name=\"" << vector.name_m
                << "_component_offsets\" NumberOfComponents=\"3\" "
                   "NumberOfTuples=\"3\" format=\"ascii\">"
                << arrayText(vector.componentOffsets_m) << "</DataArray>\n";
        }

        /** @brief Record the supplied contract without identifying it through an array name.
         * @details Reconstruction codes are R1=0, R2=1, R2Prime=2; spatial codes are
         * collocated=0 and staggered=1. phi_is_zero_declared records the caller's explicit
         * declaration.
         */
        inline void contractMetadata(std::ostream& out, const FrameMetadata& metadata) {
            scalar(out, "reconstruction_convention", "Int32",
                   static_cast<int>(metadata.contract_m.ebReconstructionConvention_m));
            scalar(out, "spatial_discretization", "Int32",
                   static_cast<int>(metadata.contract_m.spatialDiscretization_m));
            scalar(out, "phi_is_zero_declared", "Int32", metadata.phiIsZero_m ? 1 : 0);
        }

        /// Common and per-array metadata; normalization is 1 for solver units, 0 for SI.
        inline void metadata(std::ostream& out, const FrameMetadata& metadata,
                             OutputScales scales) {
            const auto& frame = metadata.frame_m;
            out << "<FieldData>\n";
            scalar(out, "TimeValue", "Float64", frame.seriesTime_m);
            scalar(out, "series_time", "Float64", frame.seriesTime_m);
            scalar(out, "dt", "Float64", frame.dt_m);
            scalar(out, "step", "Int64", frame.step_m);
            scalar(out, "normalization", "Int32", scales.si_m ? 0 : 1);
            scalar(out, "length_scale_from_solver", "Float64", scales.lengthFactor_m);
            scalar(out, "time_scale_from_solver", "Float64", scales.timeFactor_m);
            contractMetadata(out, metadata);
            for (const auto& vector : metadata.vectors_m)
                vectorMetadata(out, vector);
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

        /// Convert geometry while retaining the global cell counts and ownership indices.
        inline OutputGeometry scaledGeometry(OutputGeometry geometry, OutputScales scales) {
            for (unsigned d = 0; d < 3; ++d) {
                geometry.origin_m[d] *= scales.lengthFactor_m;
                geometry.spacing_m[d] *= scales.lengthFactor_m;
            }
            return geometry;
        }

        /// Convert the series clock and timestep without changing the completed-solve index.
        inline FrameRecord scaledFrame(FrameRecord frame, OutputScales scales) {
            frame.seriesTime_m *= scales.timeFactor_m;
            frame.dt_m *= scales.timeFactor_m;
            return frame;
        }

        /// Restrict caller-supplied names to characters safe in XML, JSON and path basenames.
        inline void validateName(const std::string& name) {
            if (name.empty()
                || name.find_first_not_of(
                       "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-")
                       != std::string::npos)
                throw std::invalid_argument("Invalid field-output name: " + name);
        }

        /// Require distinct array names and metadata keys before any rank publishes a piece.
        inline void validateNames(std::span<const VectorFieldSnapshot> fields) {
            if (fields.empty())
                throw std::invalid_argument("A field frame must contain at least one vector");
            std::set<std::string> names;
            std::set<std::string> keys{"TimeValue",
                                       "series_time",
                                       "dt",
                                       "step",
                                       "normalization",
                                       "length_scale_from_solver",
                                       "time_scale_from_solver",
                                       "reconstruction_convention",
                                       "spatial_discretization",
                                       "phi_is_zero_declared"};
            for (const auto& field : fields) {
                validateName(field.name_m);
                if (!names.insert(field.name_m).second)
                    throw std::invalid_argument("Duplicate output vector: " + field.name_m);
                for (const char* suffix :
                     {"_time", "_component_offsets", "_scale_from_solver", "_quantity"})
                    if (!keys.insert(field.name_m + suffix).second)
                        throw std::invalid_argument("Conflicting vector metadata: " + field.name_m);
            }
        }

        /// The concrete field stores three Cartesian components, not a scalar-potential slot.
        inline void validateVectorQuantity(Quantity quantity) {
            switch (quantity) {
                case Quantity::Electric:
                case Quantity::Magnetic:
                case Quantity::RawMagnetic:
                case Quantity::VectorPotential:
                case Quantity::Current:
                    return;
                default:
                    throw std::invalid_argument("Field output requires a vector Quantity");
            }
        }

        /** @brief Resolve sampling metadata through the existing LevelContract methods.
         * @details The contract owns the time and placement rules; no reconstruction is
         * performed here. In particular, labelTime preserves the R1 electric restriction.
         */
        inline VectorMetadata resolveVector(const VectorFieldSnapshot& field,
                                            const LevelContract& contract, const FrameRecord& frame,
                                            double timeScale, bool phiIsZero) {
            validateVectorQuantity(field.quantity_m);
            if (!std::isfinite(field.valueScale_m) || field.valueScale_m <= 0)
                throw std::invalid_argument("Vector output scale must be finite and positive");
            VectorMetadata result{
                field.name_m,
                field.quantity_m,
                contract.labelTime(field.quantity_m, frame.step_m, frame.dt_m, phiIsZero)
                    * timeScale,
                field.valueScale_m,
                {}};
            if (!std::isfinite(result.sampleTime_m))
                throw std::invalid_argument("Nonfinite vector sample time");
            for (unsigned c = 0; c < 3; ++c) {
                const auto offset = contract.positionOffset(field.quantity_m, c);
                for (unsigned d = 0; d < 3; ++d)
                    result.componentOffsets_m[3 * c + d] = offset[d];
            }
            return result;
        }

        /// Resolve frame metadata without retaining references to the supplied field storage.
        inline FrameMetadata resolveFrame(const FrameRecord& frame, const LevelContract& contract,
                                          std::span<const VectorFieldSnapshot> fields,
                                          OutputScales scales, bool phiIsZero) {
            validateNames(fields);
            FrameMetadata result{scaledFrame(frame, scales), contract, phiIsZero, {}};
            result.vectors_m.reserve(fields.size());
            for (const auto& field : fields)
                result.vectors_m.push_back(
                    resolveVector(field, contract, frame, scales.timeFactor_m, phiIsZero));
            return result;
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

        /// Validate matching ownership before reading any host snapshot; fields is nonempty.
        inline Extent checkedExtent(std::span<const VectorFieldSnapshot> fields,
                                    const OutputGeometry& geometry) {
            const auto local = fieldExtent(fields.front().field_m);
            for (const auto& snapshot : fields) {
                validateField(snapshot.field_m, geometry, local);
                if (local != fieldExtent(snapshot.field_m))
                    throw std::runtime_error("Output vector ownership differs");
            }
            return local;
        }

        /** @brief Write only owned vector tuples, x fastest; reject nonfinite components. */
        template <class View>
        void vectors(std::ostream& out, const char* name, const View& view, const Extent& e,
                     int ghost, double scale = 1) {
            out << "<DataArray type=\"Float64\" Name=\"" << name
                << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (int k = 0; k < e[5] - e[4]; ++k)
                for (int j = 0; j < e[3] - e[2]; ++j)
                    for (int i = 0; i < e[1] - e[0]; ++i) {
                        const auto value = view(i + ghost, j + ghost, k + ghost);
                        for (int c = 0; c < 3; ++c) {
                            const double converted = value[c] * scale;
                            if (!std::isfinite(converted))
                                throw std::runtime_error("Nonfinite output field");
                            out << converted << (c == 2 ? '\n' : ' ');
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
                out.open(temporary, std::ios::out | std::ios::trunc | std::ios::binary);
                write(out);
                out.flush();
                out.close();
                std::filesystem::rename(temporary, path);
            } catch (const std::exception& error) {
                throw std::runtime_error("Output file " + path.string() + ": " + error.what());
            }
        }

        /** @brief Append one Float64 vector block with its UInt64 byte count (VTK XML format).
         * Logical x-fastest traversal excludes halos and preserves the ASCII tuple order.
         */
        template <class View>
        void binaryVectors(std::ostream& out, const View& view, const Extent& e, int ghost,
                           double scale = 1) {
            std::vector<double> values;
            values.reserve(3ULL * (e[1] - e[0]) * (e[3] - e[2]) * (e[5] - e[4]));
            for (int k = 0; k < e[5] - e[4]; ++k)
                for (int j = 0; j < e[3] - e[2]; ++j)
                    for (int i = 0; i < e[1] - e[0]; ++i)
                        for (int c = 0; c < 3; ++c) {
                            const double value = scale * view(i + ghost, j + ghost, k + ghost)[c];
                            if (!std::isfinite(value))
                                throw std::runtime_error("Nonfinite output field");
                            values.push_back(value);
                        }
            const std::uint64_t bytes = values.size() * sizeof(double);
            out.write(reinterpret_cast<const char*>(&bytes), sizeof(bytes));
            out.write(reinterpret_cast<const char*>(values.data()), bytes);
        }

        /// Begin a serial image or a parallel wrapper with explicit Cartesian geometry.
        inline void imageHeader(std::ostream& out, const OutputGeometry& geometry,
                                const Extent& extent, bool parallel) {
            const char* type = parallel ? "PImageData" : "ImageData";
            out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"" << type
                << "\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
                << '<' << type << " WholeExtent=\"" << arrayText(extent) << '"';
            if (parallel)
                out << " GhostLevel=\"0\"";
            out << " Origin=\"" << arrayText(geometry.origin_m) << "\" Spacing=\""
                << arrayText(geometry.spacing_m) << "\">\n";
        }

        /// Declare all vectors in the same order used by the appended binary blocks.
        inline void vectorDeclarations(std::ostream& out, const std::vector<VectorMetadata>& fields,
                                       bool parallel, std::uint64_t stride = 0) {
            for (std::size_t i = 0; i < fields.size(); ++i) {
                out << '<' << (parallel ? "PDataArray" : "DataArray") << " type=\"Float64\" Name=\""
                    << fields[i].name_m << "\" NumberOfComponents=\"3\"";
                if (!parallel)
                    out << " format=\"appended\" offset=\"" << i * stride << '"';
                out << "/>\n";
            }
        }

        /// Mirror one field at a time; copies finish before host serialization reads the values.
        inline void writeVectors(std::ostream& out, const Extent& extent,
                                 std::span<const VectorFieldSnapshot> fields, bool binary) {
            for (const auto& snapshot : fields) {
                const auto& field = snapshot.field_m;
                const auto host =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field.getView());
                if (binary)
                    binaryVectors(out, host, extent, field.getNghost(), snapshot.valueScale_m);
                else
                    vectors(out, snapshot.name_m.c_str(), host, extent, field.getNghost(),
                            snapshot.valueScale_m);
            }
        }

        /// Serialize one rank's owned vectors with a common mesh and per-array metadata.
        inline void piece(std::ostream& out, const OutputGeometry& geometry, const Extent& extent,
                          std::span<const VectorFieldSnapshot> fields, const FrameMetadata& frame,
                          bool binary, OutputScales scales) {
            imageHeader(out, geometry, extent, false);
            metadata(out, frame, scales);
            out << "<Piece Extent=\"" << arrayText(extent) << "\"><CellData Vectors=\""
                << fields.front().name_m << "\">\n";
            if (binary) {
                const std::uint64_t stride = 8
                                             + 3ULL * (extent[1] - extent[0])
                                                   * (extent[3] - extent[2])
                                                   * (extent[5] - extent[4]) * sizeof(double);
                vectorDeclarations(out, frame.vectors_m, false, stride);
            } else
                writeVectors(out, extent, fields, false);
            out << "</CellData><PointData/></Piece>\n</ImageData>\n";
            if (binary) {
                out << "<AppendedData encoding=\"raw\">_";
                writeVectors(out, extent, fields, true);
                out << "</AppendedData>\n";
            }
            out << "</VTKFile>\n";
        }

        /// Publish a parallel wrapper, including the serial one-piece case.
        inline void wrapper(std::ostream& out, const OutputGeometry& geometry,
                            const std::vector<Extent>& extents, const FrameMetadata& frame,
                            const std::string& basename, OutputScales scales) {
            imageHeader(out, geometry, wholeExtent(geometry), true);
            metadata(out, frame, scales);
            out << "<PPointData/><PCellData Vectors=\"" << frame.vectors_m.front().name_m
                << "\">\n";
            vectorDeclarations(out, frame.vectors_m, true);
            out << "</PCellData>\n";
            for (std::size_t rank = 0; rank < extents.size(); ++rank)
                out << "<Piece Extent=\"" << arrayText(extents[rank]) << "\" Source=\""
                    << frameName(basename, frame.frame_m.step_m, static_cast<int>(rank))
                    << "\"/>\n";
            out << "</PImageData>\n</VTKFile>\n";
        }

        /// Explicit series clocks; no inference from filenames, array times or frame ordinals.
        inline void series(std::ostream& out, const std::vector<FrameRecord>& frames,
                           const std::string& basename) {
            out << "{\"file-series-version\":\"1.0\",\"files\":[\n";
            for (std::size_t i = 0; i < frames.size(); ++i) {
                out << (i ? ",\n" : "") << "{\"name\":\"" << frameName(basename, frames[i].step_m)
                    << "\",\"time\":" << frames[i].seriesTime_m << "}";
            }
            out << "\n]}\n";
        }
        /// Write owned scalar permittivity in the same logical order as vector tuples.
        template <class View>
        void materialValues(std::ostream& out, const View& host, const Extent& extent, int ghost) {
            for (int k = 0; k < extent[5] - extent[4]; ++k)
                for (int j = 0; j < extent[3] - extent[2]; ++j)
                    for (int i = 0; i < extent[1] - extent[0]; ++i) {
                        const double value = host(i + ghost, j + ghost, k + ghost);
                        if (!std::isfinite(value))
                            throw std::runtime_error("Nonfinite material output");
                        out << value << '\n';
                    }
        }

        /// Material output is optional and does not configure vector sampling metadata.
        template <class View>
        void materialPiece(std::ostream& out, const OutputGeometry& geometry, const Extent& extent,
                           const View& host, int ghost) {
            imageHeader(out, geometry, extent, false);
            out << "<Piece Extent=\"" << arrayText(extent)
                << "\"><CellData><DataArray type=\"Float64\" Name=\"epsilon_r\" "
                   "format=\"ascii\">\n";
            materialValues(out, host, extent, ghost);
            out << "</DataArray></CellData><PointData/></Piece></ImageData></VTKFile>\n";
        }

        /// Declare the independently written immutable cell-permittivity pieces.
        inline void materialWrapper(std::ostream& out, const OutputGeometry& geometry,
                                    const std::vector<Extent>& extents) {
            imageHeader(out, geometry, wholeExtent(geometry), true);
            out << "<PPointData/><PCellData><PDataArray type=\"Float64\" "
                   "Name=\"epsilon_r\"/></PCellData>\n";
            for (std::size_t rank = 0; rank < extents.size(); ++rank)
                out << "<Piece Extent=\"" << arrayText(extents[rank]) << "\" Source=\"material-"
                    << rank << ".vti\"/>\n";
            out << "</PImageData>\n</VTKFile>\n";
        }
    }  // namespace output_detail

    /** @brief Small collective output owner; constructor and calls require every rank.
     * @details Directory must be new or empty on a shared filesystem. Call with
     * identical geometry, contract, vector descriptions and cadence on every rank.
     * The supplied fields must share owned extents. The caller must pass their actual mesh
     * origin/spacing: this generic adapter checks index domains, not mesh objects.
     * Host mirrors exist only during output; fields
     * stay alive and unchanged until return. Normal I/O exceptions are reduced
     * before subsequent collectives; any thrown exception is terminal for this
     * writer. Fatal MPI/process/allocation failure recovery is outside this demo.
     * Each piece closes before collective success, then rank zero publishes PVTI.
     * Only successful PVTI frames enter writeSeries(); temporary/orphan pieces
     * after a failed frame are never advertised. No filesystem cleanup is hidden.
     */
    class FieldOutput {
    public:
        /** @brief Collectively prepare a fresh output directory; borrow the communicator. */
        FieldOutput(std::filesystem::path directory, OutputGeometry geometry, MPI_Comm communicator,
                    std::string basename, bool binary = false, OutputScales scales = {})
            : directory_m(std::move(directory))
            , geometry_m(geometry)
            , communicator_m(communicator)
            , basename_m(std::move(basename))
            , binary_m(binary)
            , scales_m(scales) {
            MPI_Comm_rank(communicator_m, &rank_m);
            MPI_Comm_size(communicator_m, &size_m);
            collectively([&] {
                prepareOutput();
            });
        }

        /** @brief Write selected native vectors from one completed solver state.
         * @details LevelContract supplies all array times and component positions. The
         * optional phiIsZero declaration is passed through to its R1 electric time check;
         * this writer never infers it from a field, name or experiment. Fields are not
         * modified or retained. Series time is independent of individual sample times.
         * @pre Every rank supplies the same names, quantities, factors, contract and frame.
         */
        void writeFrame(const FrameRecord& frame, const LevelContract& contract,
                        std::span<const VectorFieldSnapshot> fields, bool phiIsZero = false) {
            output_detail::Extent local{};
            output_detail::FrameMetadata resolved{};
            collectively([&] {
                resolved =
                    output_detail::resolveFrame(frame, contract, fields, scales_m, phiIsZero);
                validateFrame(resolved.frame_m);
                local = writeLocalPiece(fields, resolved);
            });
            MPI_Allgather(local.data(), 6, MPI_INT, extents_m.data(), 6, MPI_INT, communicator_m);
            collectively([&] {
                if (rank_m == 0)
                    publishWrapper(resolved);
            });
            collectively([&] {
                frames_m.push_back(resolved.frame_m);
            });
        }

        /** @brief Atomically publish the accepted frames; collective, safe after each frame. */
        void writeSeries() const {
            collectively([&] {
                if (rank_m == 0)
                    output_detail::atomicFile(directory_m / (basename_m + ".pvti.series"),
                                              [&](std::ostream& out) {
                                                  output_detail::series(out, frames_m, basename_m);
                                              });
            });
        }

        /** @brief Write immutable cell permittivity once, independently of field frames.
         * @details geometry supplies ownership and halo width; epsilon is its scalar view.
         * All ranks call once before the first frame. Vector component positions are
         * recorded by writeFrame even when no material file is requested.
         */
        template <class Field, class View>
        void writeMaterial(const Field& geometry, View epsilon) {
            const auto local = output_detail::fieldExtent(geometry);
            collectively([&] {
                if (materialWritten_m || !frames_m.empty())
                    throw std::logic_error("Write the material once before the first frame");
                output_detail::validateField(geometry, geometry_m, local);
                writeLocalMaterial(geometry, epsilon, local);
            });
            MPI_Allgather(local.data(), 6, MPI_INT, extents_m.data(), 6, MPI_INT, communicator_m);
            collectively([&] {
                publishMaterial();
            });
            materialWritten_m = true;
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
                throw std::runtime_error("Collective field output failed in "
                                         + directory_m.string());
        }

        /// Resolve output geometry and prepare the directory within the collective error path.
        void prepareOutput() {
            scales_m.validate();
            output_detail::validateName(basename_m);
            geometry_m = output_detail::scaledGeometry(geometry_m, scales_m);
            output_detail::validateGeometry(geometry_m);
            if (binary_m && std::endian::native != std::endian::little)
                throw std::runtime_error("Binary VTK output requires a little-endian host");
            extents_m.resize(size_m);
            if (rank_m == 0)
                prepareDirectory();
        }

        /** @brief Avoid accidentally presenting stale output from an earlier run. */
        void prepareDirectory() const {
            if (std::filesystem::exists(directory_m) && !std::filesystem::is_empty(directory_m))
                throw std::runtime_error("Output directory is not empty: " + directory_m.string());
            std::filesystem::create_directories(directory_m);
        }

        /// Require a finite series clock and positive timestep without assigning array times.
        void validateFrame(const FrameRecord& frame) const {
            if (frame.step_m < 0 || !std::isfinite(frame.seriesTime_m) || !std::isfinite(frame.dt_m)
                || frame.dt_m <= 0)
                throw std::runtime_error("Invalid output step/time/dt");
            if (!frames_m.empty()
                && (frame.step_m <= frames_m.back().step_m
                    || frame.seriesTime_m <= frames_m.back().seriesTime_m))
                throw std::runtime_error("Output frames must advance in step and time");
        }

        /// Mirror selected vectors by logical index and synchronously finish this rank's file.
        output_detail::Extent writeLocalPiece(std::span<const VectorFieldSnapshot> fields,
                                              const output_detail::FrameMetadata& frame) const {
            const auto local = output_detail::checkedExtent(fields, geometry_m);
            Kokkos::fence();
            output_detail::atomicFile(
                directory_m / output_detail::frameName(basename_m, frame.frame_m.step_m, rank_m),
                [&](std::ostream& out) {
                    output_detail::piece(out, geometry_m, local, fields, frame, binary_m, scales_m);
                });
            return local;
        }

        /** @brief Rank-zero wrapper publication follows collective piece success. */
        void publishWrapper(const output_detail::FrameMetadata& frame) const {
            output_detail::atomicFile(
                directory_m / output_detail::frameName(basename_m, frame.frame_m.step_m),
                [&](std::ostream& out) {
                    output_detail::wrapper(out, geometry_m, extents_m, frame, basename_m, scales_m);
                });
        }

        /// Mirror the optional scalar material view and write only its owned cells.
        template <class Field, class View>
        void writeLocalMaterial(const Field& geometry, View epsilon,
                                const output_detail::Extent& local) const {
            const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), epsilon);
            output_detail::atomicFile(directory_m / ("material-" + std::to_string(rank_m) + ".vti"),
                                      [&](std::ostream& out) {
                                          output_detail::materialPiece(out, geometry_m, local, host,
                                                                       geometry.getNghost());
                                      });
        }

        /// Publish material ownership after every rank has completed its piece.
        void publishMaterial() const {
            if (rank_m == 0)
                output_detail::atomicFile(directory_m / "material.pvti", [&](std::ostream& out) {
                    output_detail::materialWrapper(out, geometry_m, extents_m);
                });
        }

        std::filesystem::path directory_m;
        OutputGeometry geometry_m;
        MPI_Comm communicator_m;
        std::string basename_m;
        int rank_m = 0;
        int size_m = 1;
        std::vector<output_detail::Extent> extents_m;
        std::vector<FrameRecord> frames_m;
        bool materialWritten_m = false;
        bool binary_m          = false;
        OutputScales scales_m;
    };
}  // namespace chdr::output
#endif
