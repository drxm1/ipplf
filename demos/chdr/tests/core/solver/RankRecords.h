/** @file RankRecords.h
 * @brief Owned-cell snapshots and independent reductions for decomposition comparisons.
 * @details A snapshot is an ordered list of values, not a rank-local image.
 * Every rank copies its field view to a host mirror, packs owned cells only,
 * and sends plain doubles with their global component indices. Rank zero
 * restores x/y/z/component order and checks ownership before writing a record.
 */
#pragma once
#include "Ippl.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "OwnedStatistics.h"
#include "core/base/SpatialDimension.h"

static_assert(chdr::solver::SpatialDim == 3,
              "RankRecords.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    namespace rankDetail {
        struct PackedValues {
            std::vector<int> indices_m;
            std::vector<double> values_m;
        };
        struct Counts {
            std::vector<int> counts_m, offsets_m;
            int total_m = 0;
        };

        /// Reject oversized snapshots collectively before narrowing to MPI's int counts.
        template <class FieldType>
        std::size_t checkedGlobalSize(const FieldType& field) {
            const auto domain = field.getLayout().getDomain();
            std::size_t count = FieldType::value_type::dim;
            for (unsigned d = 0; d < SpatialDim; ++d) {
                const auto length = std::size_t(domain[d].length());
                if (length > std::size_t(std::numeric_limits<int>::max()) / count)
                    throw std::length_error("Rank snapshot exceeds MPI int count capacity");
                count *= length;
            }
            return count;
        }

        /// Append a cell in component order; no arithmetic is applied to its values.
        template <class CellValue>
        void appendCell(PackedValues& packed, int cell, const CellValue& value) {
            for (unsigned c = 0; c < CellValue::dim; ++c) {
                packed.indices_m.push_back(cell * CellValue::dim + c);
                packed.values_m.push_back(value[c]);
            }
        }

        /// Convert local owned indices into the global flattened cell index.
        inline int globalCell(int i, int j, int k, const ippl::NDIndex<SpatialDim>& local,
                              const ippl::NDIndex<SpatialDim>& global) {
            return ((i + local[0].first() - global[0].first()) * global[1].length() + j
                    + local[1].first() - global[1].first())
                       * global[2].length()
                   + k + local[2].first() - global[2].first();
        }

        /** @brief Copy the view before host access and omit every ghost entry.
         * @details The copied mirror owns its bytes while the source may be in a
         * different memory space. Integer global indexing makes the output order
         * independent of the layout's split direction and rank count.
         */
        template <class FieldType>
        PackedValues packOwned(const FieldType& field) {
            auto host = field.getHostMirror();
            Kokkos::deep_copy(host, field.getView());
            const auto local  = field.getLayout().getLocalNDIndex();
            const auto global = field.getLayout().getDomain();
            const int ghost   = field.getNghost();
            PackedValues packed;
            for (int i = 0; i < static_cast<int>(local[0].length()); ++i)
                for (int j = 0; j < static_cast<int>(local[1].length()); ++j)
                    for (int k = 0; k < static_cast<int>(local[2].length()); ++k)
                        appendCell(packed, globalCell(i, j, k, local, global),
                                   host(i + ghost, j + ghost, k + ghost));
            return packed;
        }

        /// Gather message lengths before allocating rank zero's receive buffers.
        inline Counts gatherCounts(int count, MPI_Comm comm) {
            int size;
            MPI_Comm_size(comm, &size);
            Counts result{std::vector<int>(size), std::vector<int>(size)};
            MPI_Gather(&count, 1, MPI_INT, result.counts_m.data(), 1, MPI_INT, 0, comm);
            for (int rank = 0; rank < size; ++rank) {
                result.offsets_m[rank] = result.total_m;
                result.total_m += result.counts_m[rank];
            }
            return result;
        }

        /// Both arrays use identical counts so indices remain paired with their values.
        inline PackedValues collect(const PackedValues& local, const Counts& counts,
                                    MPI_Comm comm) {
            PackedValues result{std::vector<int>(counts.total_m),
                                std::vector<double>(counts.total_m)};
            const int count = static_cast<int>(local.values_m.size());
            MPI_Gatherv(local.indices_m.data(), count, MPI_INT, result.indices_m.data(),
                        counts.counts_m.data(), counts.offsets_m.data(), MPI_INT, 0, comm);
            MPI_Gatherv(local.values_m.data(), count, MPI_DOUBLE, result.values_m.data(),
                        counts.counts_m.data(), counts.offsets_m.data(), MPI_DOUBLE, 0, comm);
            return result;
        }

        /// A missing, repeated or out-of-range global index is an ownership error.
        inline std::vector<double> orderValues(const PackedValues& packed, std::size_t expected) {
            const auto total = packed.values_m.size();
            if (total != expected)
                throw std::runtime_error("Gather count mismatch");
            std::vector<double> ordered(total, std::numeric_limits<double>::quiet_NaN());
            std::vector<bool> seen(total, false);
            for (std::size_t n = 0; n < total; ++n) {
                const int index = packed.indices_m[n];
                if (index < 0 || static_cast<std::size_t>(index) >= total || seen[index])
                    throw std::runtime_error("Gather ownership mismatch");
                ordered[index] = packed.values_m[n];
                seen[index]    = true;
            }
            return ordered;
        }

        /// Shared stable statistics retain the existing record ordering and MPI definitions.
        template <class FieldType>
        std::array<double, 3> statistics(const FieldType& field) {
            return ownedStatistics(field).values();
        }
    }  // namespace rankDetail

    /** @brief Return globally ordered owned values on rank zero, empty elsewhere.
     * @pre All ranks of the field communicator participate in matching order.
     * The global component count must fit MPI's int counts; oversized requests
     * are rejected on every rank before any gather or narrowing conversion.
     * @details No complex MPI datatype or field norm implementation is used.
     * Statistics are gathered independently, rather than recomputed from this
     * host array, so the tests compare both field values and the actual parallel reductions.
     * @tparam FieldType IPPL field with component-indexable cell values of any quantity.
     * @param field Read-only source; halo values are omitted and no time level is assumed.
     */
    template <class FieldType>
    std::vector<double> gatherOwned(const FieldType& field) {
        const auto comm     = field.getLayout().comm;
        const auto expected = rankDetail::checkedGlobalSize(field);
        const auto local    = rankDetail::packOwned(field);
        const auto counts = rankDetail::gatherCounts(static_cast<int>(local.values_m.size()), comm);
        const auto gathered = rankDetail::collect(local, counts, comm);
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank != 0)
            return {};
        return rankDetail::orderValues(gathered, expected);
    }

    /// The CTest fixture supplies distinct fresh directories for each rank count.
    inline std::filesystem::path recordDirectory() {
        const char* directory = std::getenv("CHDR_SOLVER_RECORD_DIR");
        return directory ? directory : "runs";
    }

    /// Write only on rank zero after all collectives have completed.
    inline void writeFieldRecord(const std::string& name, const std::vector<double>& values,
                                 const std::array<double, 3>& statistics) {
        std::filesystem::create_directories(recordDirectory());
        std::ofstream out(recordDirectory() / (name + ".json"));
        out << nlohmann::json{{"values", values},
                              {"max", statistics[0]},
                              {"sum_abs", statistics[1]},
                              {"sum_squares", statistics[2]}}
                   .dump()
            << '\n';
        if (!out)
            throw std::runtime_error("Cannot write rank record");
    }

    /** @brief Write one global field and independently reduced statistics on rank zero.
     * @pre Every rank in the field communicator calls in the same order.
     * @details gatherOwned orders owned cells by global index; ownedStatistics
     * uses its separate device reductions, so neither measurement hides the
     * decomposition sensitivity of the other. Ghosts are excluded from both.
     */
    template <class FieldType>
    void recordField(const std::string& name, const FieldType& field) {
        const auto statistics = rankDetail::statistics(field);
        const auto values     = gatherOwned(field);
        int rank;
        MPI_Comm_rank(field.getLayout().comm, &rank);
        if (rank == 0)
            writeFieldRecord(name, values, statistics);
    }

    /** @brief Record both histories and both outputs to expose time-shift or halo defects.
     * @tparam SolverType Solver exposing potential histories and E/B output accessors.
     * @param name Record prefix shared by the four saved fields.
     * @param solver Current state; A_n and A_nm1 are saved under current/previous
     * names, with their physical times determined by the solver's LevelContract.
     */
    template <class SolverType>
    void recordState(const std::string& name, SolverType& solver) {
        recordField(name + "-current", solver.A_n);
        recordField(name + "-previous", solver.A_nm1);
        recordField(name + "-E", solver.electricField());
        recordField(name + "-B", solver.magneticField());
    }
}  // namespace chdr::solver::test
