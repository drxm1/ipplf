/** @file output_fixture.cpp
 * @brief Output-only sentinel using the production writer, no Maxwell step.
 * @details The 72-cell dataset is documented in demos/chdr/apps/vacuum-wave/ViewingGuide.org.
 * Coordinates and vector data are synthetic, not SI or Maxwell-solution claims.
 *
 * The writer is core/output/FieldOutput.h, which all three programs use: chdr-vacuum-wave
 * directly, planar-interface and interface-movies through opt-in/output/FieldRunOutput.h.
 * The fixture has 4 by 3 by 6 cells, origin (-2, 3, 5) and spacing (0.5, 1.25, 2); on two
 * ranks each rank owns 3 z layers. Every halo cell holds the sentinel -987654321.
 * Under a comma-decimal global locale, the case "valid" writes three vector arrays at
 * steps 0, 1 and 3 (series times 0, 0.125, 0.375). R2 supplies staggered positions and
 * separate electric/averaged-magnetic and raw-magnetic sample times. The synthetic
 * averaged-magnetic values have an independent offset and output scale to distinguish
 * all three arrays. Serial readback exercises ASCII; two-rank readback exercises binary.
 * The failure cases "piece", "wrapper", "series" and
 * "nonfinite" write frame 0 and then make frame 1 fail on one rank; every rank must catch
 * std::runtime_error. read_output_fixture.py reads the files.
 *
 * CTest (through run_output_fixture.py): chdr.vacuum.output.serial and
 * chdr.vacuum.output.mpi2 (case "valid"), and chdr.vacuum.output.failure.piece, .wrapper,
 * .series and .nonfinite on two ranks.
 */
#include <Kokkos_MathematicalConstants.hpp>

#include "core/output/FieldOutput.h"

using Vector = chdr::EMFieldsType::value_type;

/** @brief Test locale independent of the host's installed locale catalogue. */
struct CommaPunctuation : std::numpunct<char> {
    char do_decimal_point() const override { return ','; }
};

/** @brief Prove newly constructed streams inherit a comma before testing writer overrides. */
void useCommaLocale() {
    std::locale::global(std::locale(std::locale::classic(), new CommaPunctuation));
    std::ostringstream probe;
    probe << 1.5;
    if (probe.str() != "1,5")
        throw std::runtime_error("Comma-decimal locale probe failed");
}

/** @brief Own real IPPL storage at a stable address, split only in the z direction.
 * @details Synthetic values and expected global indices remain independent of the writer.
 * No solver or boundary exchange is used; all halo entries remain forbidden sentinels.
 */
struct FixtureFields {
    ippl::NDIndex<3> domain_m{ippl::Index(4), ippl::Index(3), ippl::Index(6)};
    Mesh_t<3> mesh_m{domain_m, Vector(0.5, 1.25, 2.0), Vector(-2.0, 3.0, 5.0)};
    ippl::FieldLayout<3> layout_m{MPI_COMM_WORLD, domain_m, {false, false, true}, false};
    chdr::EMFieldsType electric_m, rawMagnetic_m, averagedMagnetic_m;

    FixtureFields() {
        electric_m.initialize(mesh_m, layout_m);
        rawMagnetic_m.initialize(mesh_m, layout_m);
        averagedMagnetic_m.initialize(mesh_m, layout_m);
    }
};

/** @brief Signed nextafter/pi values demand binary64 round trips, with a varying E offset. */
Vector fixtureValue(int i, int j, int k, double time, unsigned channel) {
    constexpr double Pi            = Kokkos::numbers::pi;
    constexpr double InfinityValue = std::numeric_limits<double>::infinity();
    if (channel == 0)
        return {std::nextafter(1.0 + i + 100.0 * time, InfinityValue),
                std::nextafter(Pi * (10.0 + j), InfinityValue),
                std::nextafter((100.0 + k) / Pi, 0.0)};
    Vector value{std::nextafter(-10.0 - 2.0 * i, -InfinityValue),
                 std::nextafter(-Pi * (100.0 + 3.0 * j), -InfinityValue),
                 std::nextafter((-1000.0 - 5.0 * k) / Pi, 0.0)};
    if (channel == 2)
        for (unsigned c = 0; c < 3; ++c)
            value[c] += 0.25;
    return value;
}

/** @brief Supply only the series clock; native array times come from the declared contract. */
chdr::output::FrameRecord fixtureFrame(int step) {
    return {step, step * 0.125, 0.125};
}

/** @brief Fill every halo with a forbidden finite sentinel, then only owned cells. */
void fillFixture(chdr::EMFieldsType& field, double time, unsigned channel) {
    auto host       = Kokkos::create_mirror_view(field.getView());
    const int start = field.getLayout().getLocalNDIndex()[2].first();
    for (std::size_t k = 0; k < host.extent(2); ++k)
        for (std::size_t j = 0; j < host.extent(1); ++j)
            for (std::size_t i = 0; i < host.extent(0); ++i) {
                host(i, j, k) = {-987654321.0, -987654321.0, -987654321.0};
                if (i > 0 && i + 1 < host.extent(0) && j > 0 && j + 1 < host.extent(1) && k > 0
                    && k + 1 < host.extent(2))
                    host(i, j, k) = fixtureValue(i - 1, j - 1, start + k - 1, time, channel);
            }
    Kokkos::deep_copy(field.getView(), host);
}

/// Fill three independently distinguishable arrays without any Maxwell update.
void fillFixture(FixtureFields& fields, int step) {
    fillFixture(fields.electric_m, step * 0.125, 0);
    fillFixture(fields.rawMagnetic_m, step * 0.125, 1);
    fillFixture(fields.averagedMagnetic_m, step * 0.125, 2);
}

/// Bind the three arrays to native quantity labels and one nontrivial conversion factor.
auto fixtureSnapshots(const FixtureFields& fields) {
    using chdr::output::VectorFieldSnapshot;
    return std::array{
        VectorFieldSnapshot{"E", fields.electric_m, chdr::Quantity::Electric, 1.0},
        VectorFieldSnapshot{"B_raw", fields.rawMagnetic_m, chdr::Quantity::RawMagnetic, 1.0},
        VectorFieldSnapshot{"B_R2", fields.averagedMagnetic_m, chdr::Quantity::Magnetic, 0.5}};
}

/// Exercise native staggered R2 metadata without writing a material map first.
void writeFixture(chdr::output::FieldOutput& writer, const FixtureFields& fields, int step) {
    const chdr::LevelContract contract{chdr::EBReconstructionConvention::R2,
                                       chdr::SpatialDiscretization::Staggered};
    const auto snapshots = fixtureSnapshots(fields);
    writer.writeFrame(fixtureFrame(step), contract, snapshots);
    writer.writeSeries();
}

/** @brief Three unequal steps 0, 1 and 3 (times 0, 0.125, 0.375), on the same writer path. */
void successfulFrames(chdr::output::FieldOutput& writer, FixtureFields& fields) {
    for (const int step : {0, 1, 3}) {
        fillFixture(fields, step);
        writeFixture(writer, fields, step);
    }
}

/** @brief Prepare one failure: a NaN in one owned E value on the last rank ("nonfinite"),
 * or a directory that blocks the temporary file of the last rank's piece ("piece"), of
 * rank zero's PVTI wrapper ("wrapper") or of rank zero's series file ("series"). */
void injectFailure(const std::string& mode, const std::filesystem::path& directory,
                   chdr::EMFieldsType& electric, int rank, int ranks) {
    if (mode == "nonfinite" && rank == ranks - 1) {
        auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), electric.getView());
        host(1, 1, 1)[0] = std::numeric_limits<double>::quiet_NaN();
        Kokkos::deep_copy(electric.getView(), host);
    }
    if (mode == "piece" && rank == ranks - 1)
        std::filesystem::create_directory(directory
                                          / ("wave-1-rank-" + std::to_string(rank) + ".vti.tmp"));
    if (mode == "wrapper" && rank == 0)
        std::filesystem::create_directory(directory / "wave-1.pvti.tmp");
    if (mode == "series" && rank == 0)
        std::filesystem::create_directory(directory / "wave.pvti.series.tmp");
    MPI_Barrier(MPI_COMM_WORLD);
}

/** @brief Try to publish frame 1; return 1 if this rank catches std::runtime_error, else 0. */
int rejectedFrame(chdr::output::FieldOutput& writer, const FixtureFields& fields) {
    try {
        writeFixture(writer, fields, 1);
    } catch (const std::runtime_error&) {
        return 1;
    }
    return 0;
}

/** @brief Establish one good frame, then require every rank to receive the injected failure. */
void failureFrame(chdr::output::FieldOutput& writer, FixtureFields& fields,
                  const std::filesystem::path& directory, const std::string& mode, int rank,
                  int ranks) {
    fillFixture(fields, 0);
    writeFixture(writer, fields, 0);
    fillFixture(fields, 1);
    injectFailure(mode, directory, fields.electric_m, rank, ranks);
    const int caught = rejectedFrame(writer, fields);
    int total        = 0;
    MPI_Allreduce(&caught, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (total != ranks)
        throw std::runtime_error("Failure did not reach every rank");
    if (rank == 0)
        std::cout << "PASS coordinated " << mode << " failure on " << ranks << " rank(s)\n";
}

/** @brief Dispatch bounded fixture modes; output directory must initially be empty. */
void runFixture(int argc, char** argv) {
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    if (argc != 4 || (ranks != 1 && ranks != 2))
        throw std::runtime_error("Use OUTPUT_DIR MODE ascii|binary|binary32, with 1/2 ranks");
    useCommaLocale();
    if (rank == 0)
        std::cout << "PASS custom comma-decimal global locale\n";
    const chdr::output::OutputGeometry geometry{{4, 3, 6}, {-2, 3, 5}, {0.5, 1.25, 2}};
    const std::string mode = argv[2], encoding = argv[3];
    chdr::output::FieldOutput writer(argv[1], geometry, MPI_COMM_WORLD, "wave", encoding != "ascii",
                                     {}, encoding == "binary32");
    FixtureFields fields;
    if (mode == "valid")
        successfulFrames(writer, fields);
    else
        failureFrame(writer, fields, argv[1], mode, rank, ranks);
    if (rank == 0 && mode == "valid")
        std::cout << "PASS wrote three sentinel frames on " << ranks << " rank(s)\n";
}

/** @brief IPPL fixture lifetime; unexpected errors use collective abort policy. */
int main(int argc, char** argv) {
    ippl::initialize(argc, argv);
    try {
        runFixture(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    ippl::finalize();
    return 0;
}
