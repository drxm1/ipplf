/** @file output_fixture.cpp
 * @brief Output-only sentinel using the production writer, no Maxwell step.
 * @details The 72-cell dataset is documented in ../ViewingGuide.org.
 * Coordinates and vector data are synthetic, not SI or Maxwell-solution claims.
 */
#include <Kokkos_MathematicalConstants.hpp>
#include <algorithm>

#include "VacuumOutput.h"

using Vector = Kokkos::Array<double, 3>;
using View   = Kokkos::View<Vector***, Kokkos::LayoutLeft>;

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

/** @brief Minimal IPPL-shaped index API; fixture ownership is specified independently. */
struct Range {
    int begin_m, end_m;
    int first() const { return begin_m; }
    int last() const { return end_m - 1; }
    std::size_t length() const { return static_cast<std::size_t>(end_m - begin_m); }
};

/** @brief Minimal IPPL-shaped layout API for a full-x/y, serial or two-z-piece grid. */
struct Layout {
    std::array<Range, 3> global_m{{{0, 4}, {0, 3}, {0, 6}}};
    std::array<Range, 3> local_m;
    const auto& getDomain() const { return global_m; }
    const auto& getLocalNDIndex() const { return local_m; }
};

/** @brief Real Kokkos allocation with one sentinel halo cell per side. */
struct FixtureField {
    View view_m;
    Layout layout_m;
    FixtureField(const char* name, int rank, int ranks)
        : view_m(name, 6, 5, 6 / ranks + 2)
        , layout_m{{{{0, 4}, {0, 3}, {0, 6}}},
                   {{{0, 4}, {0, 3}, {rank * (6 / ranks), (rank + 1) * (6 / ranks)}}}} {}
    const auto& getView() const { return view_m; }
    const auto& getLayout() const { return layout_m; }
    int getNghost() const { return 1; }
};

/** @brief Signed nextafter/pi values demand binary64 round trips, with a varying E offset. */
Vector fixtureValue(int i, int j, int k, double time, bool electric) {
    constexpr double Pi            = Kokkos::numbers::pi;
    constexpr double InfinityValue = std::numeric_limits<double>::infinity();
    if (electric)
        return {std::nextafter(1.0 + i + 100.0 * time, InfinityValue),
                std::nextafter(Pi * (10.0 + j), InfinityValue),
                std::nextafter((100.0 + k) / Pi, 0.0)};
    return {std::nextafter(-10.0 - 2.0 * i, -InfinityValue),
            std::nextafter(-Pi * (100.0 + 3.0 * j), -InfinityValue),
            std::nextafter((-1000.0 - 5.0 * k) / Pi, 0.0)};
}

/** @brief Explicit synthetic array times deliberately do not equal a writer-derived offset. */
chdr::vacuum::FrameRecord fixtureFrame(int step) {
    const double time     = step * 0.125;
    const double infinity = std::numeric_limits<double>::infinity();
    return {step, time, std::nextafter(time - 0.0625, -infinity),
            std::nextafter(time + 0.03125, infinity), 0.125};
}

/** @brief Fill every halo with a forbidden finite sentinel, then only owned cells. */
void fillFixture(FixtureField& field, double time, bool electric) {
    auto host       = Kokkos::create_mirror_view(field.view_m);
    const int start = field.layout_m.local_m[2].first();
    for (std::size_t k = 0; k < host.extent(2); ++k)
        for (std::size_t j = 0; j < host.extent(1); ++j)
            for (std::size_t i = 0; i < host.extent(0); ++i) {
                host(i, j, k) = {-987654321.0, -987654321.0, -987654321.0};
                if (i > 0 && i + 1 < host.extent(0) && j > 0 && j + 1 < host.extent(1) && k > 0
                    && k + 1 < host.extent(2))
                    host(i, j, k) = fixtureValue(i - 1, j - 1, start + k - 1, time, electric);
            }
    Kokkos::deep_copy(field.view_m, host);
}

/** @brief Three unequal actual times and steps, always writing the same writer path. */
void writeFixtureFrame(chdr::vacuum::VacuumOutput& writer, FixtureField& e, FixtureField& b,
                       int step) {
    fillFixture(e, step * 0.125, true);
    fillFixture(b, step * 0.125, false);
    writer.writeFrame(e, b, fixtureFrame(step));
    writer.writeSeries();
}

/** @brief Three unequal actual times and steps, using the same writer path. */
void successfulFrames(chdr::vacuum::VacuumOutput& writer, FixtureField& e, FixtureField& b) {
    for (const int step : {0, 1, 3})
        writeFixtureFrame(writer, e, b, step);
}

/** @brief Deliberately fail one rank's piece, or rank zero's manifest open. */
void injectFailure(const std::string& mode, const std::filesystem::path& directory,
                   FixtureField& electric, int rank, int ranks) {
    if (mode == "nonfinite" && rank == ranks - 1) {
        auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), electric.view_m);
        host(1, 1, 1)[0] = std::numeric_limits<double>::quiet_NaN();
        Kokkos::deep_copy(electric.view_m, host);
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

/** @brief Establish one good frame, then require coordinated publication failure. */
int rejectedFrame(chdr::vacuum::VacuumOutput& writer, FixtureField& e, FixtureField& b) {
    try {
        writer.writeFrame(e, b, fixtureFrame(1));
        writer.writeSeries();
    } catch (const std::runtime_error&) {
        return 1;
    }
    return 0;
}

/** @brief Establish one good frame, then require every rank to receive the injected failure. */
void failureFrame(chdr::vacuum::VacuumOutput& writer, FixtureField& e, FixtureField& b,
                  const std::filesystem::path& directory, const std::string& mode, int rank,
                  int ranks) {
    writeFixtureFrame(writer, e, b, 0);
    fillFixture(e, 0.125, true);
    fillFixture(b, 0.125, false);
    injectFailure(mode, directory, e, rank, ranks);
    const int caught = rejectedFrame(writer, e, b);
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
    if (argc < 2 || (ranks != 1 && ranks != 2))
        throw std::runtime_error("Use OUTPUT_DIR [piece|wrapper|series|nonfinite], 1/2 ranks");
    useCommaLocale();
    if (rank == 0)
        std::cout << "PASS custom comma-decimal global locale\n";
    const chdr::vacuum::OutputGeometry geometry{{4, 3, 6}, {-2, 3, 5}, {0.5, 1.25, 2}};
    chdr::vacuum::VacuumOutput writer(argv[1], geometry, MPI_COMM_WORLD);
    FixtureField e("fixture_E", rank, ranks), b("fixture_B", rank, ranks);
    if (argc == 2)
        successfulFrames(writer, e, b);
    else
        failureFrame(writer, e, b, argv[1], argv[2], rank, ranks);
    if (rank == 0 && argc == 2)
        std::cout << "PASS wrote three sentinel frames on " << ranks << " rank(s)\n";
}

/** @brief MPI/Kokkos fixture lifetime; unexpected errors use collective abort policy. */
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);
    try {
        runFixture(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
