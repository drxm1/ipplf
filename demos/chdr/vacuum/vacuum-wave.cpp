#include <exception>
#include <iostream>

#include "VacuumFailure.h"
#include "VacuumSetup.h"

namespace {
    int runVacuumWave();
}  // namespace

/** @brief Run the demo and finalize IPPL after simulation objects are destroyed.
 * @param argc Number of command-line arguments, including the program name.
 * @param argv Command-line strings passed to IPPL initialization.
 * @return Zero on success; initialization arguments are handled by IPPL.
 */
int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    const int result = runVacuumWave();
    ippl::finalize();
    return result;
}

namespace {
    /** @brief Own one initial reconstruction and one Standard solve.
     * @details For phi=0, raw E is at dt/2 and B at dt after solve; this follows
     * Fallahi, MITHRA 2.0, Eqs. (3.8)-(3.15), and FDTDSolverBase::solve/evaluate_EB
     * in FDTDSolverBase.hpp:23-26,109-127.
     * @return Zero on success; unexpected runtime exceptions abort MPI.
     */
    int runVacuumWave() {
        using namespace chdr::vacuum;
        try {
            auto fields = makePeriodicFieldContainer();
            checkLocalDomain(fields.getFL());
            initializeFieldStorage(fields);
            VacuumSolver solver(fields.getJ(), fields.getE(), fields.getB());
            initializePotentials(solver);
            solver.evaluate_EB();
            solver.solve();  // step(), timeShift(), evaluate_EB(): E at dt/2, B at dt.
            return 0;
        } catch (const IpplException& error) {
            return abortFailure(error.where(), error.what());
        } catch (const std::exception& error) {
            return abortFailure("runVacuumWave", error.what());
        }
    }
}  // namespace
