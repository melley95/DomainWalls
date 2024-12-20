#include "mpi.h"
#include <iostream>

#include "CTTK.hpp"
#include "CTTKHybrid.hpp"
#include "GRSolver.hpp"
#include "ScalarField.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    int status = 0;
#ifdef _OPENMP
    std::cout << "#threads = " << omp_get_max_threads() << std::endl;
#endif
#ifdef CH_MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        cout << "Running with MPI" << endl;
#endif

    if (argc < 1)
    {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
    }

    // Read params here for now

    // Create solver, specifying method and matter
    char *inFile = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, inFile);

#if defined(USE_CTTK)
    GRSolver<CTTK<ScalarField>, ScalarField> solver(pp);
    pout() << "Using CTTK" << endl;
#elif defined(USE_CTTKHybrid)
    GRSolver<CTTKHybrid<ScalarField>, ScalarField> solver(pp);
    pout() << "Using CTTK Hybrid" << endl;
#elif defined(USE_OtherMethod)
    GRSolver<OtherMethod<ScalarField>, ScalarField> solver(pp);
    pout() << "Using Other Method" << endl;
#else
#error                                                                         \
    "No valid method defined. Please define one of the methods in the GNUMakefile."
#endif

    solver.setup();
    status = solver.run();

#ifdef CH_MPI
    MPI_Finalize();
#endif
    return status;
}