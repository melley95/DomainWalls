#include "REAL.H"
#include <iomanip>
#include <iostream>

inline void writeFile(std::string name, Real value1, Real value2)
{
    int rank = 0;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#endif
    if (rank == 0)
    {
        ofstream myfile;
        myfile.open(name + ".txt", std::ios_base::app);

        // myfile << std::setprecision(5) << value1 << setw(20) <<
        // std::setprecision(5) << value2 << endl;
        myfile << std::fixed << std::scientific << std::setprecision(6)
               << std::left << std::setw(20) << value1 << std::left
               << std::setw(20) << value2 << std::endl;
    }
}

inline void openFile(std::string name)
{
    int rank = 0;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#endif
    if (rank == 0)
    {
        ofstream myfile;
        myfile.open(name + ".txt");
        // Write the headers
        myfile << "Ham error [%]" << std::setw(20) << "Mom error [%]" << endl;
    }
}
