#ifndef SIMULATIONPARAMETERS_HPP_
#error "This file should only be included through SimulationParameters.hpp"
#endif

#include "CoarseAverage.H"
#include "FilesystemTools.hpp"
#include "GRParmParse.hpp"
#include "ProblemDomain.H"
#include "REAL.H"
#include "RealVect.H"

#include <fstream>
#include <iostream>

template <class method_t, class matter_t>
struct SimulationParameters<method_t, matter_t>::BaseParams
{

    Real HamMom_tolerance;
    int max_NL_iter;
    bool write_iter;
    Real iter_tolerance;
    int max_iter;
    int numMGIter;
    int numMGSmooth;
    int preCondSolverDepth;
    Real alpha;
    Real beta;
    bool readin_matter_data;
    std::string input_filename;
    std::string output_filename;
    std::string output_path;
    std::string pout_filename;
    std::string pout_path;
    std::string error_filename;
    Real G_Newton;
    int verbosity;
};

template <class method_t, class matter_t>
void SimulationParameters<method_t, matter_t>::read_base_params(GRParmParse &pp)
{
#ifdef CH_MPI
    // setPoutBaseName must be called before pout is used
    if (pp.contains("pout_path"))
    {
        pp.get("pout_path", base_params.pout_path);
        // Create pout directory
        if (!FilesystemTools::directory_exists(base_params.pout_path))
            FilesystemTools::mkdir_recursive(base_params.pout_path);
        setPoutBaseName(base_params.pout_path + "pout");
    }
    else
    {
        base_params.pout_path = "";
    }
#endif

    pp.query("HamMom_tolerance", base_params.HamMom_tolerance);
    pp.query("max_NL_iterations", base_params.max_NL_iter);
    pp.query("write_iter", base_params.write_iter);

    // Setup multigrid params, most of them defaulted
    if (pp.contains("iter_tolerance"))
    {
        pp.query("iter_tolerance", base_params.iter_tolerance);
    }
    else
    {
        base_params.iter_tolerance = 1e-10;
    }

    if (pp.contains("max_iterations"))
    {
        pp.query("max_iterations", base_params.max_iter);
    }
    else
    {
        base_params.max_iter = 100;
    }

    if (pp.contains("numMGIterations"))
    {
        pp.query("numMGIterations", base_params.numMGIter);
    }
    else
    {
        base_params.numMGIter = 4;
    }

    if (pp.contains("numMGsmooth"))
    {
        pp.query("numMGsmooth", base_params.numMGSmooth);
    }
    else
    {
        base_params.numMGSmooth = 4;
    }

    if (pp.contains("preCondSolverDepth"))
    {
        pp.query("preCondSolverDepth", base_params.preCondSolverDepth);
    }
    else
    {
        base_params.preCondSolverDepth = -1;
    }

    // Params for variable coefficient multigrid solver, solving the eqn
    // alpha*aCoef(x)*I - beta*bCoef(x) * laplacian = rhs
    // spatially-varying aCoef and bCoef are set in Methods
    // (for pure laplacian, alpha = 0, beta=-1)
    if (pp.contains("alpha")){
        pp.get("alpha", base_params.alpha);
    }
    else{
        base_params.alpha = 1.0;
    }
    if (pp.contains("beta")){
        pp.get("beta", base_params.beta);
    }
    else{
        base_params.beta = -1.0;
    }

    // Read from hdf5 file
    if (pp.contains("input_filename"))
    {
        pp.get("input_filename", base_params.input_filename);
        base_params.readin_matter_data = true;
    }
    else
    {
        base_params.input_filename = "";
        base_params.readin_matter_data = false;
    }
    if (pp.contains("output_path"))
    {
        pp.get("output_path", base_params.output_path);
        if (!FilesystemTools::directory_exists(base_params.output_path))
            FilesystemTools::mkdir_recursive(base_params.output_path);
    }
    else
    {
        base_params.output_path = "";
    }

    if (pp.contains("output_filename"))
    {
        string filename;
        pp.get("output_filename", filename);
        base_params.output_filename = base_params.output_path + filename;
    }
    else
    {
        base_params.output_filename =
            base_params.output_path + "OutputDataFinal.3d.hdf5";
    }

    pp.get("G_Newton", base_params.G_Newton);
    base_params.verbosity = 3;
    pp.query("verbosity", base_params.verbosity);
}