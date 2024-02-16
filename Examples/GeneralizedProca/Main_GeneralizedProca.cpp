/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "parstream.H" //Gives us pout()

// System includes
#include <iostream>

// Our general includes
#include "DefaultLevelFactory.hpp"
#include "BHAMR.hpp"
#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "ProcaFieldLevel.hpp"

// Chombo namespace
#include "UsingNamespace.H"

int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
        return 0;

    // The line below selects the problem that is simulated
    // (To simulate a different problem, define a new child of AMRLevel
    // and an associated LevelFactory)
    BHAMR bh_amr;
    DefaultLevelFactory<ProcaFieldLevel> proca_field_level_fact(bh_amr, sim_params);
    setupAMRObject(bh_amr, proca_field_level_fact);

    //setup interpolating object
    AMRInterpolator<Lagrange<4>> interpolator(
        bh_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    bh_amr.set_interpolator(&interpolator); // also sets puncture_tracker interpolator

    using Clock = std::chrono::steady_clock;
    using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;

    std::chrono::time_point<Clock> start_time = Clock::now();

#ifdef USE_AHFINDER //Chombo make flags
    if (sim_params.AH_activate)
    {
        pout() << "Defining surface geometry"<<endl;
        AHSurfaceGeometry sph(sim_params.kerr_params.center);

#ifdef USE_CHI_CONTOURS //located in UserVariables
        std::string str_chi = std::to_string(
            sim_params.AH_params.func_params.look_for_chi_contour);
        sim_params.AH_params.stats_prefix = "stats_chi_" + str_chi + "_";
        sim_params.AH_params.coords_prefix = "coords_chi_" + str_chi + "_";
        bh_amr.m_ah_finder.add_ah(sph, sim_params.AH_initial_guess, sim_params.AH_params);
#else 
        pout() << "adding apparent horizon"<<endl;
        bh_amr.m_ah_finder.add_ah(sph, sim_params.AH_initial_guess, sim_params.AH_params);
#endif //USE_CHI_CONTOURS
    }
#endif //USE_AHFINDER


    //call the PostTimeStep right now!!!!
    auto task = [](GRAMRLevel *level)
    {
        if (level->time() == 0.)
            level->specificPostTimeStep();
    };

    // call 'now' really now
    MultiLevelTaskPtr<> call_task(task);
    call_task.execute(bh_amr); 

    //go go go !!!!!! Run simulation
    bh_amr.run(sim_params.stop_time, sim_params.max_steps);
 
    auto now = Clock::now();
    auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
    pout() << "Total simulation time (mins): " << duration.count() << ".\n";

    bh_amr.conclude();

    CH_TIMER_REPORT(); // Report results when running with Chombo timers.

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChombo(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
