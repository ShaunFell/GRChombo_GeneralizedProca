/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
/* #include "SimulationParametersBase.hpp" */

// Problem specific includes:
#include "KerrSchild.hpp"
#include "KerrQI.hpp"
#include "Potential.hpp"
#include "ProcaField.hpp"
#include "InitialProcaData.hpp"

class SimulationParameters : public ChomboParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
    {
        read_params(pp);
        check_params();
    }

    /// Read parameters from the parameter file
    void read_params(GRParmParse &pp)
    {
        //filenames
        pp.load("integrals_filename", integrals_filename);
        integrals_filename = output_path + "/" + integrals_filename;

        // Initial Kerr data
        pp.load("kerr_mass", kerrSchild_params.mass);
        pp.load("kerr_spin", kerrSchild_params.spin);
        pp.load("kerr_center", kerrSchild_params.center, center);
         pp.load("kerr_mass", kerrQI_params.mass);
        pp.load("kerr_spin", kerrQI_params.spin);
        pp.load("kerr_center", kerrQI_params.center, center);
        //pp.load("kerr_spindir", kerr_params.spin_direction);
        pp.load("KerrSchild", background_KerrSchild,true);
        pp.load("KerrQI", background_QI, false);

        //proca data
        pp.load("proca_mass", potential_params.mass);
        pp.load("proca_self_interaction", potential_params.self_interaction);

        pp.load("proca_damping", proca_params.vector_damping);

        pp.load("initial_proca_amplitude",initialdata_params.amplitude);

        //constants
        pp.load("G_Newton", G_Newton);

        pp.load("calculate_constraint_norms", calculate_constraint_norms, false);

        //extraction params
        pp.load("inner_r", inner_r, 0.0);
        pp.load("outer_r", outer_r, 200.0);

        //tagging
        pp.load("activate_ham_tagging", activate_ham_tagging, false);
        pp.load("activate_gauss_tagging", activate_gauss_tagging, false);

        //grid parameters
        pp.load("grid_scaling", grid_scaling, 1.);
        pp.load("nan_check", nan_check, 1);
        pp.load("sigma", sigma, 0.1);

        //AH Finder
        pp.load("AH_initial_guess", AH_initial_guess, 0.5*kerrSchild_params.mass);
        pp.load("AH_initial_guess", AH_initial_guess, 0.5*kerrQI_params.mass);
        pp.load("excision_with_AH", excise_with_AH, false);
    }

    void check_params()
    {
        warn_parameter("kerr_mass", kerrSchild_params.mass, kerrSchild_params.mass >= 0.0,
                       "should be >= 0.0");
        warn_parameter("kerr_mass", kerrQI_params.mass, kerrQI_params.mass >= 0.0,
                       "should be >= 0.0");
        
        check_parameter("KerrQI", background_QI, !(background_KerrSchild&&background_QI), "Either choose KerrSchild background or Quasi-Isotropic");
        
        if (background_QI){
            pout() << "Using Quasi-Isotropic background" << std::endl;
        } else if (background_KerrSchild){
            pout() << "Using Kerr-Schild background" << std::endl;
        }


        warn_parameter("inner_r", inner_r, inner_r != 0.0, "set to default parameters (0.0)");
        warn_parameter("outer_r", outer_r, outer_r != 200.0, "set to default parameter (200.0)");

        check_parameter("kerr_spin", kerrSchild_params.spin,
                        std::abs(kerrSchild_params.spin) <= kerrSchild_params.mass,
                        "must satisfy |a| <= M = " +
                            std::to_string(kerrSchild_params.mass));

        check_parameter("kerr_spin", kerrQI_params.spin,
                        std::abs(kerrQI_params.spin) <= kerrQI_params.mass,
                        "must satisfy |a| <= M = " +
                            std::to_string(kerrQI_params.mass));

        check_parameter("excision_with_AH", excise_with_AH, !(excise_with_AH && !AH_activate), "AH Finder must be turned on to use dynamical excision");

        check_parameter("grid_scaling", grid_scaling, grid_scaling>0, "Grid scaling parameter must be greater than zero");

        FOR(idir)
        {
            std::string name = "kerr_center[" + std::to_string(idir) + "]";
            warn_parameter(
                name, kerrSchild_params.center[idir],
                (kerrSchild_params.center[idir] >= 0) &&
                    (kerrSchild_params.center[idir] <= (ivN[idir] + 1) * coarsest_dx),
                "should be within the computational domain");
            
            warn_parameter(
                name, kerrQI_params.center[idir],
                (kerrQI_params.center[idir] >= 0) &&
                    (kerrQI_params.center[idir] <= (ivN[idir] + 1) * coarsest_dx),
                "should be within the computational domain");
        }
    }

    double G_Newton, outer_r, inner_r;
    double sigma;
    int nan_check;
    std::string integrals_filename;

    KerrSchild::params_t kerrSchild_params;
    KerrQI::params_t kerrQI_params;
    ProcaPotential::params_t potential_params;
    ProcaField<ProcaPotential>::params_t proca_params;

    init_params_t  initialdata_params;

    bool calculate_constraint_norms;
    bool activate_ham_tagging;
    bool activate_gauss_tagging;
    double grid_scaling;

    double AH_initial_guess;
    bool excise_with_AH;
    bool AH_activate;

    bool background_KerrSchild;
    bool background_QI;

};

#endif /* SIMULATIONPARAMETERS_HPP_ */
