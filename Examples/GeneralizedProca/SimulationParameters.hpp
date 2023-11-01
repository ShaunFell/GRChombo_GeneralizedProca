/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "KerrBH.hpp"
#include "Potential.hpp"
#include "ProcaField.hpp"
#include "InitialProcaData.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        read_params(pp);
        check_params();
    }

    /// Read parameters from the parameter file
    void read_params(GRParmParse &pp)
    {

        // Initial Kerr data
        pp.load("kerr_mass", kerr_params.mass);
        pp.load("kerr_spin", kerr_params.spin);
        pp.load("kerr_center", kerr_params.center, center);
        pp.load("kerr_spindir", kerr_params.spin_direction);

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
    }

    void check_params()
    {
        warn_parameter("kerr_mass", kerr_params.mass, kerr_params.mass >= 0.0,
                       "should be >= 0.0");
        check_parameter("kerr_spin", kerr_params.spin,
                        std::abs(kerr_params.spin) <= kerr_params.mass,
                        "must satisfy |a| <= M = " +
                            std::to_string(kerr_params.mass));
        FOR(idir)
        {
            std::string name = "kerr_center[" + std::to_string(idir) + "]";
            warn_parameter(
                name, kerr_params.center[idir],
                (kerr_params.center[idir] >= 0) &&
                    (kerr_params.center[idir] <= (ivN[idir] + 1) * coarsest_dx),
                "should be within the computational domain");
        }
    }

    double G_Newton, outer_r, inner_r;

    KerrBH::params_t kerr_params;
    ProcaPotential::params_t potential_params;
    ProcaField<ProcaPotential>::params_t proca_params;
    InitialProcaData::init_params_t initialdata_params;

    bool calculate_constraint_norms;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
