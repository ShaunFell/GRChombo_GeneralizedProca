

//general includes common to most GR problems
#include "ProcaFieldLevel.hpp"
#include "AMRReductions.hpp"
#include "ComputePack.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"



//RHS Update
#include "FixedBGEvolution.hpp"


//cell tagging
#include "FixedGridsTaggingCriterion.hpp"

//problem specific includes
#include "InitialProcaData.hpp"
#include "KerrSchild.hpp"
#include "KerrQI.hpp"
#include "Potential.hpp"
#include "ProcaField.hpp"
#include "Diagnostics.hpp"
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"


/* #include "MatterWeyl4.hpp"
#include "WeylExtraction.hpp"  */


//do things at end of advance step, after RK4 calculation
void ProcaFieldLevel::specificAdvance()
{
    //check for nans
    if (m_p.nan_check){
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                        EXCLUDE_GHOST_CELLS, disable_simd());
    }
};

//initialize data for the field
void ProcaFieldLevel::initialData()
{
    CH_TIME("ProcaFieldLevel::initialData");
    if (m_verbosity){
        pout()<<"ProcaFieldLevel::initialData " << m_level << endl;
    }

    if (m_p.background_KerrSchild) {
        KerrSchild kerr_init { m_p.kerrSchild_params, m_dx };
        InitialProcaData<KerrSchild> proca_init { m_p.initialdata_params, m_p.potential_params, m_p.kerrSchild_params, m_dx };
        BoxLoops::loop(kerr_init, m_state_new, m_state_new, FILL_GHOST_CELLS);
        BoxLoops::loop(proca_init,m_state_new, m_state_new, FILL_GHOST_CELLS);

        BoxLoops::loop(
                ExcisionProcaEvolution<ProcaFieldWithPotential, KerrSchild>(kerr_init, m_dx, m_p.center, m_p.inner_r),
                m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd()
            );

    } else if (m_p.background_QI) {
        KerrQI kerr_init { m_p.kerrQI_params, m_dx };
        InitialProcaData<KerrQI> proca_init { m_p.initialdata_params, m_p.potential_params, m_p.kerrQI_params, m_dx };
        BoxLoops::loop(kerr_init, m_state_new, m_state_new, FILL_GHOST_CELLS);
        BoxLoops::loop(proca_init,m_state_new, m_state_new, FILL_GHOST_CELLS);

        BoxLoops::loop(
                ExcisionProcaEvolution<ProcaFieldWithPotential, KerrQI>(kerr_init, m_dx, m_p.center, m_p.inner_r),
                m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd()
            );
    }


};

#ifdef CH_USE_HDF5
//things to do before outputting a checkpoint file
void ProcaFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    ProcaPotential potential(m_p.potential_params);
    ProcaFieldWithPotential proca_field(potential, m_p.proca_params);
    ProcaConstraint<ProcaPotential> proca_constraint(m_dx, m_p.potential_params.mass, m_p.proca_params.vector_damping, potential);
    EffectiveMetric<ProcaPotential> proca_eff_met(m_dx, m_p.potential_params.mass, m_p.proca_params.vector_damping, potential);
    ProcaSquared Asquared(m_dx);

    //compute diagnostics on each cell of current level
    BoxLoops::loop(
        make_compute_pack(
            Asquared,
            proca_constraint,
            proca_eff_met
            ),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
        );

    if (m_p.background_KerrSchild) {
        KerrSchild kerr_bh { m_p.kerrSchild_params, m_dx };
        BoxLoops::loop(
            ExcisionDiagnostics<KerrSchild>(kerr_bh, m_dx, m_p.center, m_p.inner_r),
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd()
        );
    } else if (m_p.background_QI) {
        KerrQI kerr_bh { m_p.kerrQI_params, m_dx };
        BoxLoops::loop(
            ExcisionDiagnostics<KerrQI>(kerr_bh, m_dx, m_p.center, m_p.inner_r),
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd()
        );
    }
};
#endif //CH_USE_HDF5

//RHS routines used at each RK4 step
void ProcaFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                const double a_time)
{

    //Calculate right hand side with matter_t = ProcaField
    ProcaPotential potential(m_p.potential_params);
    ProcaFieldWithPotential proca_field(potential, m_p.proca_params);
    if (m_p.background_KerrSchild) {
        KerrSchild kerr_init { m_p.kerrSchild_params, m_dx };
        FixedBGEvolution<ProcaFieldWithPotential, KerrSchild> matter_class(proca_field, kerr_init, m_p.sigma, m_dx, m_p.center);
        BoxLoops::loop(matter_class, a_soln, a_rhs, SKIP_GHOST_CELLS);
        BoxLoops::loop(
            ExcisionProcaEvolution<ProcaFieldWithPotential, KerrSchild>(kerr_init, m_dx, m_p.center, m_p.inner_r),
            a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd()
        );
    
    } else if (m_p.background_QI) {
        KerrQI kerr_init { m_p.kerrQI_params, m_dx };
        FixedBGEvolution<ProcaFieldWithPotential, KerrQI> matter_class(proca_field, kerr_init, m_p.sigma, m_dx, m_p.center);
        BoxLoops::loop(matter_class, a_soln, a_rhs, SKIP_GHOST_CELLS);
        BoxLoops::loop(
            ExcisionProcaEvolution<ProcaFieldWithPotential, KerrQI>(kerr_init, m_dx, m_p.center, m_p.inner_r),
            a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd()
        );
    }
};

void ProcaFieldLevel::specificUpdateODE(GRLevelData &a_soln, const GRLevelData &a_rhs,
                                Real a_dt)
{
};


//things to do before tagging cells (e.g. filling ghosts)
void ProcaFieldLevel::preTagCells()
{
};


//compute tagging criteria for grid
void ProcaFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                             const FArrayBox &current_state,
                                             const FArrayBox &current_state_diagnostics)
{
    CH_TIME("ProcaFieldLevel::computeTaggingCriterion");

    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level,
                                                    m_p.grid_scaling*m_p.L, m_p.center),
                       current_state, tagging_criterion, disable_simd());
}


void ProcaFieldLevel::specificPostTimeStep()
{
    CH_TIME("ProcaFieldLevel::specificPostTimeStep");

    bool first_step = (m_time == m_dt); //is this the first call of posttimestep?

    int min_level = 0;
    bool calculate_charges = at_level_timestep_multiple(min_level);

    //calculate energy and angular momentum densities
    if (calculate_charges)
    {
        fillAllEvolutionGhosts();

        //setup proca field
        ProcaPotential potential(m_p.potential_params);
        ProcaFieldWithPotential proca_field(potential, m_p.proca_params);

        if (m_p.background_KerrSchild) {
            KerrSchild kerr_init { m_p.kerrSchild_params, m_dx };

            //compute energy and angular momentum
            EnergyAndAngularMomentum<ProcaFieldWithPotential, KerrSchild> energy_and_angularmomentum(kerr_init, m_dx, proca_field, m_p.center);

            //compute background first
            BoxLoops::loop(kerr_init, m_state_new, m_state_new, FILL_GHOST_CELLS);
            
            //compute diagnostics on each cell of current level
            BoxLoops::loop(
                energy_and_angularmomentum,
                m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
            );

            //excise within horizon
            BoxLoops::loop(
                ExcisionDiagnostics<KerrSchild>(kerr_init, m_dx, m_p.center, m_p.inner_r),
                m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
                disable_simd()
            );
        } else if (m_p.background_QI) {
            KerrQI kerr_init { m_p.kerrQI_params, m_dx };

            //compute energy and angular momentum
            EnergyAndAngularMomentum<ProcaFieldWithPotential, KerrQI> energy_and_angularmomentum(kerr_init, m_dx, proca_field, m_p.center);

            //compute background first
            BoxLoops::loop(kerr_init, m_state_new, m_state_new, FILL_GHOST_CELLS);
            
            //compute diagnostics on each cell of current level
            BoxLoops::loop(
                energy_and_angularmomentum,
                m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
            );

            //excise within horizon
            BoxLoops::loop(
                ExcisionDiagnostics<KerrQI>(kerr_init, m_dx, m_p.center, m_p.inner_r),
                m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
                disable_simd()
            );
        }
    };

    //Calculate integrals of densities
    if ( m_level == min_level )
    {
        //Initialize the integrator
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);

        //integrate densities
        double rho_sum { amr_reductions.sum(c_rho) };
        double rhoJ_sum { amr_reductions.sum(c_rhoJ) };

        //setup integral file writer
        SmallDataIO integral_file(m_p.integrals_filename, m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step);

        //remove any duplicates
        integral_file.remove_duplicate_time_data();
        std::vector<double> data_for_writing = {rho_sum, rhoJ_sum};

        //save to disk
        if (first_step) {
            integral_file.write_header_line({"rho", "rhoJ"});
        }
        integral_file.write_time_data_line(data_for_writing);
    }

  /*   if (m_p.activate_extraction)
    {

        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {

            fillAllGhosts();
            ProcaPotential potential(m_p.potential_params);
            ProcaFieldWithPotential proca_field(potential, m_p.proca_params);

            //populate Weyl scalar values on grid
            
            BoxLoops::loop(
                MatterWeyl4<ProcaFieldWithPotential>(
                    proca_field,
                    m_p.extraction_params.center, 
                    m_dx, 
                    m_p.formulation, 
                    m_p.G_Newton),
                m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
            );

            //excise within horizon
            BoxLoops::loop(
                ExcisionDiagnostics(m_dx, m_p.center, 
                    0.0,m_p.extraction_params.extraction_radii[1]),
                m_state_diagnostics, 
                m_state_diagnostics,
                SKIP_GHOST_CELLS,
                disable_simd()
            );

            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                //refresh interpolator
                //fill ghosts manually
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt, m_time, first_step, m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    } */
/* 
    if (m_p.calculate_constraint_norms)
    {
        fillAllGhosts();
        BoxLoops::loop(
            Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
        );

        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double L2_Ham = amr_reductions.norm(c_Ham);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
            SmallDataIO constraint_file(m_p.data_path + "constraint_norms",
                                        m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step
            );
            constraint_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraint_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            }
            constraint_file.write_time_data_line({L2_Ham, L2_Mom});
        }
    } */

}