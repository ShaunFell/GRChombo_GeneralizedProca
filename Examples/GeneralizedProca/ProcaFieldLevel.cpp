

//general includes common to most GR problems
#include "ProcaFieldLevel.hpp"
#include "AMRReductions.hpp"
#include "ComputePack.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"

//RHS Update
#include "MatterCCZ4RHS.hpp"

//constraint calculation
#include "NewMatterConstraints.hpp"

//Weyl extraction
#include "MatterWeyl4.hpp"
#include "WeylExtraction.hpp"
#include "SmallDataIO.hpp"


//cell tagging
#include "TaggingCriterion.hpp"

//problem specific includes
#include "GammaCalculator.hpp"
#include "InitialProcaData.hpp"
#include "KerrBH.hpp"
#include "Potential.hpp"
#include "ProcaField.hpp"
#include "SetValue.hpp"
#include "Diagnostics.hpp"
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"


//do things at end of advance step, after RK4 calculation
void ProcaFieldLevel::specificAdvance()
{
    //enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
        PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);


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

    BoxLoops::loop(
        make_compute_pack(
                        SetValue(0.), 
                        KerrBH(m_p.kerr_params, m_dx),
                        InitialProcaData(m_p.initialdata_params, m_p.potential_params, m_p.kerr_params, m_dx)
                        ),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);

    //excise inside horizon
    BoxLoops::loop(
        ExcisionEvolution<ProcaFieldWithPotential>(m_dx, m_p.center, m_p.inner_r),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS,disable_simd()
    );

    //check for nans in initial data
    if (m_p.nan_check){
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                        EXCLUDE_GHOST_CELLS, disable_simd());
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
            MatterConstraints<ProcaFieldWithPotential>(proca_field, 
                                                        m_dx, m_p.G_Newton, 
                                                        c_Ham, 
                                                        Interval(c_Mom1, c_Mom3),
                                                        c_Ham_abs_sum,
                                                        Interval(c_Mom_abs_sum1,c_Mom_abs_sum3)
                                                        ),
            Asquared,
            proca_constraint,
            proca_eff_met
            ),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
        );
    
    BoxLoops::loop(
        ExcisionDiagnostics(m_dx, m_p.center, m_p.inner_r, m_p.outer_r),
        m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
        disable_simd()
    );
};
#endif //CH_USE_HDF5

//RHS routines used at each RK4 step
void ProcaFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                const double a_time)
{
    //enforce trace free A_ij and positive conformal factor and lapse
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(), PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS
    );

    //Calculate MatterCCZ4 right hand side with matter_t = ProcaField
    ProcaPotential potential(m_p.potential_params);
    ProcaFieldWithPotential proca_field(potential, m_p.proca_params);
    if (m_p.max_spatial_derivative_order == 4){
        MatterCCZ4RHS<ProcaFieldWithPotential, MovingPunctureGauge, FourthOrderDerivatives> my_ccz4_matter(proca_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    } else if (m_p.max_spatial_derivative_order == 6 ) {
        MatterCCZ4RHS<ProcaFieldWithPotential, MovingPunctureGauge, SixthOrderDerivatives> my_ccz4_matter(proca_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }

    //excise inside horizon for evolution variables
    BoxLoops::loop(
        ExcisionEvolution<ProcaFieldWithPotential>(m_dx, m_p.center, m_p.inner_r),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS, disable_simd()
    );

};

void ProcaFieldLevel::specificUpdateODE(GRLevelData &a_soln, const GRLevelData &a_rhs,
                                Real a_dt)
{
    //Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
};


//things to do before tagging cells (e.g. filling ghosts)
void ProcaFieldLevel::preTagCells()
{
    //fill all ghosts
    fillAllGhosts();

    //setup class instances
    ProcaPotential potential(m_p.potential_params);
    ProcaFieldWithPotential proca_field(potential, m_p.proca_params);
    ProcaConstraint<ProcaPotential> proca_constraint(m_dx, m_p.potential_params.mass, m_p.proca_params.vector_damping, potential);

    //compute Hamiltonian and Guass diagnostics on each cell of current level as these are required for tagging
    BoxLoops::loop(
        make_compute_pack(
            MatterConstraints<ProcaFieldWithPotential>(proca_field, 
                                                        m_dx, m_p.G_Newton, 
                                                        c_Ham, 
                                                        Interval(c_Mom1, c_Mom3),
                                                        c_Ham_abs_sum,
                                                        Interval(c_Mom_abs_sum1,c_Mom_abs_sum3)
                                                        ),
            ProcaConstraint<ProcaPotential>(m_dx, 
                                            m_p.potential_params.mass, 
                                            m_p.proca_params.vector_damping, 
                                            potential)
            ),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
        );
    
    //excise diagnostics according to parameters set in parameter file
    BoxLoops::loop(
        ExcisionDiagnostics(m_dx, m_p.center, m_p.inner_r, m_p.outer_r),
        m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
        disable_simd()
    );

};


//compute tagging criteria for grid
void ProcaFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                             const FArrayBox &current_state,
                                             const FArrayBox &current_state_diagnostics)
{
    //tag cells based on Hamiltonian constraint
    BoxLoops::loop(
        CustomTaggingCriterion(
                                m_dx, m_level, 2.0*m_p.L, 
                                m_p.extraction_params, 
                                m_p.max_level-1,
                                m_p.activate_extraction
                            ),
        current_state_diagnostics, 
        tagging_criterion
    );
}


void ProcaFieldLevel::specificPostTimeStep()
{
    CH_TIME("ProcaFieldLevel::specificPostTimeStep");
    bool first_step = (m_time == m_dt); //is this the first call of posttimestep?

    if (m_p.activate_extraction == 1)
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
    }

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
    }

}
