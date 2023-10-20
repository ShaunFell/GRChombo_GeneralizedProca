

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

//cell tagging
#include "FixedGridsTaggingCriterion.hpp"

//problem specific includes
#include "GammaCalculator.hpp"
#include "InitialProcaData.hpp"
#include "KerrBH.hpp"
#include "Potential.hpp"
#include "ProcaField.hpp"
#include "SetValue.hpp"
#include "Diagnostics.hpp"


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
            MatterConstraints<ProcaFieldWithPotential>(proca_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom1, c_Mom3)),
            Asquared,
            proca_constraint,
            proca_eff_met
            ),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
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

};

//tagging criterion for AMR
void ProcaFieldLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    //tag cells based on distance to boundary of computational domain
    BoxLoops::loop(
        FixedGridsTaggingCriterion(m_dx, m_level, 2.0*m_p.L, m_p.center),
        current_state, tagging_criterion
    );
};