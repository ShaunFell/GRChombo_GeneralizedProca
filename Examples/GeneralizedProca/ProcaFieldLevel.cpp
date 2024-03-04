

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

//Chi Relaxation
#include "ChiRelaxation.hpp"


//cell tagging
#include "TaggingCriterion.hpp"
#include "FixedGridsTaggingCriterion.hpp"

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

#include "AHInitialGuess.hpp"


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

    #ifdef USE_AHFINDER
    //apparently this is needed for the AHFinder
    if (m_p.AH_activate) {
        BoxLoops::loop(
            Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
        );
    };
    #endif //USE_AHFINDER

};

#ifdef CH_USE_HDF5
//things to do before outputting a checkpoint file
void ProcaFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    ProcaPotential potential(m_p.potential_params);
    ProcaFieldWithPotential proca_field(potential, m_p.proca_params);
    GaussConstraint<ProcaPotential> proca_constraint(m_dx, m_p.potential_params.mass, m_p.proca_params.vector_damping, potential);
    SecondClassConstraint<ProcaPotential> proca_eff_met(m_dx, m_p.potential_params.mass, m_p.proca_params.vector_damping, potential);
    ProcaSquared Asquared(m_dx);
    EnergyAndAngularMomentum<ProcaFieldWithPotential> EM(m_dx, proca_field, m_p.center);
    MatterConstraints<ProcaFieldWithPotential> matter_constraints(proca_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom1, c_Mom3), c_Ham_abs_sum,Interval(c_Mom_abs_sum1,c_Mom_abs_sum3));

    //compute diagnostics on each cell of current level
    BoxLoops::loop(
        make_compute_pack(
            matter_constraints,
            Asquared,
            proca_constraint,
            proca_eff_met,
            EM
            ),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
        );
    

    #ifdef USE_AHFINDER
    //excision of diagnostics
    if (m_p.excise_with_AH && m_p.AH_activate)
    {
        //shouldnt need to resolve AH

        //query AH 
        auto AH { m_bh_amr.m_ah_finder.get(0) };
        double minRadius { AH -> get_min_F() };

        ExcisionDiagnosticsWithAH<ProcaFieldWithPotential> excisor (m_dx, m_p.center, minRadius, m_p.AH_buffer);

        BoxLoops::loop(
            excisor,
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd()
        ); 

    } else if (m_p.excise_with_chi && m_p.AH_activate) 
    {
        
        //Excise with conformal factor
        auto AH_mass { m_bh_amr.m_ah_finder.get(0) -> m_mass };
        auto AH_spin { m_bh_amr.m_ah_finder.get(0) -> m_spin };
        double AH_dimless_spin { AH_spin / AH_mass };

        ExcisionDiagnosticsWithChi<ProcaFieldWithPotential> excision_init(m_dx, m_p.kerr_params.center, AH_dimless_spin, m_p.outer_excision);

        BoxLoops::loop(
            excision_init,
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS, disable_simd()
        ) ;

    }
    #endif //USE_AHFINDER

    if (m_p.excise_with_cutoff)
    {
    
        ExcisionDiagnostics<ProcaFieldWithPotential> excisor (m_dx, m_p.center, m_p.inner_r);
        //excise diagnostics according to parameters set in parameter file
        BoxLoops::loop(
            excisor,
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd()
        ); 

    }// end of excision

    #ifdef USE_AHFINDER
    //already calculated in specific PostTimeStep
    if(m_p.AH_activate && m_bh_amr.m_ah_finder.need_diagnostics(m_dt, m_time)){
        return;
    }
    #endif //USE_AHFINDER
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

    //Moving puncture gauge to handle spacetime singularites
    ProcaPotential potential(m_p.potential_params);
    ProcaFieldWithPotential proca_field(potential, m_p.proca_params);

    if (a_time > m_p.relaxation_time) 
    {
        //carry out normal CCZ4 evolution
        if (m_p.max_spatial_derivative_order == 4){
            MatterCCZ4RHS<ProcaFieldWithPotential, MovingPunctureGauge, FourthOrderDerivatives> my_ccz4_matter(proca_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation, m_p.G_Newton);
            BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
        } else if (m_p.max_spatial_derivative_order == 6 ) {
            MatterCCZ4RHS<ProcaFieldWithPotential, MovingPunctureGauge, SixthOrderDerivatives> my_ccz4_matter(proca_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation, m_p.G_Newton);
            BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
        }
    } else {
        //perform chi relaxation
        ChiRelaxation<ProcaFieldWithPotential> chi_relaxation (proca_field,  m_dx, m_p.relaxation_speed, m_p.G_Newton );
        BoxLoops::loop(
            chi_relaxation,
            a_soln, a_rhs, EXCLUDE_GHOST_CELLS
        );
    } //end of relaxation if-statement

    
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
    CH_TIME("ProcaFieldLevel::preTagCells");

    if (m_p.activate_gauss_tagging || m_p.activate_ham_tagging || m_p.activate_effmetric_tagging)
    {

        //fill all ghosts
        fillAllGhosts();

        //setup class instances
        ProcaPotential potential(m_p.potential_params);
        ProcaFieldWithPotential proca_field(potential, m_p.proca_params);
        GaussConstraint<ProcaPotential> proca_constraint(m_dx, m_p.potential_params.mass, m_p.proca_params.vector_damping, potential);
        SecondClassConstraint<ProcaPotential> proca_eff_met(m_dx, m_p.potential_params.mass, m_p.proca_params.vector_damping, potential);
        MatterConstraints<ProcaFieldWithPotential> matter_constraints(proca_field, m_dx, m_p.G_Newton, 
                                                            c_Ham,  Interval(c_Mom1, c_Mom3),
                                                            c_Ham_abs_sum, Interval(c_Mom_abs_sum1,c_Mom_abs_sum3)
                                                            );


        //compute Hamiltonian and Guass diagnostics on each cell of current level as these are required for tagging
        BoxLoops::loop(
            make_compute_pack(proca_constraint, proca_eff_met, matter_constraints), 
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
        );
        
        if (m_time==0.)
        { //run excision at t==0 using knowledge of horizon from initial data
            
            double kerr_spin = m_p.kerr_params.spin;
            double kerr_mass = m_p.kerr_params.mass;
            double cutoff { 1./4. * kerr_mass * (1 + sqrt(1 - kerr_spin*kerr_spin)) };
            //excise diagnostics according to parameters set in parameter file
            ExcisionDiagnostics<ProcaFieldWithPotential> excisor (m_dx, m_p.center, cutoff);
            BoxLoops::loop(
                excisor,
                m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS, disable_simd()
            ) ;
        } else
        {
            //excision of diagnostics
            #ifdef USE_AHFINDER
            if (m_p.excise_with_AH && m_p.AH_activate)
            {
                //query AH 
                auto AH { m_bh_amr.m_ah_finder.get(0) };
                double minRadius { AH -> get_min_F() };

                ExcisionDiagnosticsWithAH<ProcaFieldWithPotential> excisor (m_dx, m_p.center, minRadius, m_p.AH_buffer);
                BoxLoops::loop(
                    excisor,
                    m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS, disable_simd()
                ) ;

            } else if (m_p.excise_with_chi && m_p.AH_activate) 
            {
                //Excise with conformal factor
                auto AH_mass { m_bh_amr.m_ah_finder.get(0) -> m_mass };
                auto AH_spin { m_bh_amr.m_ah_finder.get(0) -> m_spin };
                double AH_dimless_spin { AH_spin / AH_mass };

                ExcisionDiagnosticsWithChi<ProcaFieldWithPotential> excisor(m_dx, m_p.kerr_params.center, AH_dimless_spin, m_p.outer_excision);

                BoxLoops::loop(
                    excisor,
                    m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS, disable_simd()
                ) ;

            } 
            #endif //USE_AHFINDER

            if (m_p.excise_with_cutoff)
            {
                //excise diagnostics according to parameters set in parameter file
                ExcisionDiagnostics<ProcaFieldWithPotential> excisor (m_dx, m_p.center, m_p.inner_r);
                BoxLoops::loop(
                    excisor,
                    m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS, disable_simd()
                ) ;

            }// end of excision
        }

    };

};


//compute tagging criteria for grid
void ProcaFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                             const FArrayBox &current_state,
                                             const FArrayBox &current_state_diagnostics)
{
    CH_TIME("ProcaFieldLevel::computeTaggingCriterion");

    if (m_p.activate_gauss_tagging || m_p.activate_ham_tagging || m_p.activate_extraction) {
        CustomTaggingCriterion tagger(
                                    m_dx, m_level, m_p.grid_scaling*m_p.L, 
                                    m_p.center,
                                    m_p.extraction_params, 
                                    m_p.activate_extraction,
                                    m_p.activate_gauss_tagging,
                                    m_p.activate_ham_tagging
                                );

        BoxLoops::loop(tagger, current_state_diagnostics, tagging_criterion, disable_simd());
    } else {
       FixedGridsTaggingCriterion tagger(m_dx, m_level,
                                                    m_p.grid_scaling * m_p.L, m_p.center);
        BoxLoops::loop(tagger, current_state, tagging_criterion, disable_simd());

    };

    
    
}


void ProcaFieldLevel::specificPostTimeStep()
{
    CH_TIME("ProcaFieldLevel::specificPostTimeStep");
    bool first_step = (m_time == 0.); //is this the first call of posttimestep?


    //  ##### AH Finder ####
    #ifdef USE_AHFINDER
    if (m_bh_amr.m_ah_finder.need_diagnostics(m_dt, m_time) && m_p.AH_activate)
    {
        fillAllGhosts();
        BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }

    if (m_p.AH_activate)
    {
        if (m_level == m_p.AH_params.level_to_run)
        {
            //Formulate smart initial guess using previous solution
            auto new_guess { m_bh_amr.m_ah_finder.get(0) -> get_ave_F() }; //average radius of current AH

            //set new guess
            auto initguess { m_bh_amr.m_ah_finder.get(0) -> get_petsc_solver().get_initial_guess() };
            dynamic_cast<AHInitialGuessConstant&>(*initguess).m_initial_guess = new_guess;

            pout() << "Solving AH with initial guess = " << initguess -> get(0,0) << endl;
            m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);    
        }
        
        if (first_step)
        {
            pout() << "Running AH Finder at t=0" << endl;
            m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);    
        }
    }
    #endif //USE_AHFINDER


    
    //  ##### Waveform Extraction ####


    if (m_p.activate_extraction)
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

            //excision of diagnostics
            #ifdef USE_AHFINDER
            if (m_p.excise_with_AH && m_p.AH_activate)
            {
                //shouldnt need to resolve AH

                //query AH 
                auto AH { m_bh_amr.m_ah_finder.get(0) };
                double minRadius { AH -> get_min_F() };

                ExcisionDiagnosticsWithAH<ProcaFieldWithPotential> excisor (m_dx, m_p.center, minRadius, m_p.AH_buffer);

                BoxLoops::loop(
                    excisor,
                    m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
                    disable_simd()
                ); 

            } else if (m_p.excise_with_chi && m_p.AH_activate) 
            {
                
                //Excise with conformal factor
                auto AH_mass { m_bh_amr.m_ah_finder.get(0) -> m_mass };
                auto AH_spin { m_bh_amr.m_ah_finder.get(0) -> m_spin };
                double AH_dimless_spin { AH_spin / AH_mass };

                ExcisionDiagnosticsWithChi<ProcaFieldWithPotential> excision_init(m_dx, m_p.kerr_params.center, AH_dimless_spin, m_p.outer_excision);

                BoxLoops::loop(
                    excision_init,
                    m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS, disable_simd()
                ) ;

            } 
            #endif //USE_AHFINDER

            if (m_p.excise_with_cutoff)
            {
            
                ExcisionDiagnostics<ProcaFieldWithPotential> excisor (m_dx, m_p.center, m_p.inner_r);
                //excise diagnostics according to parameters set in parameter file
                BoxLoops::loop(
                    excisor,
                    m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
                    disable_simd()
                ); 

            }// end of excision

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



    //  ##### Constraint Norms ####
    int coarsest_level = 0;
    bool at_course_timestep_on_any_level = at_level_timestep_multiple(coarsest_level);
    // Dont want to do this on every level at posttimestep, only where we will actually
    // do the outputs, so on the coarsest level timestep (but still need to fill data
    // on the finer levels)

    if (m_p.calculate_norms && at_course_timestep_on_any_level)
    {
        // this operation to fill all ghosts is expensive so avoid doing it too
        // often is possible. We could specific only some variables get their
        // ghosts filled (as in the Weyl case above), but for the constraints we probably
        // need the majority so we won't gain much by specifying. But at least avoid doing it on
        // every substep.
        fillAllGhosts();

        ProcaPotential potential(m_p.potential_params);
        ProcaFieldWithPotential proca_field(potential, m_p.proca_params);
        EnergyAndAngularMomentum<ProcaFieldWithPotential> EM(m_dx, proca_field, m_p.center);
        BoxLoops::loop(
            make_compute_pack(
                Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
                EM
	        ),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
        );

        // Once values calculated on all levels, only the coarsest level does
        // the integral and the output. It will use finer level data if it exists.
        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double L2_Ham = amr_reductions.norm(c_Ham);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
            double SUM_rho = amr_reductions.sum(c_rho);
            double SUM_rhoE = amr_reductions.sum(c_rhoE);
            double SUM_rhoJ = amr_reductions.sum(c_rhoJ);

            SmallDataIO constraint_file(m_p.data_path + "constraint_norms",
                                        m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step
            );
            constraint_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraint_file.write_header_line({"L^2_Ham", "L^2_Mom", "SUM_rho", "SUM_rhoE", "SUM_rhoJ"});
            }
            constraint_file.write_time_data_line({L2_Ham, L2_Mom, SUM_rho, SUM_rhoE, SUM_rhoJ});
        }
    }

}
