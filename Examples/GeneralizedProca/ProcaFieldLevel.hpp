#ifndef PROCAFIELDLEVEL_H_INCLUDED
#define PROCAFIELDLEVEL_H_INCLUDED

#include "DefaultLevelFactory.hpp"
#include "AHInterpolation.hpp"
#include "GRAMRLevel.hpp"
#include "BHAMR.hpp"
#include "Potential.hpp"
#include "ProcaField.hpp"


class ProcaFieldLevel : public GRAMRLevel
{

    friend class DefaultLevelFactory<ProcaFieldLevel>;
    //inherit constructors from GRAMRLevel;
    using GRAMRLevel::GRAMRLevel;

    template <class SurfaceGeometry, class ApparentHorizon>
    using AHInterpolation = AHInterpolation_t<SurfaceGeometry, ApparentHorizon>;

    BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);

    //typedef the scalar field class;
    typedef ProcaField<ProcaPotential> ProcaFieldWithPotential;

    virtual void specificAdvance(); //do things at end of advance step, after RK4 calculation

    virtual void initialData(); //initialize data for the field

#ifdef CH_USE_HDF5
    virtual void prePlotLevel();
#endif

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                    const double a_time); //RHS routines used at each RK4 step
    
    virtual void specificUpdateODE(GRLevelData &a_soln, const GRLevelData &a_rhs,
                                    Real a_dt);

    virtual void preTagCells() override; //things to do before tagging cells (e.g. filling ghosts)

    virtual void computeTaggingCriterion(
        FArrayBox &tagging_criterion, const FArrayBox &current_state,
        const FArrayBox &current_state_diagnostics) override;

    virtual void specificPostTimeStep() override;
    
};

#endif //PROCAFIELDLEVEL_H_INCLUDED