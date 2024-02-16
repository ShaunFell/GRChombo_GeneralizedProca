#ifndef INITIALPROCADATA_H_INCLUDED
#define INITIALPROCADATA_H_INCLUDED

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "CoordinateTransformations.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ProcaField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"
#include "Potential.hpp"
#include "KerrBH.hpp"


#ifdef EQUATION_DEBUG_MODE
#include <cstdlib>
#endif

class InitialProcaData: public KerrBH
{
public:

    //typename all variables, CCZ4 variables + Matter
    template <class data_t>
    using Vars = typename MatterCCZ4RHS<ProcaField<ProcaPotential>>::template Vars<data_t>;
    
    //typename only matter variables
    template <class data_t>
    using MatterVars = typename ProcaField<ProcaPotential>::template Vars<data_t>;

    struct init_params_t
    {
        double amplitude;
    };

    using PotentialParams = ProcaPotential::params_t;
    using KerrParams = KerrBH::params_t;

protected:
    double m_dx;
    const init_params_t m_params;
    const PotentialParams m_paramsPotential;
    const KerrParams m_paramsKerr;

public:
    
    //constructor
    InitialProcaData(init_params_t a_params, PotentialParams b_params, KerrParams c_params, double a_dx): 
        KerrBH(c_params, a_dx), m_dx{a_dx}, m_params{a_params}, m_paramsPotential{b_params}, m_paramsKerr{c_params}
    {
    };

    template <class data_t>
    void compute(Cell<data_t> current_cell) const
    {
        //based off the initial conditions used in    http://arxiv.org/abs/1705.01544

        //location of cell
        Coordinates<data_t> coords(current_cell, m_dx, m_paramsKerr.center);
        
        //load variables from already calculated kerr black hole (See ProcaFieldLevel.cpp::59)
        const auto vars = current_cell.template load_vars<Vars>();
        
        //flush all variables on cell
        MatterVars<data_t> mattervars;
        VarsTools::assign(mattervars,0.);
        
        const data_t kerrMass = m_paramsKerr.mass;
        const data_t kerrSpin = m_paramsKerr.spin;
        const data_t kerrSpin2 = kerrSpin*kerrSpin;      
        const data_t rP_BL = kerrMass * (1 + sqrt(1 - kerrSpin2));
        const data_t rho = coords.get_radius(); //x^2 + y^2 + z^2


        //Use relation for quasi-isotropic coords to boyer-lindquist
        const data_t r_BL = rho * ( 1 + rP_BL/(4 * rho )) *  ( 1 + rP_BL/(4 * rho ));


        data_t alpha = kerrMass*m_paramsPotential.mass;
        data_t r0_BL { 1.0/(m_paramsPotential.mass*alpha) };

        // if outside horizon, set initial data, else if inside, truncate to 0
        mattervars.Avec[0] = m_params.amplitude*pow(vars.chi, 3.)*exp(-r_BL/r0_BL);
        mattervars.Avec[1] = 0.;
        mattervars.Avec[2] = 0.;
        mattervars.phi = 0.;
        mattervars.Z = 0.;
        FOR1(i){
            mattervars.Evec[i] = 0.;
        };

        current_cell.store_vars(mattervars);

    };

};





#endif //INITIALPROCADATA_H_INCLUDED