#ifndef INITIALPROCADATA_H_INCLUDED
#define INITIALPROCADATA_H_INCLUDED

#include "ADMVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
/* #include "CoordinateTransformations.hpp" */
#include "ProcaField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"
#include "Potential.hpp"
#include "KerrSchild.hpp"

class InitialProcaData: public KerrSchild
{
public:

    //typename all variables, ADM variables + Matter
    template <class data_t>
    using MetricVars = ADMVars::template Vars<data_t>;
    
    //typename only matter variables
    template <class data_t>
    using MatterVars = typename ProcaField<ProcaPotential>::template Vars<data_t>;

    struct init_params_t
    {
        double amplitude;
    };

    using PotentialParams = ProcaPotential::params_t;
    using KerrParams = KerrSchild::params_t;

protected:
    double m_dx;
    const init_params_t m_params;
    const PotentialParams m_paramsPotential;
    const KerrParams m_paramsKerr;

public:
    
    //constructor
    InitialProcaData(init_params_t a_params, PotentialParams b_params, KerrParams c_params, double a_dx): 
        KerrSchild(c_params, a_dx), m_dx{a_dx}, m_params{a_params}, m_paramsPotential{b_params}, m_paramsKerr{c_params}
    {
    };

    template <class data_t>
    void compute(Cell<data_t> current_cell) const
    {
        //based off the initial conditions used in    http://arxiv.org/abs/1705.01544

        //location of cell
        Coordinates<data_t> coords(current_cell, m_dx, m_paramsKerr.center);
        
        //load kerr variables
        const auto metric_vars = current_cell.template load_vars<MetricVars>();
        
        //flush all variables on cell
        MatterVars<data_t> mattervars;
        VarsTools::assign(mattervars,0.);
        
        const data_t kerrMass = m_paramsKerr.mass;
        const data_t kerrSpin = m_paramsKerr.spin;
        const data_t kerrSpin2 = kerrSpin*kerrSpin;      
        const data_t rP_QI = 1./4. * (kerrMass + sqrt(kerrMass*kerrMass - kerrSpin2));
        const data_t rP_BL = 4. * rP_QI; 
        const data_t rho = coords.get_radius(); //x^2 + y^2 + z^2

        //convert the quasi-isotropic radial coordinate to boyer-lindquist
        auto QI_to_BL {
            [rP_BL] (data_t r) {
                return r * (1 + rP_BL/(4*r)) * (1 + rP_BL/(4*r));
            }
        };
        const data_t r_BL { QI_to_BL(rho) };
        const data_t sinTheta_BL { sqrt(coords.x*coords.x + coords.y*coords.y) / rho };
        const data_t cosTheta_BL { coords.z / rho };
/*         const data_t sinPhi_BL { coords.y / sqrt(coords.x * coords.x + coords.y * coords.y) };
        const data_t cosPhi_BL { coords.x / sqrt(coords.x * coords.x + coords.y * coords.y) };
        const data_t x_BL { r_BL * cosPhi_BL * sinTheta_BL };
        const data_t y_BL { r_BL * sinPhi_BL * sinTheta_BL };
        const data_t z_BL { r_BL * cosTheta_BL }; */

        data_t detGamma_BL { ( pow( (r_BL*r_BL + kerrSpin2*cosTheta_BL*cosTheta_BL) ,2. ) * sinTheta_BL*sinTheta_BL * 
                                ( r_BL*r_BL + kerrSpin2 + ( 2*r_BL*kerrSpin2*kerrMass*sinTheta_BL*sinTheta_BL ) / (r_BL*r_BL + kerrSpin2*cosTheta_BL*cosTheta_BL) )
                            ) / (r_BL*r_BL + kerrSpin2 - 2*r_BL*kerrMass)
         };
        data_t conformalFactor_BL { pow(detGamma_BL, -1.0/3.0) };



        data_t alpha = kerrMass*m_paramsPotential.mass;
        data_t r0_BL { 1.0/(m_paramsPotential.mass*alpha) };


        mattervars.Avec[0] = m_params.amplitude*pow(metric_vars.chi, 3./2.)*exp(-rho/r0_BL);
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