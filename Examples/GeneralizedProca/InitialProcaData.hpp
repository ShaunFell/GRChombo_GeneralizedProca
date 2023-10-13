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
        double width;
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
        
        
        //flush all variables on cell
        MatterVars<data_t> mattervars;
        VarsTools::assign(mattervars,0.);


        //compute Kerr metric components
        Tensor<2,data_t> spherical_g;
        Tensor<2,data_t> spherical_K;
        Tensor<1,data_t> spherical_shift;
        Tensor<1,data_t> xyz;
        data_t kerr_lapse;

        //populate tensors and calculate rotation matrices
        const Tensor<1, double> spin_direction { m_paramsKerr.spin_direction[0], m_paramsKerr.spin_direction[1], m_paramsKerr.spin_direction[2] };
        static const Tensor<1,double> z_dir = {0.,0.,1.};
        Tensor<2,double> R = CoordinateTransformations::rotation_matrix(spin_direction, z_dir); //rotation matrix from standard cartesian to spin-aligned coordinates
        Tensor<1,data_t> coord_location {coords.x, coords.y, coords.z};
        xyz = CoordinateTransformations::transform_vector(coord_location, R );


        compute_kerr(spherical_g, spherical_K, spherical_shift, kerr_lapse, xyz); //compute kerr metric in spherical coords
        Tensor<2, data_t> cartesian_h = CoordinateTransformations::spherical_to_cartesian_LL(spherical_g, xyz[0], xyz[1], xyz[2]); //transform to cartesian coords
        Tensor<2, data_t> gamma_LL = CoordinateTransformations::transform_tensor_LL(cartesian_h, R); //rotate back to original coords

        const data_t det_gamma { TensorAlgebra::compute_determinant_sym(gamma_LL) };
        const data_t conformalFact { pow(det_gamma, -1./3.) };

        const data_t kerrMass = m_paramsKerr.mass;
        const data_t kerrSpin = m_paramsKerr.spin;
        const data_t kerrSpin2 = kerrSpin*kerrSpin;
        const data_t coordZ = coords.z;
        const data_t rho = coords.get_radius(); //x^2 + y^2 + z^2
        const data_t rho2 = rho*rho; //r^2

        const data_t radial2 = 0.5*(rho2 - kerrSpin2)  + sqrt(0.25*(rho2-kerrSpin2) + kerrSpin2*coordZ*coordZ);
        const data_t radius = sqrt(radial2);

        data_t alpha = kerrMass*m_paramsPotential.mass;
        data_t r0 = 1.0/(m_paramsPotential.mass*alpha); //peak of boson condensate

        
        
        mattervars.Avec[0] = m_params.amplitude*pow(conformalFact, 3./2.)*exp(-radius/r0);
        mattervars.Avec[1] = 0.;
        mattervars.Avec[2] = 0.;
        mattervars.phi = 0.;
        mattervars.Z = 0.;
        FOR1(i){
            mattervars.Evec[i] = 0.;
        }
    
        //################################################################################
#ifdef EQUATION_DEBUG_MODE
        DEBUG_HEADER;
        DEBUG_OUT2(mattervars.Z, mattervars.phi);
        DEBUG_OUT3(mattervars.Avec[0], mattervars.Avec[2], mattervars.Avec[1]);
        DEBUG_OUT3(mattervars.Evec[0], mattervars.Evec[2], mattervars.Evec[1]);
        DEBUG_END;
#endif //EQUATION_DEBUG_MODE
        //################################################################################

        current_cell.store_vars(mattervars);

    };

};





#endif //INITIALPROCADATA_H_INCLUDED