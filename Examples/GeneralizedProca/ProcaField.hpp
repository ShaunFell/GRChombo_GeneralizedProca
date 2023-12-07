#ifndef PROCAFIELD_H_INCLUDED
#define PROCAFIELD_H_INCLUDED

#include "CCZ4Geometry.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Potential.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the Proca fieldspecific elements such as the EMTensor 
//  and matter evolution

template <class potential_t>
class ProcaField
{
public:
    struct params_t{
        double vector_damping; // vector damping coefficient
    };

    const params_t m_params;
    const potential_t m_potential;

    //! Constructor, inputs are matter params
    ProcaField(const potential_t potential, params_t a_params) : 
        m_params{a_params}, m_potential{potential}
    {}

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        data_t phi;             // scalar part
        Tensor<1, data_t> Avec; // vector (spatial) part, indices DOWN
        Tensor<1, data_t> Evec; // electric part of the strength tensor, indices UP
        data_t Z;               // auxiliary damping scalar

        // function that maps between above Vars and Chombo grid variables
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function){

            using namespace VarsTools;
            define_enum_mapping(mapping_function, c_phi, phi);
            define_enum_mapping(mapping_function, c_Z, Z);
            define_enum_mapping(mapping_function, GRInterval<c_Avec1, c_Avec3>(), Avec);
            define_enum_mapping(mapping_function, GRInterval<c_Evec1, c_Evec3>(), Evec);

        }
    };

    //! Structure holding the matter field variables that require 2nd derivatives
    template <class data_t> struct Diff2Vars 
    {
        Tensor<1, data_t> Avec;

        // function that maps between above Vars and Chombo grid variables
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function){
            using namespace VarsTools;
            define_enum_mapping(mapping_function, GRInterval<c_Avec1, c_Avec3>(), Avec);
        }
    };

    //! Method that computes EM tensor, given vars and derivatives
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars, // the value of the variables
        const vars_t<Tensor<1,data_t>> &d1, // the 1st derivatives
        const Tensor<2, data_t> &h_UU, // the inverse metric with indices UP
        const Tensor<3, data_t> &chris_ULL // conformal christoffel symbols
    ) const;


    //! Method which adds in the matter field RHS, given vars and derivatives
    template <class data_t, 
        template <typename> class vars_t, 
        template <typename> class diff2_vars_t, 
        template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs, // RHS terms for all vars
        const vars_t<data_t> &vars, // the value fo the variables
        const vars_t<Tensor<1, data_t>> &d1, // value of 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, // 2nd derivs
        const vars_t<data_t> &advec // value of the beta^i d_i(var) terms
    ) const;

};

#include "ProcaField.impl.hpp"
#endif //PROCAFIELD_H_INCLUDED

