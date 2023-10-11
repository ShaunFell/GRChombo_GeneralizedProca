#ifndef PROCAFIELD_H_INCLUDED
#define PROCAFIELD_H_INCLUDED

#include "ADMVariables.hpp"
#include "CCZ4Geometry.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Potential.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"



template <class potential_t>
class ProcaField
{
public: 
    struct params_t{
        double vector_damping; //local copy of vector damping coefficient
    };

    const params_t m_params;
    const potential_t m_potential;

    //constructor, inputs are matter params
    ProcaField(const potential_t potential, params_t a_params): m_params{a_params}, m_potential{potential}{};

    template <class data_t> 
    struct Vars{
        data_t phi;
        data_t Z; //auxilliary damping scalar
        Tensor<1, data_t> Avec; //Spatial part of Proca field
        Tensor<1, data_t> Evec; //Electric part of Proca field strength tensor

        //provide function that maps between above Vars and Chombo grid variables
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function){
            using namespace VarsTools; //define_enum_mapping is part of VarsTools namespace
            define_enum_mapping(mapping_function, c_phi, phi);
            define_enum_mapping(mapping_function, c_Z, Z);
            define_enum_mapping(mapping_function, GRInterval<c_Avec1, c_Avec3>(), Avec);
            define_enum_mapping(mapping_function, GRInterval<c_Evec1, c_Evec3>(), Evec);
        }
    }; //end of struct Vars

    //structure holding the matter field variables that require 2nd derivatives
    template <class data_t>
    struct Diff2Vars {
        Tensor<1, data_t> Avec;

        //provide function that maps between above Vars and Chombo grid variables
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function){
            using namespace VarsTools; //define_enum_mapping is part of VarsTools namespace
            define_enum_mapping(mapping_function, GRInterval<c_Avec1, c_Avec3>(), Avec);
        }
    }; //end of struct Diff2Vars

    //struct for non-grid ADM vars
    template <class data_t>
    using MetricVars = typename ADMVars::template Vars<data_t>;

    //method that computes EM tensor, given vars and derivatives
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars, //the value of the variables
        const vars_t<Tensor<1,data_t>> &d1, //the 1st derivatives
        const Tensor<2, data_t> &gamma_UU, //the inverse metric
        const Tensor<3, data_t> &chris_phys_ULL //conformal christoffel symbols
    ) const;


    //method which adds in the matter field RHS, given vars and derivatives
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t, template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs, //RHS terms for all vars
        const vars_t<data_t> &vars, //the value fo the variables
        const vars_t<Tensor<1, data_t>> &d1, //value of 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //2nd derivs
        const vars_t<data_t> &advec //value of the beta^i d_i(var) terms
    ) const; 
};

#endif //PROCAFIELD_H_INCLUDED

