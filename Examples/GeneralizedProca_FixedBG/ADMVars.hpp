

#ifndef ADMVARS_HPP
#define ADMVARS_HPP

#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

namespace ADMVars
{

template <class data_t>
struct Vars
{
    Tensor<2, data_t> h;         //conformal metric
    Tensor<2, data_t> A;         //conformal traceless part of  extrinsic curvature
    data_t chi;                          //conformal factor
    data_t K;                             // trace of physical extrinsic curvature
    data_t lapse;                      // lapse
    Tensor<1,data_t> shift;    // shift

    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        // Scalars
        define_enum_mapping(mapping_function, c_chi, chi);
        define_enum_mapping(mapping_function, c_K, K);
        define_enum_mapping(mapping_function, c_lapse, lapse);

        //vectors
        define_enum_mapping(mapping_function, GRInterval<c_shift1, c_shift3>(), shift);


        // Symmetric 2-tensors
        define_symmetric_enum_mapping(mapping_function, GRInterval<c_h11, c_h33>(), h);
        define_symmetric_enum_mapping(mapping_function, GRInterval<c_A11, c_A33>(), A);
    }


};

}

#endif //ADMVARS_HPP