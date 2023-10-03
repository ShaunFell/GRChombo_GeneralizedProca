#ifndef ADMVARIABLES_H_INCLUDED
#define ADMVARIABLES_H_INCLUDED

#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

namespace ADMVars
{

    template <class data_t> struct Vars
    {
        Tensor<2, data_t> gamma; 
        Tensor<2, data_t> K_tensor;
        data_t K;
        data_t lapse;
        Tensor<1, data_t> shift;
        Tensor<2, Tensor<1, data_t>> d1_gamma;
        Tensor<1, data_t> d1_lapse;
        Tensor<2, data_t> d1_shift;
    };
}//end of ADMVars namespace
















#endif //ADMVARIABLES_H_INCLUDED