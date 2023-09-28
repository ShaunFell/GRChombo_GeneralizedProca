#include "ProcaField.hpp"


template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars, //the value of the variables
        const MetricVars<data_t> &metric_vars, //value of metric variables
        const vars_t<Tensor<1,data_t>> &d1, //the 1st derivatives
        const Tensor<2, data_t> &gamma_UU, //the inverse metric
        const Tensor<3, data_t> &chris_phys_ULL //conformal christoffel symbols
    ) const {
        emtensor_t<data_t> out;

        data_t rho_potential { 0 };
        Tensor<1, data_t> Si_potential;
        Tensor<2, data_t> Sij_potential;
        m_
    }