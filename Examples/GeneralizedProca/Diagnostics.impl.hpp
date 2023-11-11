

#if !defined(DIAGNOSTIC_H_INCLUDED)
#error "This file should only be included through Diagnostics.hpp"
#endif

#ifndef DIAGNOSTIC_IMPL_H_INCLUDED
#define DIAGNOSTIC_IMPL_H_INCLUDED


//Proca constraint constructor
template <class potential_t>
ProcaConstraint<potential_t>::ProcaConstraint(double dx, double a_vector_mass, double a_vector_damping, const potential_t potential): 
m_deriv(dx), m_vector_mass(a_vector_mass), m_vector_damping(a_vector_damping), m_potential(potential)
{
};

//Proca constraint calculator
template <class potential_t>
template <class data_t>
data_t ProcaConstraint<potential_t>::constraint_equations(Cell<data_t> current_cell) const
{
    //load variables from Chombo grid
    const auto vars = current_cell.template load_vars<Vars>();
    const auto varsd1 = m_deriv.template diff1<Vars>(current_cell);

    //compute contravariant conformal spatial metric
    const Tensor<2, data_t> h_UU { TensorAlgebra::compute_inverse_sym(vars.h) };
    
    //compute contravariant and covariant physical spatial metric 
    Tensor<2, data_t> gamma_LL;
    FOR2(i,j){
        gamma_LL[i][j] = vars.h[i][j]/vars.chi;
    };
    const Tensor<2,data_t> gamma_UU { TensorAlgebra::compute_inverse_sym(gamma_LL) };

    //compute conformal christoffel symbols
    const Tensor<3, data_t> chris_ULL { TensorAlgebra::compute_christoffel(varsd1.h, h_UU).ULL };

    //compute physical christoffel symbols
    const Tensor<3, data_t> chris_phys { TensorAlgebra::compute_phys_chris(varsd1.chi, vars.chi, vars.h, h_UU, chris_ULL) }; 
    
    //the return value
    data_t gauss_constraint;

    //compute potential
    data_t V { 0. };
    data_t dVdA { 0. };
    data_t dVddA { 0. };
    m_potential.compute_potential(V, dVdA, dVddA, vars, gamma_UU);

    gauss_constraint = 2*dVdA*vars.phi;

    FOR1(i)
    {
        gauss_constraint += varsd1.Evec[i][i];

        FOR1(j)
        {
            gauss_constraint += chris_phys[i][i][j]*vars.Evec[j];
        };
    };

    return gauss_constraint;
};

template <class potential_t>
template<class data_t>
void ProcaConstraint<potential_t>::compute(Cell<data_t> current_cell) const
{
    data_t gauss_constraint { constraint_equations(current_cell) };

    current_cell.store_vars(gauss_constraint,c_gauss);
};

template <class potential_t>
template<class data_t>
void EffectiveMetric<potential_t>::compute(Cell<data_t> current_cell) const
{
    //load variables from Chombo grid
    const auto vars = current_cell.template load_vars<Vars>();
    const auto varsd1 = m_deriv.template diff1<Vars>(current_cell);

    //compute contravariant conformal spatial metric
    const Tensor<2, data_t> h_UU { TensorAlgebra::compute_inverse_sym(vars.h) };
    
    //compute contravariant and covariant physical spatial metric 
    Tensor<2, data_t> gamma_LL;
    FOR2(i,j){
        gamma_LL[i][j] = vars.h[i][j]/vars.chi;
    };
    const Tensor<2,data_t> gamma_UU { TensorAlgebra::compute_inverse_sym(gamma_LL) };
    
    //the return value
    data_t gnn;

    //compute potential
    data_t V { 0. };
    data_t dVdA { 0. };
    data_t dVddA { 0. };
    m_potential.compute_potential(V, dVdA, dVddA, vars, gamma_UU);

    gnn = dVdA - 2*dVddA*vars.phi*vars.phi;



    current_cell.store_vars(gnn, c_gnn);



}

template <class data_t>
void ProcaSquared::compute(Cell<data_t> current_cell) const
{
    //load variables from Chombo grid
    const auto vars = current_cell.template load_vars<Vars>();
    const auto varsd1 = m_deriv.template diff1<Vars>(current_cell);

    //compute contravariant conformal spatial metric
    const auto h_UU { TensorAlgebra::compute_inverse_sym(vars.h) };
    
    //compute contravariant and covariant physical spatial metric 
    Tensor<2, data_t> gamma_LL;
    FOR2(i,j){
        gamma_LL[i][j] = vars.h[i][j]/vars.chi;
    };
    const Tensor<2,data_t> gamma_UU { TensorAlgebra::compute_inverse_sym(gamma_LL) };

    data_t Asquared;
    Asquared = -vars.phi*vars.phi;
    FOR2(i,j)
    {
        Asquared += gamma_UU[i][j] *vars.Avec[i]*vars.Avec[j];
    };

    current_cell.store_vars(Asquared, c_Asquared);
};

template <class matter_t>
FluxDensities<matter_t>::FluxDensities(double a_vector_mass, double a_vector_damping, double a_dx, std::array<double, CH_SPACEDIM> a_center, const matter_t a_matter):
m_deriv{a_dx}, m_dx{a_dx}, m_vector_damping{a_vector_damping}, m_vector_mass{a_vector_mass}, m_matter{a_matter}, m_center{a_center}
{
};

template <class matter_t>
template <class data_t>
void FluxDensities<matter_t>::compute(Cell<data_t> current_cell) const
{
    //extract data from grid
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    Coordinates<data_t> coords(current_cell, m_dx, m_center);

    //extract determinant of spatial metric
    const data_t det_gamma { pow(vars.chi, -3./2.) };

    //compute inverse conformal metric and conformal christoffel symbols
    const Tensor<2, data_t> h_UU { TensorAlgebra::compute_inverse_sym(vars.h) };
    const Tensor<3, data_t> chris_ULL { TensorAlgebra::compute_christoffel(d1.h, h_UU).ULL };

    //compute EM Tensor
    const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(vars, d1, h_UU, chris_ULL);

}



#endif //DIAGNOSTIC_IMPL_H_INCLUDED