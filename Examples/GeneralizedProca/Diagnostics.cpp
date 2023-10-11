

#include "Diagnostics.hpp"


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
    //load metric variables from Chombo
    const auto metricvars = current_cell.template load_vars<MetricVars>();
    const auto mattervars = current_cell.template load_vars<MatterVars>();
    const auto mattervarsd1 = m_deriv.template diff1<mattervars>(current_cell);
    const auto metricvarsd1 = m_deriv.template diff1<metricvars>(current_cell);

    const auto h_UU = TensorAlgebra::compute_inverse_sym(metricvars.h); //conformal spatial metric
    const auto chris = TensorAlgebra::compute_christoffel(metricvarsd1.h, h_UU); //conformal christoffel symbs

    const auto chris_phys  = TensorAlgebra::compute_phys_chris(metricvarsd1.chi, metricvars.chi, metricvars.h, h_UU, chris.Ull);
    
    //the return value
    data_t gauss_constraint;

    //compute potential
    data_t V { 0. };
    data_t dVdA { 0. };
    data_t dVddA { 0. };
    m_potential.compute_potential(V, dVdA, dVddA, mattervars, metricvarsd1, metricvars);

    gauss_constraint = 2*dVdA*mattervars.phi;

    FOR1(i)
    {
        gauss_constraint += mattervarsd1.Evec[i][i];

        FOR1(j)
        {
            gauss_constraint += chris_phys.ULL[i][i][j]*mattervars.Evec[j];
        };
    };

    return gauss_constraint;
};

template <class potential_t>
template<class data_t>
void ProcaConstraint<potential_t>::compute(Cell<data_t> current_cell) const
{
    data_t gauss_constraint {constraint_equations(current_cell) };

    current_cell.store_vars(gauss_constraint,c_gauss);
};






template <class data_t>
void ProcaSquared::compute(Cell<data_t> current_cell) const
{
    const auto metricvars = current_cell.template load_vars<MetricVars>();
    const auto mattervars = current_cell.template load_vars<MatterVars>();
    const auto mattervarsd1 = m_deriv.template diff1<mattervars>(current_cell);
    const auto metricvarsd1 = m_deriv.template diff1<metricvars>(current_cell);

    const auto h_UU = TensorAlgebra::compute_inverse_sym(metricvars.h); //conformal spatial metric

    Tensor<2, data_t> vars_gamma;
    FOR2(i,j)
    {
        vars_gamma[i][j] = metricvars.h[i][j]/metricvars.chi;
    };

    data_t Asquared;
    Asquared = -mattervars.phi*mattervars.phi;
    FOR2(i,j)
    {
        Asquared += vars_gamma[i][j] *mattervars.Avec[i]*mattervars.Avec[j];
    };

    current_cell.store_vars(Asquared, c_Asquared);
};