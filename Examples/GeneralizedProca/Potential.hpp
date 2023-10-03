#ifndef POTANTIAL_H_INCLUDED
#define POTANTIAL_H_INCLUDED


//#include "Tensor.hpp" //for tensorial objects
//#include "ADMVariables.hpp" //adm variables

class ProcaPotential{
    public:
        struct params_t{
            double mass;
            double self_interaction;
        };

        //class parameters
        params_t m_params;

        template <typename data_t>
        using MetricVars = typename ADMVars::template Vars<data_t>;

        //constructor for class
        ProcaPotential(params_t a_params): m_params(a_params) {};

        //compute potential
        template <class data_t, template <typename> class vars_t>
        void compute_potential(data_t V, data_t &dVdA, data_t &dVddA,
                                const vars_t<data_t> &vars,
                                const vars_t<Tensor<1,data_t>> &d1,
                                const MetricVars<data_t> &metric_vars) const
        {

        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse(metric_vars.gamma); //inverse of spatial metric calculated using function in TensorAlgebra.hpp
        const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU); ///compute christoffel symbols using function in TensorAlgebra

        const double lambda = m_params.self_interaction; //self interaction coupling

        //compute norm of proca field
        data_t Asquared;
        Asquared = -vars.phi*vars.phi;
        FOR2(i,j) {
            Asquared += gamma_UU[i][j] * vars.Avec[i] * vars.Avec[j];
        }

        // V = mu^2 /2 * A^2 + lambda*mu^2/4 * A^4
        // dVdA = mu^2/2 + lambda*mu^2/2 * A^2
        // dVddA = lambda*mu^2/2
        V = pow(m_params.mass,2.0)/2.0 * Asquared + lambda*pow(m_params.mass,2.0)/4.0 * Asquared*Asquared;
        dVdA = pow(m_params.mass, 2.0)/2.0 * ( 1.0 + lambda*Asquared);
        dVddA = pow(m_params.mass, 2.0)/2.0 * lambda;

        };

};



















#endif //POTANTIAL_H_INCLUDED