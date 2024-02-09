#ifndef POTANTIAL_H_INCLUDED
#define POTANTIAL_H_INCLUDED


#include "Tensor.hpp" //for tensorial objects


class ProcaPotential{
    public:
        struct params_t{
            double mass;
            double self_interaction;
        };

        //class parameters
        params_t m_params;

        //constructor for class
        ProcaPotential(params_t a_params): m_params(a_params) {};

        //compute potential
        template <class data_t, template <typename> class vars_t>
        void compute_potential(data_t &V, data_t &dVdA, data_t &dVddA,
                                const vars_t<data_t> &vars,
                                const Tensor<2, data_t> &gamma_UU) const
        {

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

    //################################################################################
#ifdef EQUATION_DEBUG_MODE
        DEBUG_OUT(V);
        DEBUG_OUT(dVdA);
        DEBUG_OUT(dVddA);
        DEBUG_OUT(Asquared);
        DEBUG_OUT3(vars.Avec[0], vars.Avec[1], vars.Avec[2]);
        DEBUG_OUT(vars.phi);
        DEBUG_OUT4(gamma_UU[0][0], gamma_UU[0][1], gamma_UU[0][2], gamma_UU[1][1]);
        DEBUG_OUT2(gamma_UU[1][2], gamma_UU[2][2]);
#endif //EQUATION_DEBUG_MODE
    //################################################################################

        };

};



















#endif //POTANTIAL_H_INCLUDED