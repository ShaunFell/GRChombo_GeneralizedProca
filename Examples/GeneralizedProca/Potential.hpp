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


        // V = mass^2 (A^2 + 2 lambda A^2 A^2)
        data_t potPrefactor {pow(m_params.mass,2.0)  };
        V =  potPrefactor * (Asquared + 2.0 * lambda * Asquared * Asquared);
        dVdA = potPrefactor * ( 1.0 + 4.0*lambda*Asquared);
        dVddA = potPrefactor * 4.0 * lambda;

        };

};



















#endif //POTANTIAL_H_INCLUDED