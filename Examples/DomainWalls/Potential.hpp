/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        double eta;
        double lambda;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        V_of_phi = m_params.lambda / 4.0 * pow((vars.phi * vars.phi - m_params.eta*m_params.eta), 2);

        // The potential gradient at phi
        dVdphi = m_params.lambda * vars.phi * (vars.phi * vars.phi - m_params.eta*m_params.eta);
    }
};

#endif /* POTENTIAL_HPP_ */
