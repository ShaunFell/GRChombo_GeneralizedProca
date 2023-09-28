/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "Potential.hpp"

//! Class which sets the initial scalar field matter config
class InitialScalarData
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //!< Amplitude of bump in initial SF bubble
        std::array<double, CH_SPACEDIM>
            center;   //!< Centre of perturbation in initial SF bubble
        double width; //!< Width of bump in initial SF bubble
    };

    //! The constructor
    InitialScalarData(params_t a_params, Potential::params_t b_params, double a_dx)
        : m_dx(a_dx), m_params(a_params), m_paramsPotential{b_params}
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        data_t rr = coords.get_radius();
        data_t rr2 = rr * rr;

        // calculate the field value
        data_t x = coords.x;
        data_t lambda = m_paramsPotential.lambda;
        data_t eta = m_paramsPotential.eta;
        data_t phi = eta*tanh( pow(lambda/2.0,0.5) * eta * x );

        // store the vars
        current_cell.store_vars(phi, c_phi);
        current_cell.store_vars(0.0, c_Pi);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
    const Potential::params_t m_paramsPotential;
};

#endif /* INITIALSCALARDATA_HPP_ */
