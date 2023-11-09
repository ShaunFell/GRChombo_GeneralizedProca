/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONEVOLUTION_HPP_
#define EXCISIONEVOLUTION_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Does excision for fixed BG BH solutions
//! Note that it is does not using simd so one must set disable_simd()
template <class matter_t>
class ExcisionEvolution
{
    using Vars = typename matter_t::template Vars<double>;

  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const double m_inner_r;



  public:
    ExcisionEvolution(const double a_dx,
                        const std::array<double, CH_SPACEDIM> a_center,
                        const double a_inner_r)
        : m_dx(a_dx), m_center(a_center), m_inner_r(a_inner_r)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        if (coords.get_radius() < m_inner_r)
        {
            Vars rhs;
            VarsTools::assign(rhs,0.0);

            current_cell.store_vars(rhs);
        } // else do nothing
    }
};

#endif /* EXCISIONEVOLUTION_HPP_ */
