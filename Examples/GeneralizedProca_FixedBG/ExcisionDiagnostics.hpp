/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONDIAGNOSTICS_HPP_
#define EXCISIONDIAGNOSTICS_HPP_

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
template <class background_t>
class ExcisionDiagnostics
{
  protected:
        const double m_dx; //grid spacing
        const double m_excision_width;
        const std::array<double, CH_SPACEDIM> m_center; //center of BH
        
        background_t m_background;

  public:

        //constructor
        ExcisionDiagnostics(background_t a_background, const double a_dx, const std::array<double, CH_SPACEDIM> a_center, double a_excision_cut = 1): m_background{a_background}, m_dx{a_dx}, m_center{a_center}, m_excision_width{a_excision_cut}
        {
        };

        template <class data_t>
        void compute(const Cell<data_t> current_cell) const
        {
            data_t horizon_distance { m_background.excise(current_cell) };

            data_t cell_Inside_Cutoff { (double)simd_compare_lt(horizon_distance, 1.0) };

            if (cell_Inside_Cutoff)
            {
              current_cell.store_vars(0.0, c_gauss);
              current_cell.store_vars(0.0, c_Asquared);
              current_cell.store_vars(0.0, c_gnn);
              current_cell.store_vars(0.0, c_Ham);

            } //excision

        }//end of method def

};

#endif /* EXCISIONDIAGNOSTICS_HPP_ */
