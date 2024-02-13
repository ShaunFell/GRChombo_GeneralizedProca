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
template <class matter_t> 
class ExcisionDiagnostics
{
    // Use matter_t class
    template <class data_t>
    using Vars = typename matter_t::template Vars<data_t>;

    protected:
        const double m_dx; //grid spacing
        const double m_excision_width;
        const std::array<double, CH_SPACEDIM> m_center; //center of BH

    public:

        //constructor
        ExcisionDiagnostics(const double a_dx, const std::array<double, CH_SPACEDIM> a_center, double a_excision_cut = 1): m_dx{a_dx}, m_center{a_center}, m_excision_width{a_excision_cut}
        {
        };

        template <class data_t>
        void compute(const Cell<data_t> current_cell) const
        {
            const Coordinates<data_t> coords(current_cell, m_dx, m_center);

            //Next line is probably a bug. Coordinates should already be black hole-centered
            //Tensor<1,data_t> coords_BHCentered { coords.x - m_center[0], coords.y-m_center[1], coords.z - m_center[2] };

            //data_t cell_radius_BHCentered { sqrt ( TensorAlgebra::compute_dot_product(coords_BHCentered, coords_BHCentered) ) };
            data_t cell_radius_BHCentered { coords.get_radius() };

            data_t cell_Inside_Cutoff { (double)simd_compare_lt(cell_radius_BHCentered, m_excision_width) };

            if (cell_Inside_Cutoff)
            {
              current_cell.store_vars(0.0, c_gauss);
              current_cell.store_vars(0.0, c_Asquared);
              current_cell.store_vars(0.0, c_gnn);
              current_cell.store_vars(0.0, c_Ham);
              current_cell.store_vars(0.0, c_rho);
              current_cell.store_vars(0.0, c_rhoJ);
              current_cell.store_vars(0.0, c_rhoE);


            } //excision

        }//end of method def

};

#endif /* EXCISIONDIAGNOSTICS_HPP_ */
