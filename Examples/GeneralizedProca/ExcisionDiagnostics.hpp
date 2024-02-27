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

            //can we have a more sophisticated way to excise inside horizon? Maybe based on conformal factor?
            // We can maybe use the AH-finder to extract spin and use fitted formal in Radia thesis for conformal factor of horizon, then excise cells
            //          whose conformal factor value is less than that of the horizon

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






#ifdef USE_AHFINDER

//Excise matter vars using conformal factor
template <class matter_t> 
class ExcisionDiagnosticsWithChi
{
    // Use matter_t class
    template <class data_t>
    using Vars = typename matter_t::template Vars<data_t>;

    template <class data_t>
    using MetricVars = CCZ4Vars::VarsWithGauge<data_t>;

    protected:
        const double m_dx; //grid spacing
        const std::array<double, CH_SPACEDIM> m_center; //center of BH
        double m_kerr_chi;

    public:

        //constructor
        ExcisionDiagnosticsWithChi(const double a_dx, const std::array<double, CH_SPACEDIM> a_center, double a_kerr_spin = 0.): m_dx{a_dx}, m_center{a_center}
        {
            // Find the approximate location of the Horizon using fit to conformal factor in 
            // https://api.repository.cam.ac.uk/server/api/core/bitstreams/320ef77b-f6ff-426a-852d-00e9c2007940/content
            m_kerr_chi =  0.2666 * sqrt(1 - a_kerr_spin * a_kerr_spin);
        };

        void compute(const Cell<double> current_cell) const
        {
            const Coordinates<double> coords(current_cell, m_dx, m_center);

            double chi_value { current_cell.template load_vars<MetricVars>().chi };

            bool cell_Inside_Horizon {  chi_value < m_kerr_chi }; 

            if (cell_Inside_Horizon)
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


        
       /*  double HorizonChi(double kerr_spin)
        {
            double Horizon_chi { 0.2666 * sqrt(1 - kerr_spin * kerr_spin) };

            return Horizon_chi;
        } */

};//end of class def

#endif //USE_AHFINDER





#ifdef USE_AHFINDER
//Excise matter vars using conformal factor
template <class matter_t> 
class ExcisionDiagnosticsWithAH
{
    // Use matter_t class
    template <class data_t>
    using Vars = typename matter_t::template Vars<data_t>;

    template <class data_t>
    using MetricVars = CCZ4Vars::VarsWithGauge<data_t>;

    protected:
        const double m_dx; //grid spacing
        const std::array<double, CH_SPACEDIM> m_center; //center of BH
        double m_minimal_radius;
        double m_buffer;

    public:

        //constructor
        ExcisionDiagnosticsWithAH(const double a_dx, const std::array<double, CH_SPACEDIM> a_center, double a_minimal_radius = 1., double a_buffer = 1.): m_dx{a_dx}, m_center{a_center}, m_minimal_radius{a_minimal_radius}, m_buffer{a_buffer} {};

        void compute(const Cell<double> current_cell) const
        {
            const Coordinates<double> coords(current_cell, m_dx, m_center);
            double cell_radius { coords.get_radius() };
            double buffered_radius { m_buffer * m_minimal_radius };

            bool cell_Inside_Horizon {  cell_radius < buffered_radius }; 

            if (cell_Inside_Horizon)
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

};//end of class def
#endif //USE_AHFINDER

#endif /* EXCISIONDIAGNOSTICS_HPP_ */
