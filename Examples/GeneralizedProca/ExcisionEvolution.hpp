#ifndef EXCISIONPROCAEVOLUTION_H_INCLUDED
#define EXCISIONPROCAEVOLUTION_H_INCLUDED

#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "AHInterpolation.hpp"
#include <algorithm>
#include "simd.hpp"



template <class matter_t> 
class ExcisionProcaEvolution
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
        ExcisionProcaEvolution(const double a_dx, const std::array<double, CH_SPACEDIM> a_center, double a_excision_cut = 1): m_dx{a_dx}, m_center{a_center}, m_excision_width{a_excision_cut}
        {
        };

        template <class data_t>
        void compute(const Cell<data_t> current_cell) const
        {
            const Coordinates<data_t> coords(current_cell, m_dx, m_center);

            data_t cell_radius_BHCentered { coords.get_radius() };

            data_t cell_Inside_Cutoff { (double)simd_compare_lt(cell_radius_BHCentered, m_excision_width) };

            if (cell_Inside_Cutoff)
            {
                //matter vars within excision zone
                Vars<data_t> rhs_vars;
                VarsTools::assign(rhs_vars, 0.0);

                //assign values of variables to cell
                current_cell.store_vars(rhs_vars);

            } //excision

        }//end of method def
};//end of class def






#ifdef USE_AHFINDER

//Excise matter vars using conformal factor
template <class matter_t> 
class ExcisionProcaEvolutionWithChi
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
        ExcisionProcaEvolutionWithChi(const double a_dx, const std::array<double, CH_SPACEDIM> a_center, double a_kerr_spin = 0.): m_dx{a_dx}, m_center{a_center}
        {
            m_kerr_chi =  0.2666 * sqrt(1 - a_kerr_spin * a_kerr_spin);
        };

        void compute(const Cell<double> current_cell) const
        {
            const Coordinates<double> coords(current_cell, m_dx, m_center);

            double chi_value { current_cell.template load_vars<MetricVars>().chi };

            bool cell_Inside_Horizon {  chi_value < m_kerr_chi }; 

            if (cell_Inside_Horizon)
            {
                //matter vars within excision zone
                Vars<double> rhs_vars;
                VarsTools::assign(rhs_vars, 0.0);

                //assign values of variables to cell
                current_cell.store_vars(rhs_vars);

            } //excision

        }//end of method def


        // Find the approximate location of the Horizon using fit to conformal factor in 
        // https://api.repository.cam.ac.uk/server/api/core/bitstreams/320ef77b-f6ff-426a-852d-00e9c2007940/content
       /*  double HorizonChi(double kerr_spin)
        {
            double Horizon_chi { 0.2666 * sqrt(1 - kerr_spin * kerr_spin) };

            return Horizon_chi;
        } */
};//end of class def

#endif //USE_AHFINDER



#ifdef USE_AHFINDER
template <class matter_t, class AHinterp_t> 
class ExcisionProcaEvolutionWithAH
{
    // Use matter_t class
    using Vars = typename matter_t::template Vars<double>;


    protected:
        const double m_dx; //grid spacing
        const double m_excision_width;
        const std::array<double, CH_SPACEDIM> m_center; //center of BH
        std::vector<Tensor<1,double>> m_AH_coords;
        AHinterp_t& m_ah_interp;
        const double m_num_AH_points;


    public:

        //constructor
        ExcisionProcaEvolutionWithAH(const double a_dx, const std::array<double, CH_SPACEDIM> a_center, AHinterp_t& a_ah_interp, const double a_num_AH_points, double a_excision_width=1.0): m_dx{a_dx}, m_center{a_center}, m_excision_width{a_excision_width}, m_ah_interp{a_ah_interp}, m_num_AH_points{a_num_AH_points}
        {
            for (int i{0}; i < m_num_AH_points; ++i){
                auto point { m_ah_interp.get_cartesian_coords(i) };
                m_AH_coords.push_back(
                    point
                );
            };
        };

        void compute(const Cell<double> current_cell) const
        {
            const Coordinates<double> coords(current_cell, m_dx, m_center);

            std::vector<double> distance_to_AH_points;
            double min_distance;
            Tensor<1,double> closest_AH_point;

            for (int i{0}; i < m_num_AH_points; ++i)
            {
                distance_to_AH_points.push_back(
                    sqrt(
                        (m_AH_coords[i][0]  - coords.x)*(m_AH_coords[i][0]  - coords.x) +
                        (m_AH_coords[i][1]  - coords.y)*(m_AH_coords[i][1]  - coords.y) +
                        (m_AH_coords[i][2]  - coords.z)*(m_AH_coords[i][2]  - coords.z)
                    )
                );

                if ( *min_element(distance_to_AH_points.begin(), distance_to_AH_points.end()) == distance_to_AH_points.back())
                {
                    min_distance = distance_to_AH_points.back();
                    closest_AH_point = m_AH_coords[i];
                }
            }

            double AH_coord_BH_Centered_Norm { 
                                                sqrt(
                                                    closest_AH_point[0]*closest_AH_point[0] +
                                                    closest_AH_point[1]*closest_AH_point[1] +
                                                    closest_AH_point[2]*closest_AH_point[2]
                                                )
                                                };
            double cell_BH_centered_Norm {coords.get_radius()};
            
            bool cell_Inside_Horizon { cell_BH_centered_Norm <= AH_coord_BH_Centered_Norm };
            bool cell_Inside_Buffered_Horizon { cell_BH_centered_Norm <= m_excision_width * AH_coord_BH_Centered_Norm };

            if (cell_Inside_Horizon)
            {
                current_cell.store_vars(0.0, c_Z);

                if (cell_Inside_Buffered_Horizon)
                {
                    //matter vars within excision zone
                    Vars rhs_vars;
                    VarsTools::assign(rhs_vars, 0.0);

                    //assign values of variables to cell
                    current_cell.store_vars(rhs_vars);

                } //matter excision
            } //auxiliary field excision

        }//end of method def

};//end of class def
#endif //USE_AHFINDER


#endif //EXCISIONPROCAEVOLUTION_H_INCLUDED