#ifndef EXCISIONPROCAEVOLUTION_H_INCLUDED
#define EXCISIONPROCAEVOLUTION_H_INCLUDED

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "AHInterpolation.hpp"
#include <algorithm>
#include "simd.hpp"


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
            pout() << "Num AH points: " << m_num_AH_points <<endl;
            for (int i{0}; i < m_num_AH_points; ++i){
                m_AH_coords.push_back(
                    m_ah_interp.get_cartesian_coords(i)
                );
            };
            pout() << "Points extracted" << endl;
        };

        void compute(const Cell<double> current_cell) const
        {
            const Coordinates<double> coords(current_cell, m_dx, m_center);

            std::vector<double> distance_to_AH_points;
            double min_distance;
            Tensor<1,double> closest_AH_point;

            pout() << "Finding closest AH point"<<endl;
            pout() << "m_AH_coords size: " << m_AH_coords.size() << endl;
            pout() << "m_num_AH_points: " << m_num_AH_points << endl;
            for (int i{0}; i < m_num_AH_points; ++i)
            {
                pout () << i << " pushing to distance_to_AH_point" << endl;
                distance_to_AH_points.push_back(
                    sqrt(
                        (m_AH_coords[i][0] + m_center[0] - coords.x)*(m_AH_coords[i][0] + m_center[0] - coords.x) +
                        (m_AH_coords[i][1] + m_center[1]  - coords.y)*(m_AH_coords[i][1] + m_center[1]  - coords.y) +
                        (m_AH_coords[i][2] + m_center[2]  - coords.z)*(m_AH_coords[i][2] + m_center[2]  - coords.z)
                    )
                );
                pout() << i << " minimum check" << endl;
                if ( *min_element(distance_to_AH_points.begin(), distance_to_AH_points.end()) == distance_to_AH_points.back())
                {
                    min_distance = distance_to_AH_points.back();
                    closest_AH_point = m_AH_coords[i];
                }
                pout() << i << "continue" <<endl;
            }
            pout() << "Found point. Determining norms"<<endl;
            double AH_coord_BH_Centered_Norm { 
                                                sqrt(
                                                    closest_AH_point[0]*closest_AH_point[0] +
                                                    closest_AH_point[1]*closest_AH_point[1] +
                                                    closest_AH_point[2]*closest_AH_point[2]
                                                )
                                                };
            double cell_BH_centered_Norm {
                                            sqrt(
                                                (coords.x - m_center[0])*(coords.x - m_center[0]) +
                                                (coords.y - m_center[1])*(coords.y - m_center[1]) +
                                                (coords.z - m_center[2])*(coords.z - m_center[2])
                                            )
                                            };
            pout() << "AH_coord_BH_Centered_Norm: " << AH_coord_BH_Centered_Norm << endl;
            pout() << "cell_BH_centered_Norm: " << cell_BH_centered_Norm << endl;

            bool cell_Inside_Horizon { cell_BH_centered_Norm <= AH_coord_BH_Centered_Norm};
            bool cell_Inside_Buffered_Horizon { cell_BH_centered_Norm <= m_excision_width * AH_coord_BH_Centered_Norm};
            pout() << "cell_Inside_Horizon " << cell_Inside_Horizon <<endl;
            pout() << "cell_Inside_Buffered_Horizon " << cell_Inside_Buffered_Horizon <<endl;
            pout() << "Perform Excision"<<endl;

            if (cell_Inside_Horizon)
            {
                pout() << "Excising Z" <<endl;
                current_cell.store_vars(0.0, c_Z);

                if (cell_Inside_Buffered_Horizon)
                {
                    //matter vars within excision zone
                    Vars rhs_vars;
                    pout() << "Excising Vars" <<endl;
                    VarsTools::assign(rhs_vars, 0.0);

                    //assign values of variables to cell
                    pout() << "Storing vars" <<endl;
                    current_cell.store_vars(rhs_vars);

                } //matter variables
            } //excision

            pout() <<"Finished excision." <<endl;
        }//end of method def
};//end of class def


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

            //Next line is probably a bug. Coordinates should already be black hole-centered
            //Tensor<1,data_t> coords_BHCentered { coords.x - m_center[0], coords.y-m_center[1], coords.z - m_center[2] };

            //data_t cell_radius_BHCentered { sqrt ( TensorAlgebra::compute_dot_product(coords_BHCentered, coords_BHCentered) ) };
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

#endif //EXCISIONPROCAEVOLUTION_H_INCLUDED