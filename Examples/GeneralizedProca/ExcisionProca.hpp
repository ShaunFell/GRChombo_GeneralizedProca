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


template <class matter_t, AHinterp_t> 
class ExcisionProcaEvolution
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
        ExcisionProcaEvolution(const double a_dx, const std::array<double, CH_SPACEDIM> a_center, AHinterp_t& a_ah_interp, const double a_num_AH_points, double a_excision_width=1.0): m_dx{a_dx}, m_center{a_center}, m_excision_width{a_excision_width}, m_ah_interp{a_ah_interp}, m_num_AH_points{a_num_AH_points}
        {
            for (int i{0}; i < m_bh_amr.m_ah_finder.get(0).m_params.num_points_u * m_bh_amr.m_ah_finder.get(0).m_params.num_points_v; ++i){
                m_AH_coords.push_back(
                    m_bh_amr.m_ah_finder.get(0).get_ah_interp().get_cartesian_coords(i);
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
                        (m_AH_coords[i][0] - coords.x)*(m_AH_coords[i][0] - coords.x) +
                        (m_AH_coords[i][1] - coords.y)*(m_AH_coords[i][1] - coords.y) +
                        (m_AH_coords[i][2] - coords.z)*(m_AH_coords[i][2] - coords.z)
                    )
                );
                if (distance_to_AH_points.min() == distance_to_AH_points.back())
                {
                    min_distance = distance_to_AH_points.back();
                    closest_AH_point = m_AH_coords[i];
                }
            }
            double AH_coord_BH_Centered_Norm { 
                                                sqrt(
                                                    closest_AH_point[0] - m_center[0] +
                                                    closest_AH_point[1] - m_center[1] +
                                                    closest_AH_point[2] - m_center[2]
                                                )
                                                };
            double cell_BH_centered_Norm {
                                            sqrt(
                                                coords.x - m_center[0] +
                                                coords.y - m_center[1] +
                                                coords.z - m_center[2]
                                            )
                                            };
            bool cell_Inside_Horizon { cell_BH_centered_Norm <= AH_coord_BH_Centered_Norm};
            bool cell_Inside_Buffered_Horizon { cell_BH_centered_Norm <= m_excision_width * AH_coord_BH_Centered_Norm};

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

                } //matter variables
            } //excision
        }//end of method def
};//end of class def
#endif //EXCISIONPROCAEVOLUTION_H_INCLUDED