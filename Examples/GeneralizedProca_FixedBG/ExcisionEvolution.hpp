#ifndef EXCISIONPROCAEVOLUTION_H_INCLUDED
#define EXCISIONPROCAEVOLUTION_H_INCLUDED

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <algorithm>




template <class matter_t, class background_t> 
class ExcisionProcaEvolution
{
    // Use matter_t class
    template <class data_t>
    using Vars = typename matter_t::template Vars<data_t>;

    protected:
        const double m_dx; //grid spacing
        const double m_excision_width;
        const std::array<double, CH_SPACEDIM> m_center; //center of BH
        
        background_t m_background;

    public:

        //constructor
        ExcisionProcaEvolution(background_t a_background, const double a_dx, const std::array<double, CH_SPACEDIM> a_center, double a_excision_cut = 1): m_background{a_background}, m_dx{a_dx}, m_center{a_center}, m_excision_width{a_excision_cut}
        {
        };

        template <class data_t>
        void compute(const Cell<data_t> current_cell) const
        {
            Vars<data_t> rhs_vars;

            data_t horizon_distance { m_background.excise(current_cell) };

            data_t cell_Inside_Cutoff { (double)simd_compare_lt(horizon_distance, 1.0) };
            data_t cell_Inside_Padded_Cutoff { (double)simd_compare_lt(horizon_distance, 0.95*m_excision_width) };

            if (cell_Inside_Cutoff)
            {
                //matter vars within excision zone
                current_cell.store_vars(0.0,c_Z);

                if (cell_Inside_Padded_Cutoff)
                {
                    VarsTools::assign(rhs_vars, 0.0);

                    //assign values of variables to cell
                    current_cell.store_vars(rhs_vars);
                }

            } //excision

        }//end of method def
};//end of class def

#endif //EXCISIONPROCAEVOLUTION_H_INCLUDED