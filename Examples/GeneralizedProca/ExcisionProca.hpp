#ifndef EXCISIONPROCAEVOLUTION_H_INCLUDED
#define EXCISIONPROCAEVOLUTION_H_INCLUDED

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"


template <class matter_t> 
class ExcisionProcaEvolution
{
    // Use matter_t class
    using Vars = typename matter_t::template Vars<double>;

    protected:
        const double m_dx; //grid spacing
        const double m_excision_width;
        const std::array<double, CH_SPACEDIM> m_center; //center of BH
        double m_horizon_size;


    public:

        //constructor
        ExcisionProcaEvolution(const double a_dx, const std::array<double, SPACEDIM> a_center, double a_excision_width=1.0, double a_horizon_size): m_dx{a_dx}, m_center{a_center}, m_excision_width{a_excision_width}, m_horizon_size{a_horizon_size} {};

        void compute(const Cell<double> current_cell) const
        {
            const Coordinates<double> coords(current_cell, m_dx, m_center);

            double R { sqrt(coords.x*coords.x + coords.y*coords.y + coords.z*coords.z) };

            if (R < m_horizon_size)
            {
                current_cell.store_vars(0.0, Vars.Z);

                if (R < m_excision_width*m_horizon_size)
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