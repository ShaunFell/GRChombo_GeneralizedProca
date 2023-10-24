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
        const FourthOrderDerivatives m_deriv;
        const double m_excision_width;

    public:

        //constructor
        ExcisionProcaEvolution(const double a_dx, double a_excision_width=1.0): m_dx{a_dx}, m_deriv{a_dx}, m_excision_width{a_excision_width} {};

        void compute(const Cell<double> current_cell) const
        {


            
        }








};








#endif //EXCISIONPROCAEVOLUTION_H_INCLUDED