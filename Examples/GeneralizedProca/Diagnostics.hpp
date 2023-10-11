#ifndef DIAGNOSTIC_H_INCLUDED
#define DIAGNOSTIC_H_INCLUDED

#include "ADMVariables.hpp"
#include "CCZ4Geometry.hpp"
#include "ProcaField.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Potential.hpp"
#include "ProcaField.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

//calculate constraint violation of the Proca field

template <class potential_t>
class ProcaConstraint
{

    protected:
        double m_vector_mass;
        double m_vector_damping;
        const FourthOrderDerivatives m_deriv;
        
        const potential_t m_potential;

    public:
        //constructor
        ProcaConstraint(double dx, double a_vector_mass, double a_vector_damping, const potential_t potential);

        // use variable definition in Proca class
        template<class data_t>
        using MatterVars = typename ProcaField<potential_t>::template Vars<data_t>;

        //ADM vars
        template<class data_t>
        using MetricVars = typename ADMVars::template Vars<data_t>;

        
        //calculate constraint equations
        template<class data_t>
        data_t constraint_equations(Cell<data_t> current_cell) const;

        template<class data_t>
        void compute(Cell<data_t> current_cell) const;
};



class ProcaSquared
{
    protected:

        const FourthOrderDerivatives m_deriv;
        template<class data_t>
        using MetricVars = ADMVars::Vars<data_t>;

        template<class data_t>
        using MatterVars = ProcaField<ProcaPotential>::template Vars<data_t>;

    public:
        ProcaSquared(double a_dx): m_deriv(a_dx){};

        template<class data_t>
        void compute(Cell<data_t> current_cell) const;
};



#endif //DIAGNOSTIC_H_INCLUDED

