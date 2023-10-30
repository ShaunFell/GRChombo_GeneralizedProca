/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 /*
 Code retrieved from 
 https://github.com/KAClough/GRChombo_public/blob/main/Examples/Oscillaton/HamTaggingCriterion.hpp

 Extraction piece from ChiExtractionTaggingCriterion.hpp
 */

#ifndef HAMEXTRACTIONTAGGINGCRITERION_HPP_
#define HAMEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

class HamExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const spherical_extraction_params_t m_params;
    const int m_level;
    const bool m_activate_extraction;

  public:
    HamExtractionTaggingCriterion(double dx, const int a_level,
                                  const spherical_extraction_params_t a_params,
                                  const bool activate_extraction = false) : m_dx(dx), m_params{a_params}, m_level{a_level}, m_activate_extraction{activate_extraction} {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto Ham_abs_sum = current_cell.load_vars(c_Ham_abs_sum);
        data_t criterion = sqrt(Ham_abs_sum) * m_dx;

        //If extractin weyl data at given radius, enforce given resolution there
        if (m_activate_extraction)
        {
            for (int iradius {0}; iradius < m_params.num_extraction_radii; ++iradius)
            {
                //regrid if within extraction level
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
                    const data_t r = coords.get_radius();
                    // add 20% buffer to extraction zone to avoid boundary
                    auto regrid = simd_compare_lt(r, 1.2*m_params.extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
                }
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* HAMEXTRACTIONTAGGINGCRITERION_HPP_ */
