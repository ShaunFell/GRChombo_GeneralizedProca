/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 /*
 Code retrieved from 
 https://github.com/KAClough/GRChombo_public/blob/main/Examples/Oscillaton/HamTaggingCriterion.hpp

 Extraction piece from ChiExtractionTaggingCriterion.hpp
 */

#ifndef CUSTOMTAGGINGCRITERION_HPP_
#define CUSTOMTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

class CustomTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const spherical_extraction_params_t m_params;
    const int m_max_level;
    const int m_level;
    const bool m_activate_extraction;
    const double m_radius;

  public:
    CustomTaggingCriterion(double dx, const int a_level, const double a_L,    
                                  const spherical_extraction_params_t a_params,
                                  const int a_max_level,
                                  const bool activate_extraction = false,
                                  const double a_radius=2) :m_max_level{a_max_level},m_radius{a_radius}, m_dx(dx), m_L{a_L}, m_params{a_params}, m_level{a_level}, m_activate_extraction{activate_extraction} {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        data_t FixedGridCriterion = 0.0;
        data_t ConstraintCriterion = 0.0;
        data_t ExtractionCriterion = 0.0;

        //Fixed Grid tagging
        double ratio = pow(2.0, -(m_level+2.0));
        const data_t max_abs_xy = simd_max(abs(coords.x), abs(coords.y));
        const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
        auto regrid = simd_compare_lt(max_abs_xyz, m_L*ratio);
        FixedGridCriterion = simd_conditional(regrid, 100.0, FixedGridCriterion);

        //Make sure horizon always covered
        if (m_level == m_max_level-1){
            auto regrid_bh = simd_compare_lt(max_abs_xyz, m_radius);
            FixedGridCriterion = simd_conditional(regrid_bh, 100.0, FixedGridCriterion);
        }


        //Extraction radius Tagging
        //If extractin weyl data at given radius, enforce given resolution there
        if (m_activate_extraction)
        {
            for (int iradius {0}; iradius < m_params.num_extraction_radii; ++iradius)
            {
                //regrid if within extraction level
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const data_t r = coords.get_radius();
                    // add 20% buffer to extraction zone to avoid boundary
                    auto regrid = simd_compare_lt(r, 1.2*m_params.extraction_radii[iradius]);
                    ExtractionCriterion = simd_conditional(regrid, 100.0, ExtractionCriterion);
                }
            }
        }

        //Hamiltonian and Gauss tagging
        auto Ham_abs_sum = current_cell.load_vars(c_Ham_abs_sum);
        ConstraintCriterion = sqrt(Ham_abs_sum) * m_dx;

        auto Gauss_abs_sum = current_cell.load_vars(c_gauss);
        ConstraintCriterion *= abs(Gauss_abs_sum)*m_dx;

	//#####################################Temporarily turn off constraint criterion;
	ConstraintCriterion = 0;
	//##############################################################################;

        data_t maxFixedGridExtraction { simd_max(FixedGridCriterion, ExtractionCriterion) };
        data_t criterion { simd_max(maxFixedGridExtraction, ConstraintCriterion) };

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* CUSTOMTAGGINGCRITERION_HPP_ */
