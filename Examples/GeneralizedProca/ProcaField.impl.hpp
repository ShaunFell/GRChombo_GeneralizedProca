/*
implementation file for ProcaField.hpp
*/

#if !defined(PROCAFIELD_H_INCLUDED)
#error "This file should only be included through ProcaField.hpp"
#endif

#ifndef PROCAFIELD_IMPL_H_INCLUDEDP
#define PROCAFIELD_IMPL_H_INCLUDED

// Remove these after debugging
#include "DebuggingTools.hpp"

template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> ProcaField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars,          // the value of the variables
    const vars_t<Tensor<1, data_t>> &d1, // the 1st derivatives
    const Tensor<2, data_t> &h_UU,       // the inverse metric
    const Tensor<3, data_t> &chris_ULL   // conformal christoffel symbols
) const
{
    emtensor_t<data_t> out;

    // compute contravariant physical spatial metric
    Tensor<2, data_t> gamma_UU;
    FOR2(i, j) { gamma_UU[i][j] = h_UU[i][j] * vars.chi; };

    // compute covariant physical spatial metric
    Tensor<2, data_t> gamma_LL{TensorAlgebra::compute_inverse_sym(gamma_UU)}; // inverse of inverse metric -> pure covariant spatial metric

    // compute physical christoffel symbols
    Tensor<3, data_t> chris_phys_ULL{TensorAlgebra::compute_phys_chris(
        d1.chi, vars.chi, vars.h, h_UU, chris_ULL)};

    // Compute potential and its derivatives
    data_t V{0.};     // value of potential
    data_t dVdA{0.};  // first derivative of Potential function w.r.t argument
    data_t dVddA{0.}; // second derivative...
    m_potential.compute_potential(V, dVdA, dVddA, vars, gamma_UU);

    // D_i A_j  3-covariant derivative of spatial covector
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys_ULL[k][i][j] * vars.Avec[k]; };
    };

    // D_i A_j - D_j A_i
    Tensor<2, data_t> DA_antisym;
    FOR2(i, j) { DA_antisym[i][j] = d1.Avec[j][i] - d1.Avec[i][j]; };

    // Electric Field Norm
    data_t Enorm{0};
    FOR2(i, j) { Enorm += gamma_LL[i][j] * vars.Evec[i] * vars.Evec[j]; };

    /////Components of EM tensor

    // Eulerian Energy //Checked.
    out.rho = 1. / 2. * Enorm + dVdA * vars.phi * vars.phi + 0.5*V;
    FOR4(i, j, k, l)
    {
        out.rho += 1. / 2. * gamma_UU[k][i] * gamma_UU[l][j] * DA[i][j] *
                   DA_antisym[k][l];
    };

    // Eulerian Momentum //Checked.
    FOR1(i)
    {
        out.Si[i] = 0; // zero initialize
        out.Si[i] +=  vars.phi * dVdA * vars.Avec[i];

        FOR1(j) { out.Si[i] += vars.Evec[j] * DA_antisym[i][j]; };
    };

    // Eulerian Stress //Checked. 
    FOR2(i, j)
    {
        out.Sij[i][j] = 0; // zero initialize

        out.Sij[i][j] +=  dVdA * vars.Avec[i] * vars.Avec[j] -
                         gamma_LL[i][j] * V + 0.5 * gamma_LL[i][j] * Enorm;

        FOR2(l, k)
        {
            out.Sij[i][j] +=
                -gamma_LL[i][l] * gamma_LL[j][k] * vars.Evec[l] * vars.Evec[k] +
                gamma_UU[k][l] * DA_antisym[i][l] * DA_antisym[j][k];

            FOR2(m, n)
            {
                out.Sij[i][j] += - 0.5 * gamma_LL[i][j] * gamma_UU[m][l] *
                                 gamma_UU[n][k] * DA[m][n] * DA_antisym[l][k];
            };
        };
    };

    // Eulerian Stress scalar
    out.S = 0.0;
    FOR2(i, j) { out.S += out.Sij[i][j] * gamma_UU[i][j]; };

    return out;
};

template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void ProcaField<potential_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs,             // RHS terms for all vars
    const vars_t<data_t> &vars,                // the value of the variables
    const vars_t<Tensor<1, data_t>> &d1,       // value of 1st derivs
    const diff2_vars_t<Tensor<2, data_t>> &d2, // 2nd derivs
    const vars_t<data_t> &advec // value of the beta^i d_i(var) terms
) const
{

    // calculate conformal contravariant metric and conformal christoffel
    // symbols
    const Tensor<2, data_t> h_UU = TensorAlgebra::compute_inverse(vars.h);
    const Tensor<3, data_t> chris_ULL =
        TensorAlgebra::compute_christoffel(d1.h, h_UU).ULL;

    // compute physical christoffel symbols
    Tensor<3, data_t> chris_phys{TensorAlgebra::compute_phys_chris(
        d1.chi, vars.chi, vars.h, h_UU, chris_ULL)};

    // calulate physical contravariant spatial metric
    Tensor<2, data_t> gamma_LL;
    FOR2(i, j) { gamma_LL[i][j] = vars.h[i][j] / vars.chi; }

    // physical covariant spatial metric
    Tensor<2, data_t> gamma_UU{TensorAlgebra::compute_inverse_sym(gamma_LL)};

    // compute potential and its derivatives
    data_t V{0.};
    data_t dVdA{0.};
    data_t dVddA{0.};
    m_potential.compute_potential(V, dVdA, dVddA, vars, gamma_UU);

    // covariant derivative of spatial part of Proca field
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys[k][i][j] * vars.Avec[k]; }
    }

    // Extrinsic curvature
    Tensor<2, data_t> ExCurv;
    FOR2(i, j)
    {
        ExCurv[i][j] =
            (1. / vars.chi) * (vars.A[i][j] + 1. / 3. * vars.h[i][j] * vars.K);
    }

    // evolution equation for the scalar part of the Proca field
    data_t gnn{dVdA - 2.0 * dVddA * vars.phi * vars.phi};
    data_t mass{m_potential.m_params.mass};

    //Calculations we want to do only once
    Tensor<1,data_t> Elow { 0., 0., 0. };
    Tensor<2,data_t> DA_antisym;
    FOR2(i,j)
    {
        Elow[i] += gamma_LL[i][j] * vars.Evec[j];
        DA_antisym[i][j] = d1.Avec[j][i] - d1.Avec[i][j];
    }

    // evolution equations for spatial part of vector field (index down)
    FOR1(i)
    {
        total_rhs.Avec[i] =
            -vars.lapse * d1.phi[i] - vars.phi * d1.lapse[i] + advec.Avec[i];

        FOR1(j)
        {
            total_rhs.Avec[i] += -vars.lapse * Elow[i] +
                                 vars.Avec[j] * d1.shift[j][i];
        };
    };

    // evolution equations for Electric vector field (index up)
    FOR1(i)
    {
        total_rhs.Evec[i] = vars.lapse * vars.K * vars.Evec[i] + advec.Evec[i];

        FOR1(j)
        {
            total_rhs.Evec[i] +=
                vars.lapse * gamma_UU[i][j] * d1.Z[j] +
                2 * vars.lapse * dVdA * gamma_UU[i][j] * vars.Avec[j] -
                vars.Evec[j] * d1.shift[i][j];
        }

        FOR3(j, k, l)
        {
            total_rhs.Evec[i] += gamma_UU[j][k] * gamma_UU[i][l] * (
                                    d1.lapse[j] * DA_antisym[l][k] +
                                    vars.lapse * ( d2.Avec[k][l][j] - d2.Avec[l][k][j] )
                                );

            FOR1(m)
            {
                total_rhs.Evec[i] +=
                    -vars.lapse * gamma_UU[j][k] * gamma_UU[i][l] *
                    (
                        chris_phys[m][j][l] * DA_antisym[m][k] +
                        chris_phys[m][j][k] * DA_antisym[l][m]
                     );
            };
        };
    };

    // evolution equation for auxiliary constraint-damping scalar field Z
    total_rhs.Z = 2 * vars.lapse * dVdA * vars.phi -
                  m_params.vector_damping * vars.lapse * vars.Z + advec.Z;
    FOR1(i)
    {
        total_rhs.Z += vars.lapse * d1.Evec[i][i];
        FOR1(j)
        {
            total_rhs.Z += vars.lapse * chris_phys[i][i][j] * vars.Evec[j];
        }
    }

    total_rhs.phi =  - vars.lapse * vars.Z * mass * mass / (2 * gnn) +
                    vars.lapse * dVdA * vars.phi * vars.K / (gnn) + advec.phi;
    FOR1(i)
    {
        total_rhs.phi += 2 * vars.lapse * dVddA * vars.phi * vars.Avec[i] *
                         vars.Evec[i] / gnn;

        FOR1(j)
        {
            total_rhs.phi +=
                gamma_UU[i][j] * (-vars.lapse * dVdA / gnn * DA[i][j] -
                                  vars.Avec[i] * d1.lapse[j] +
                                  2 * vars.lapse * dVddA / gnn * 2 * vars.phi *
                                      vars.Avec[i] * d1.phi[j]);

            FOR2(k, l)
            {
                total_rhs.phi -=
                    2 * vars.lapse * dVddA / gnn * gamma_UU[i][k] * gamma_UU[j][l] *
                    ( 
                        vars.phi * vars.Avec[i] * vars.Avec[j] * ExCurv[k][l]
                      + vars.Avec[i] * vars.Avec[j] * DA[k][l]
                      );
            }
        }
    }
}

#endif // PROCAFIELD_IMPL_H_INCLUDED
