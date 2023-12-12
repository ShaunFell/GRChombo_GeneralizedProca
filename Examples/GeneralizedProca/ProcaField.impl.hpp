/*
 * Implementation file for ProcaField.hpp
 */

#if !defined(PROCAFIELD_H_INCLUDED)
#error "This file should only be included through ProcaField.hpp"
#endif

#ifndef PROCAFIELD_IMPL_H_INCLUDEDP
#define PROCAFIELD_IMPL_H_INCLUDED

// Remove these after debugging
#include "DebuggingTools.hpp"

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> ProcaField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars,          // the value of the variables
    const vars_t<Tensor<1, data_t>> &d1, // the 1st derivatives
    const Tensor<2, data_t> &h_UU,       // the inverse metric
    const Tensor<3, data_t> &chris_ULL   // conformal christoffel symbols
) const
{
    // output energy-momentum
    emtensor_t<data_t> out;

    // compute contravariant and covariant physical spatial metric
    Tensor<2, data_t> gamma_UU;
    Tensor<2, data_t> gamma_LL;
    FOR2(i, j)
    {
        gamma_UU[i][j] = h_UU[i][j] * vars.chi;
        gamma_LL[i][j] = vars.h[i][j] / vars.chi;
    }

    // compute physical christoffel symbols
    auto chris_phys_ULL = TensorAlgebra::compute_phys_chris(
        d1.chi, vars.chi, vars.h, h_UU, chris_ULL);

    // compute potential and its derivatives
    data_t V = 0.0;     // value of potential
    data_t dVdA = 0.0;  // first derivative
    data_t dVddA = 0.0; // second derivative
    m_potential.compute_potential(V, dVdA, dVddA, vars, gamma_UU);

    // D_i A_j : 3-covariant derivative of spatial covector
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys_ULL[k][i][j] * vars.Avec[k]; }
    }

    // antisymmetric tensor D_i A_j - D_j A_i
    Tensor<2, data_t> DA_asym;
    FOR2(i, j) { DA_asym[i][j] = DA[i][j] - DA[j][i]; }

    // the Levi-Civita 3D symbol
    auto LC_symbol = TensorAlgebra::epsilon();

    // the epsilon tensor (3D volume form), indices DOWN and UP
    // remember chi = det(gamma)^(-1/3)
    // then eps = sqrt(det(gamma)) LC_symbol
    // (gamma is the physical 3-metric, h the conformal one)
    data_t sqrt_detgamma = pow(vars.chi, -3. / 2.);
    Tensor<3, data_t> epsilon_LLL;
    Tensor<3, data_t> epsilon_UUU;
    FOR3(i, j, k)
    {
        epsilon_LLL[i][j][k] = sqrt_detgamma * LC_symbol[i][j][k];
        epsilon_UUU[i][j][k] = (1. / sqrt_detgamma) * LC_symbol[i][j][k];
    }

    // TODO: maybe better to forget about the B field?

    // B vector, index UP
    Tensor<1, data_t> Bvec;
    FOR1(i)
    {
        Bvec[i] = 0.0;
        FOR2(j, k) { Bvec[i] += epsilon_UUU[i][j][k] * DA[j][k]; }
    }

    // electric field squared
    data_t Esq = 0.0;
    FOR2(i, j) { Esq += gamma_LL[i][j] * vars.Evec[i] * vars.Evec[j]; }

    // magnetic field squared
    data_t Bsq = 0.0;
    FOR2(i, j) { Bsq += gamma_LL[i][j] * Bvec[i] * Bvec[j]; }

    // ### components of EM tensor ###

    // Eulerian energy \rho
    out.rho = 0.5 * Esq + 0.5 * Bsq + 2 * dVdA * vars.phi * vars.phi + V;

    // Eulerian momentum P_i
    FOR1(i)
    {
        out.Si[i] = 2.0 * vars.phi * dVdA * vars.Avec[i];
        FOR1(j) { out.Si[i] += vars.Evec[j] * DA_asym[i][j]; }
    }

    // Eulerian stress S_ij
    FOR2(i, j)
    {
        out.Sij[i][j] = gamma_LL[i][j] * (0.5 * Esq + 0.5 * Bsq) +
                        2 * dVdA * vars.Avec[i] * vars.Avec[j] -
                        gamma_LL[i][j] * V;
        FOR2(m, n)
        {
            out.Sij[i][j] -= gamma_LL[i][m] * gamma_LL[j][n] *
                             (vars.Evec[m] * vars.Evec[n] + Bvec[m] * Bvec[n]);
        }
    }

    // Eulerian Stress scalar
    out.S = 0.0;
    FOR2(i, j) { out.S += out.Sij[i][j] * gamma_UU[i][j]; }

    return out;
}

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
    const Tensor<2, data_t> h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    const Tensor<3, data_t> chris_ULL =
        TensorAlgebra::compute_christoffel(d1.h, h_UU).ULL;

    // compute physical christoffel symbols
    auto chris_phys_ULL = TensorAlgebra::compute_phys_chris(
        d1.chi, vars.chi, vars.h, h_UU, chris_ULL);

    // calculate physical contravariant and covariant spatial metrics
    Tensor<2, data_t> gamma_UU;
    Tensor<2, data_t> gamma_LL;
    FOR2(i, j)
    {
        gamma_UU[i][j] = h_UU[i][j] * vars.chi;
        gamma_LL[i][j] = vars.h[i][j] / vars.chi;
    }

    // compute potential and its derivatives
    data_t V = 0.0;     // value of potential
    data_t dVdA = 0.0;  // first derivative
    data_t dVddA = 0.0; // second derivative
    m_potential.compute_potential(V, dVdA, dVddA, vars, gamma_UU);

    // the acceleration vector, index DOWN
    Tensor<1, data_t> acceleration;
    FOR1(i) { acceleration[i] = d1.lapse[i] / vars.lapse; }

    // remember: d/dt = lapse * LieD[m] + LieD[shift]

    // evolution equations for spatial part of vector field (index down)
    FOR1(i)
    {
        total_rhs.Avec[i] =
            -vars.lapse * (acceleration[i] * vars.phi + d1.phi[i]) +
            advec.Avec[i];
        FOR1(j)
        {
            total_rhs.Avec[i] += -vars.lapse * gamma_LL[i][j] * vars.Evec[j] +
                                 vars.Avec[j] * d1.shift[j][i];
        }
    }

    // evolution equations for electric vector field (index up)
    // I need:

    // D_i A_j : 3-covariant derivative of spatial covector
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys_ULL[k][i][j] * vars.Avec[k]; }
    }

    // antisymmetric tensor D_i A_j - D_j A_i
    Tensor<2, data_t> DA_asym;
    FOR2(i, j) { DA_asym[i][j] = DA[i][j] - DA[j][i]; }

    // I also calculate the Christoffel LLL and LLU here
    Tensor<3, data_t> chris_phys_LLU;
    Tensor<3, data_t> chris_phys_LLL;

    FOR3(i, j, k)
    {
        chris_phys_LLL[i][j][k] = 0.0;
        FOR1(m)
        {
            chris_phys_LLL[i][j][k] += gamma_LL[i][m] * chris_phys_ULL[m][j][k];
        }
    }

    FOR3(i, j, k)
    {
        chris_phys_LLU[i][j][k] = 0.0;
        FOR1(m)
        {
            chris_phys_LLU[i][j][k] += gamma_UU[m][k] * chris_phys_LLL[i][j][m];
        }
    }

    // D_k (D_i A_j - D_j A_i) : covd of DA_asym, asymm in (ij)
    // (the derivatives of Christoffel symbols cancel!)
    // (always remember the confusing order of indices in d1 and d2 tensors)
    Tensor<3, data_t> DDA_asym;
    FOR3(k, i, j)
    {
        DDA_asym[k][i][j] = d2.Avec[j][i][k] - d2.Avec[i][j][k];

        FOR1(m)
        {
            DDA_asym[k][i][j] += chris_phys_ULL[m][i][k] * d1.Avec[m][j] -
                                 chris_phys_ULL[m][j][k] * d1.Avec[m][i];
            DDA_asym[k][i][j] += chris_phys_ULL[m][j][k] * d1.Avec[i][m] -
                                 chris_phys_ULL[m][i][k] * d1.Avec[j][m];
        }

        FOR2(m, n)
        {
            DDA_asym[k][i][j] += chris_phys_ULL[m][i][k] *
                                     chris_phys_LLU[m][j][n] * vars.Avec[n] -
                                 chris_phys_ULL[m][j][k] *
                                     chris_phys_LLU[m][i][n] * vars.Avec[n];
            DDA_asym[k][i][j] += chris_phys_ULL[m][i][k] *
                                     chris_phys_ULL[n][j][m] * vars.Avec[n] -
                                 chris_phys_ULL[m][j][k] *
                                     chris_phys_ULL[n][i][m] * vars.Avec[n];
        }
    }

    // assembling...
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

        FOR3(a, b, c)
        {
            total_rhs.Evec[i] += vars.lapse * gamma_UU[i][a] * gamma_UU[b][c] *
                                 DDA_asym[b][a][c];
            total_rhs.Evec[i] += vars.lapse * gamma_UU[i][a] * gamma_UU[b][c] *
                                 DA_asym[a][b] * acceleration[c];
        }
    }

    // evolution equation for auxiliary constraint-damping scalar field Z
    total_rhs.Z =
        vars.lapse * (2 * dVdA * vars.phi - m_params.vector_damping * vars.Z) +
        advec.Z;
    FOR1(i)
    {
        total_rhs.Z += vars.lapse * d1.Evec[i][i];
        FOR1(j)
        {
            total_rhs.Z += vars.lapse * chris_phys_ULL[i][i][j] * vars.Evec[j];
        }
    }

    // extrinsic curvature tensor
    Tensor<2, data_t> Kt;
    FOR2(i, j)
    {
        Kt[i][j] = (1. / vars.chi) *
                   (vars.A[i][j] + (1. / 3.) * vars.h[i][j] * vars.K);
    }

    // evolution equation for the scalar part of the Proca field
    data_t gnn = dVdA - 2.0 * dVddA * vars.phi * vars.phi;
    data_t lapse_gnn = vars.lapse / gnn;

    total_rhs.phi =
        -lapse_gnn * (vars.Z - dVdA * vars.phi * vars.K) + advec.phi;
    FOR1(i)
    {
        total_rhs.phi +=
            -lapse_gnn * (-2. * dVddA * vars.phi * vars.Avec[i] * vars.Evec[i]);
        FOR1(j)
        {
            total_rhs.phi -=
                vars.lapse * gamma_UU[i][j] * acceleration[i] * vars.Avec[j];
            total_rhs.phi += -lapse_gnn * dVdA * DA[i][j];
            total_rhs.phi +=
                -lapse_gnn * (-4. * dVddA * vars.phi * gamma_UU[i][j] *
                              vars.Avec[i] * d1.phi[j]);
            FOR2(k, l)
            {
                total_rhs.phi +=
                    -lapse_gnn * 2 * dVddA * gamma_UU[i][k] * gamma_UU[j][l] *
                    (vars.Avec[i] * vars.Avec[j] * DA[l][k] +
                     vars.phi * vars.Avec[i] * vars.Avec[j] * Kt[k][l]);
            }
        }
    }

    // ################################################################################
#ifdef EQUATION_DEBUG_MODE
    DEBUG_HEADER;
    DEBUG_OUT3(V, dVdA, dVddA);
    DEBUG_OUT(vars.lapse);
    DEBUG_OUT(vars.K);
    DEBUG_OUT(total_rhs.Z);
    DEBUG_OUT3(total_rhs.Avec[0], total_rhs.Avec[1], total_rhs.Avec[2]);
    DEBUG_OUT3(total_rhs.Evec[0], total_rhs.Evec[1], total_rhs.Evec[2]);
    DEBUG_OUT2(vars.phi, vars.Z);
    DEBUG_OUT3(vars.Avec[0], vars.Avec[1], vars.Avec[2]);
    DEBUG_OUT3(vars.Evec[0], vars.Evec[1], vars.Evec[2]);
    DEBUG_OUT4(gamma_UU[0][0], gamma_UU[0][1], gamma_UU[0][2], gamma_UU[1][1]);
    DEBUG_OUT2(gamma_UU[1][2], gamma_UU[2][2]);
    FOR2(i, j)
    {
        DEBUG_OUT(d1.Evec[i][j]);
        FOR1(k) { DEBUG_OUT(chris_phys[i][j][k]); }
    }
    DEBUG_END;
#endif // EQUATION_DEBUG_MODE
    // ################################################################################
}

#endif // PROCAFIELD_IMPL_H_INCLUDED
