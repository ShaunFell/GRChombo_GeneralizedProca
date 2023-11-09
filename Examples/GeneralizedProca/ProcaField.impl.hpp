/*
implementation file for ProcaField.hpp
*/

#if !defined(PROCAFIELD_H_INCLUDED)
#error "This file should only be included through ProcaField.hpp"
#endif 

#ifndef PROCAFIELD_IMPL_H_INCLUDEDP
#define PROCAFIELD_IMPL_H_INCLUDED

//Remove these after debugging
#include "DebuggingTools.hpp"


template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> ProcaField<potential_t>::compute_emtensor(
        const vars_t<data_t> &vars, //the value of the variables
        const vars_t<Tensor<1,data_t>> &d1, //the 1st derivatives
        const Tensor<2, data_t> &h_UU, //the inverse metric
        const Tensor<3, data_t> &chris_ULL //conformal christoffel symbols
    ) const {
        emtensor_t<data_t> out;

        //compute contravariant physical spatial metric
        Tensor<2,data_t> gamma_UU;
        FOR2(i,j){
            gamma_UU[i][j] = h_UU[i][j]/vars.chi;
        };

        //compute covariant physical spatial metric
        Tensor<2, data_t> gamma_LL { TensorAlgebra::compute_inverse_sym(gamma_UU) }; //inverse of inverse metric -> pure covariant spatial metric

        //compute physical christoffel symbols
        Tensor<3, data_t> chris_phys_ULL {TensorAlgebra::compute_phys_chris(d1.chi, vars.chi, vars.h, h_UU, chris_ULL) };

        //Compute potential and its derivatives
        data_t V {0.}; //value of potential
        data_t dVdA {0.}; //first derivative of Potential function w.r.t argument
        data_t dVddA {0.}; //second derivative...
        m_potential.compute_potential(V, dVdA, dVddA, vars, gamma_UU);
        
        // D_i A_j  3-covariant derivative of spatial covector
        Tensor<2, data_t> DA;
        FOR2(i,j){
            DA[i][j] = d1.Avec[j][i];
            FOR1(k){
                DA[i][j] -= chris_phys_ULL[k][i][j]*vars.Avec[k];
            };
        };

        //D_i A_j - D_j A_i
        Tensor<2,data_t> DA_antisym;
        FOR2(i,j){
            DA_antisym[i][j] = d1.Avec[j][i] - d1.Avec[i][j];
        };


        //Electric Field Norm
        data_t Enorm {0};
        FOR2(i,j){
            Enorm += gamma_LL[i][j]*vars.Evec[i]*vars.Evec[j];
        };


        /////Components of EM tensor

        //Eulerian Energy 
        out.rho = 1/2 * Enorm + 2*dVdA*vars.phi*vars.phi + V;
        FOR4(i,j,k,l){
            out.rho += 1/2 * gamma_UU[k][i]*gamma_UU[l][j]*DA[i][j]*DA_antisym[k][l];
        };


        //Eulerian Momentum
        FOR1(i){
            out.Si[i] = 0; //zero initialize
            out.Si[i] += 2*vars.phi*dVdA*vars.Avec[i];
            
            FOR1(j){
                out.Si[i] += -vars.Evec[j]*DA_antisym[j][i];
            };
        };

        //Eulerian Stress
        FOR2(i,j){
            out.Sij[i][j] = 0; //zero initialize

            out.Sij[i][j] +=   2*dVdA*vars.Avec[i]*vars.Avec[j] - gamma_LL[i][j]*V + 1/2 * gamma_LL[i][j]*Enorm;

            FOR2(l,k){
                out.Sij[i][j] += -gamma_LL[i][l]*gamma_LL[j][k]*vars.Evec[l]*vars.Evec[k] - gamma_UU[k][l]*DA_antisym[i][l]*DA_antisym[k][j];

                FOR2(m,n){
                    out.Sij[i][j] += -1/2 * gamma_LL[i][j]*gamma_UU[m][l]*gamma_UU[n][k]*DA[m][n]*DA_antisym[l][k];
                };
            };
        };

        //Eulerian Stress scalar
        out.S = 0.0;
        FOR2(i,j){
            out.S += out.Sij[i][j]*gamma_UU[i][j];
        };

	    //################################################################################
#ifdef EQUATION_DEBUG_MODE
        DEBUG_OUT(out.S);
        DEBUG_OUT(out.rho);
        FOR1(i){
		    DEBUG_OUT(out.Si[i]);
        };
        FOR2(i,j){
            DEBUG_OUT(out.Sij[i][j]);
        };
#endif //EQUATION_DEBUG_MODE
    //################################################################################

        return out;
    };


template <class potential_t>
template <class data_t, 
            template <typename> class vars_t, 
            template <typename> class diff2_vars_t, 
            template <typename> class rhs_vars_t>
void ProcaField<potential_t>::add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs, //RHS terms for all vars
        const vars_t<data_t> &vars, //the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //value of 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //2nd derivs
        const vars_t<data_t> &advec //value of the beta^i d_i(var) terms
        ) const
{



    //calculate conformal contravariant metric and conformal christoffel symbols
    const Tensor<2,data_t> h_UU = TensorAlgebra::compute_inverse(vars.h);
    const Tensor<3, data_t> chris_ULL = TensorAlgebra::compute_christoffel(d1.h, h_UU).ULL;

    //compute physical christoffel symbols
    Tensor<3, data_t> chris_phys {TensorAlgebra::compute_phys_chris(d1.chi, vars.chi, vars.h, h_UU, chris_ULL) };

    //calulate physical contravariant spatial metric
    Tensor<2, data_t> gamma_LL;
    FOR2(i,j){
        gamma_LL[i][j] = vars.h[i][j]/vars.chi;
    }

    //physical covariant spatial metric
    Tensor<2, data_t> gamma_UU { TensorAlgebra::compute_inverse_sym(gamma_LL) };


    //compute potential and its derivatives
    data_t V {0.};
    data_t dVdA {0.};
    data_t dVddA {0.};
    m_potential.compute_potential(V, dVdA, dVddA, vars, gamma_UU);


    const double c4 = m_potential.self_interaction;

    FOR1(i)
    {
        total_rhs.Avec[i] = -vars.lapse * d1.phi[i] -
                            vars.phi * d1.lapse[i] + advec.Avec[i];
        FOR1(j)
        {
            total_rhs.Avec[i] +=
                -vars.lapse * gamma_LL[i][j] * vars.Evec[j] +
                vars.Avec[j] * d1.shift[j][i];
        }
    }

    // variable for term (D_i A_j - D_j A_i)
    // NB Christoffel symbols cancel and take care with
    // indices - the second index is the derivative index
    Tensor<2, data_t> diff_DA;
    FOR2(i, j) { diff_DA[i][j] = d1.Avec[j][i] - d1.Avec[i][j]; }


    // Xsquared = X^/mu X_/mu
    data_t Xsquared;
    Xsquared = -vars.phi * vars.phi;
    FOR2(i, j) { Xsquared += gamma_UU[i][j] * vars.Avec[j] * vars.Avec[i]; }

    // NB This is for E^i with indices up
    FOR1(i)
    {
        total_rhs.Evec[i] =
            vars.lapse * vars.K * vars.Evec[i] + advec.Evec[i];
        FOR1(j)
        {
            // dVdA = mu^2 ( 1 + 4 c4 (A^k A_k - phi^2))
            total_rhs.Evec[i] +=
                gamma_UU[i][j] * (vars.lapse * d1.Z[j] +
                                  vars.lapse * dVdA * vars.Avec[j]) -
                vars.Evec[j] * d1.shift[i][j];
        }

        FOR3(j, k, l)
        {
            total_rhs.Evec[i] +=
                gamma_UU[k][j] * gamma_UU[i][l] *
                (d1.lapse[k] * diff_DA[l][j] +
                 vars.lapse * (d2.Avec[j][l][k] - d2.Avec[l][j][k]));

            FOR1(m)
            {
                total_rhs.Evec[i] += -gamma_UU[k][j] * gamma_UU[i][l] *
                                     vars.lapse *
                                     (chris_phys[m][k][l] * diff_DA[m][j] +
                                      chris_phys[m][k][j] * diff_DA[l][m]);
            }
        }
    }

    // evolution equation for the damping term Z
    // dVdA = mu^2 ( 1 + 4 c4 A^k A_k - 12 c4 phi^2)
    // (ie the second part of the constraint, eqn 27)
    total_rhs.Z =
        vars.lapse * (dVdA * vars.phi - m_vector_damping * vars.Z) +
        advec.Z;

    FOR1(i)
    {
        total_rhs.Z += vars.lapse * d1.Evec[i][i];
        FOR1(j)
        {
            total_rhs.Z +=
                vars.lapse * chris_phys[i][i][j] * vars.Evec[j];
        }
    }

    //Extrinsic curvature
    Tensor<2, data_t> ExCurv;
    FOR2(i,j){
        ExCurv = (1./vars.chi)*(vars.A[i][j] + 1./3.*vars.h[i][j]*vars.K);
    }

    // C = 1 + 4 c4 A^k A_k - 12 c4 phi^2
    data_t C = 1.0 - 12.0 * c4 * vars.phi * vars.phi;
    FOR2(i, j)
    {
        C += 4.0 * c4 * gamma_UU[i][j] * vars.Avec[j] * vars.Avec[i];
    }

    total_rhs.phi = vars.lapse / C * (1.0 + 4.0 * c4 * Xsquared) *
                     (vars.K * vars.phi - DA_scalar)
                 // QUESTION: Should this be lapse * Z / C  or lapse * Z??
                 - vars.lapse * vars.Z / C + advec.phi;
    FOR1(i)
    {
        total_rhs.phi += 8.0 * c4 * vars.phi * vars.lapse / C *
                    (vars.Evec[i] * vars.Avec[i]);

        FOR1(j)
        {
            total_rhs.phi +=
                -gamma_UU[i][j] * vars.Avec[i] * d1.lapse[j] +
                8.0 * c4 * vars.phi * vars.lapse / C *
                    (2.0 * vars.Avec[i] * d1.phi[j] * gamma_UU[i][j]);

            FOR2(k, l)
            {
                total_rhs.phi +=
                    -8.0 * c4 * vars.lapse / C * gamma_UU[i][k] *
                        gamma_UU[j][l] * vars.Avec[i] * vars.Avec[j] *
                        DA[k][l] +
                    8.0 * c4 * vars.phi * vars.lapse / C *
                        (-ExCurv[i][j] * vars.Avec[k] *
                            vars.Avec[l] * gamma_UU[i][k] * gamma_UU[j][l]);
            }
        }
    }
    
    
    //################################################################################
#ifdef EQUATION_DEBUG_MODE
    DEBUG_HEADER;
    DEBUG_OUT3(V, dVdA, dVddA);
    DEBUG_OUT(vars.lapse);
    DEBUG_OUT(vars.K);
    DEBUG_OUT(total_rhs.phi);
    DEBUG_OUT(total_rhs.Z);
    DEBUG_OUT3(total_rhs.Avec[0], total_rhs.Avec[1], total_rhs.Avec[2]);
    DEBUG_OUT3(total_rhs.Evec[0], total_rhs.Evec[1], total_rhs.Evec[2]);
    DEBUG_OUT2(vars.phi, vars.Z);
    DEBUG_OUT3(vars.Avec[0], vars.Avec[1], vars.Avec[2]);
    DEBUG_OUT3(vars.Evec[0], vars.Evec[1], vars.Evec[2]);
    DEBUG_OUT(gnn);
    DEBUG_OUT4(gamma_UU[0][0], gamma_UU[0][1], gamma_UU[0][2], gamma_UU[1][1]);
    DEBUG_OUT2(gamma_UU[1][2], gamma_UU[2][2]);
    FOR2(i,j){
        DEBUG_OUT(d1.Evec[i][j]);
        //FOR1(k){
           // DEBUG_OUT(chris_phys[i][j][k]);
        //}
    }
    DEBUG_END;
#endif //EQUATION_DEBUG_MODE
    //################################################################################


}



#endif //PROCAFIELD_IMPL_H_INCLUDED


