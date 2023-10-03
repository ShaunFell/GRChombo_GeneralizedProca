#include "ProcaField.hpp"



template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> ProcaField<potential_t>::compute_emtensor(
        const vars_t<data_t> &vars, //the value of the variables
        const ADMVars::template Vars<data_t> &metric_vars, //value of metric variables
        const vars_t<Tensor<1,data_t>> &d1, //the 1st derivatives
        const Tensor<2, data_t> &gamma_UU, //the inverse metric
        const Tensor<3, data_t> &chris_phys_ULL //conformal christoffel symbols
    ) const {
        emtensor_t<data_t> out;

        //Compute potential and its derivatives
        data_t V; //value of potential
        data_t dVdA; //first derivative of Potential function w.r.t argument
        data_t dVddA; //second derivative...
        m_potential.compute_potential(V, dVdA, dVddA, vars, d1, metric_vars);

        //initialize components of EM tensor
        data_t rho_potential { 0 };
        Tensor<1, data_t> Si_potential;
        Tensor<2, data_t> Sij_potential;

        
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
            Enorm += metric_vars.gamma[i][j]*vars.Evec[i]*vars.Evec[j];
        };


        /////Components of EM tensor

        //Eulerian Energy 
        out.rho = Enorm + dVdA*vars.phi*vars.phi + 1/2*V;
        FOR4(i,j,k,l){
            out.rho += gamma_UU[k][i]*gamma_UU[l][j]*DA[i][j]*DA_antisym[k][l];
        };


        //Eulerian Momentum
        FOR1(i){
            out.Si[i] = 0; //zero initialize
            out.Si += vars.phi*dVdA*vars.Avec[i];
            
            FOR1(j){
                out.Si += 1./2. * (vars.Evec[j]*DA_antisym[i][j]);
            };
        };

        //Eulerian Stress
        FOR2(i,j){
            out.Sij[i][j] = 0; //zero initialize

            out.Sij[i][j] += 1./2. * ( 2*dVdA*vars.Avec[i]*vars.Avec[j] - metric_vars.gamma[i][j]*V + 2*metric_vars.gamma[i][j]*Enorm);

            FOR2(l,k){
                out.Sij[i][j] += 1./2. * (-metric_vars.gamma[i][l]*metric_vars.gamma[j][k]*vars.Evec[l]*vars.Evec[k] + gamma_UU[k][l]*DA_antisym[i][l]*DA_antisym[j][k]);

                FOR2(m,n){
                    out.Sij[i][j] += -metric_vars.gamma[i][j]*gamma_UU[m][l]*gamma_UU[n][k]*DA[m][n]*DA_antisym[l][k];
                };
            };
        };

        //Eulerian Stress scalar
        out.S = 0.0;
        FOR2(i,j){
            out.S += out.Sij[i][j]*gamma_UU[i][j];
        };

    };


template <class potential_t>
template <class data_t, 
            template <typename> class vars_t, 
            template <typename> class diff2_vars_t, 
            template <typename> class rhs_vars_t>
void ProcaField<potential_t>::matter_rhs(
        rhs_vars_t<data_t> &total_rhs, //RHS terms for all vars
        const vars_t<data_t> &vars, //the value fo the variables
        const ADMVars::template Vars<data_t> &metric_vars, //value of metric variables
        const vars_t<Tensor<1, data_t>> &d1, //value of 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //2nd derivs
        const vars_t<data_t> &advec //value of the beta^i d_i(var) terms
){


    //calculate full spatial christoffel symbols
    const auto gamma_UU = TensorAlgebra::compute_inverse(metric_vars.gamma);
    const auto chris_phys = TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    data_t V, dVdA, dVddA;
    m_potential.compute_potential(V, dVdA, dVddA, vars, d1, metric_vars);


    //evolution equations for spatial part of vector field (index down)
    FOR1(i){
        total_rhs.Avec[i] = - metric_vars.lapse*d1.phi[i] - vars.phi*metric_vars.d1_lapse[i] + advec.Avec[i];

        FOR1(j){
            total_rhs.Avec[i] += - metric_vars.lapse*metric_vars.gamma[i][j]*vars.Evec[j] + vars.Avec[j]*metric_vars.d1_shift[j][i];
        };
    };


    //evolution equations for Electric vector field (index up)
    FOR1(i){
        total_rhs.Evec[i] = metric_vars.lapse*metric_vars.K*vars.Evec[i] + advec.Evec[i];

        FOR1(j){
            total_rhs.Evec[i] += -metric_vars.lapse*gamma_UU[i][j]*d1.Z[j] + 2*metric_vars.lapse*dVdA*gamma_UU[i][j]*vars.Avec[j] - vars.Evec[j]*metric_vars.d1_shift[i][j];
        }

        FOR3(j,k,l){
            total_rhs.Evec[i] += -metric_vars.d1_lapse[j] * gamma_UU[j][k]*gamma_UU[i][l]*(d1.Avec[l][k] - d1.Avec[k][l]) - metric_vars.lapse*gamma_UU[j][k]*gamma_UU[i][l]*(d2.Avec[l][k][j] - d2.Avec[k][l][j]);

            FOR1(m){
                total_rhs.Evec[i] += -metric_vars.lapse*gamma_UU[j][k]*gamma_UU[i][l]*(chris_phys[m][j][l]*(d1.Avec[k][m] - d1.Avec[m][k]) + chris_phys[m][j][k]*(d1.Avec[m][l] - d1.Avec[l][m]));
            };
        };
    };


    //evolution equation for auxiliary constraint-damping scalar field Z
    total_rhs.Z = -2*metric_vars.lapse*dVdA*vars.phi - m_vector_damping*metric_vars.lapse*vars.Z + advec.Z;

    FOR1(i){
        total_rhs.Z += -metric_vars.lapse*d1.Evec[i][i];
        FOR1(j){
            total_rhs.Z += -metric_vars.lapse*chris_phys.ULL[i][i][j]*vars.Evec[j];
        }
    }


}


