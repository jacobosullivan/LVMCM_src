////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// The parallelizable Lotka-Volterra Metacommunity assembly Model (LVMCM) ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Jacob Dinner O'Sullivan -- j.l.dinner@qmul.ac.uk | j.osullivan@zoho.com ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// A metacommunity assembly model parallelised via domain decomposition using MPI ////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
    Copyright (C) 2022  Jacob D. O'Sullivan, Axel G. Rossberg

    This file is part of LVMCM

    LVMCM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/



/*
 * This class contains the members and methods required interfacing with the SUNDIALS numerical ODE solver.
 */

#include "CommunityDynamics.h"
#include "Metacommunity.h"
#include "LVMCM_rng.h"

int CommunityDynamics::number_of_variables() const {
    // define number of state variables
    return (xMat->n_rows*xMat->n_cols);
}

int CommunityDynamics::private_number_of_derived_quantities() {
    // define number of state derivatives
    ////FATAL_ERROR("Should be number of derived quantities!");
    return (xMat->n_rows*xMat->n_cols);
}

void CommunityDynamics::read_state_from(const ODE_vector & state) {
    // converts 1D vector to N dimensional biomassMat
    int k = 0;
    for (int j=0; j<xMat->n_cols; j++) {
        for (int i=0; i<xMat->n_rows; i++) {
            (*xMat)(i,j) = state[k];
            k++;
        }
    }
}

void CommunityDynamics::write_state_to(ODE_vector & state) const {
    // converts N dimensional biomassMat to 1D vector
    // note B_p and B_c are concatenated into single ODE_vector state
    int k = 0;
    for (int j=0; j<xMat->n_cols; j++) {
        for (int i=0; i<xMat->n_rows; i++) {
            state[k] = (*xMat)(i,j);
            k++;
        }
    }
}

void CommunityDynamics::prepare_for_integration(){
    dMat_sp = *dMat; // cast as sparse matrix
    if(number_of_derived_quantities()==0){
	// we use the two-parameter form of dynamcis(..)
	// Not nice! Better to defined .resize for ODE_vector!
	derived_quantity_dump=ODE_vector(private_number_of_derived_quantities());
    }
}

void CommunityDynamics::cleanup_after_integration(){
    // We might want to clean up dMat_n_sp, dMat_sp, prodDisp,
    // consDisp here, but this will happen anyway once the object is
    // removed, and if it is not cleaning them up might be
    // inefficient.  So we don't.
}

int CommunityDynamics::number_of_root_functions() const{
//    cout << "NO ROOT FUNCTIONS = " << 2*(indices_S->n_rows) + indices_DP->n_rows << endl;
    if (g_form_of_dynamics) {
        return(2*(indices_S->n_rows) + indices_DP->n_rows);
    } else {
        return(1);
    }

}

void CommunityDynamics::compute_intrinsic_growth_rates(const mat B, mat & Gt) const{
    ASSERT(Gt.n_rows == cMat->n_rows && Gt.n_cols == xMat->n_cols);
    ASSERT(B.n_rows == cMat->n_rows && B.n_cols == xMat->n_cols);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Modular construction of model objects ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // producer intrinsic growth rates
    // Set: Gt_p = R
    Gt.rows(0, rMat->n_rows-1) = *rMat;

    if (efMat) { // environmental fluctuations added to growth rate matrix
        Gt.rows(0, rMat->n_rows-1) += *efMat;
    }

    // producer density dependence
    // Set: Gt = (R - C x B_p)
    // Set: Gt = (R - C x B_p)
    if (cMat) { // intra-/inter-specific
        if (scVec) { // scaling of local interaction terms - assume cMat hollow
            mat intersp = (cMat->submat(0, 0, rMat->n_rows - 1, rMat->n_rows - 1)) * B.rows(0, rMat->n_rows-1);
            intersp.each_row() %= conv_to<rowvec>::from(*scVec);
	    mat intrasp = B.rows(0, rMat->n_rows-1); // assume self interaction term = 1
            intrasp.each_row() %= conv_to<rowvec>::from(*scVec_prime);
            Gt.rows(0, rMat->n_rows-1) -= (intersp + intrasp); // full model
        } else {
            Gt.rows(0, rMat->n_rows-1) -= (cMat->submat(0, 0, rMat->n_rows - 1, rMat->n_rows - 1)) * B.rows(0, rMat->n_rows-1);
        }
    } else {
        Gt.rows(0, rMat->n_rows-1) -= B.rows(0, rMat->n_rows-1); // intra-specific only
    }

    // trophic interactions -- spatial scaling???
    // Set: Gt_c = (R - C x B_p - A x B_c)
    // Set: Gt_c = rho * (AB_p - 1)
    if (B.n_rows > rMat->n_rows) {
        Gt.rows(0, rMat->n_rows-1) -= (cMat->submat(0, rMat->n_rows, rMat->n_rows - 1, cMat->n_rows - 1)) * B.rows(rMat->n_rows, B.n_rows-1);
        Gt.rows(rMat->n_rows, Gt.n_rows-1) = *rho * (cMat->submat(rMat->n_rows, 0, cMat->n_rows - 1, rMat->n_rows - 1) * B.rows(0, rMat->n_rows-1) - 1);
    }

    return;
}


template <class Container, class UnaryFunction>
inline Container apply_all(Container c, UnaryFunction f)
{
    std::transform(c.begin(), c.end(), c.begin(), f);
    return c;
}

int CommunityDynamics::root_functions(const ODE_vector & state, const ODE_vector & derived, ODE_vector & gout) {

    if (g_form_of_dynamics) {
        typedef std::vector<double> stdvec;
        for(int i=0; i<indices_S->n_rows; i++) { // pass state of Poisson clocks to root finder
            gout[i]=state[conv_to< stdvec >::from(*indices_S)[i]]+1;
        }

        mat X(&state[0], cMat->n_rows, xMat->n_cols, false, true);

#ifdef IDA
        mat Gt(&derived[0], cMat->n_rows, xMat->n_cols, false, true);
#endif

#ifndef IDA
        mat Gt(X.n_rows, X.n_cols);
        {
            mat B(X);
            B.elem(*indices_S).fill(0);
            compute_intrinsic_growth_rates(B,Gt);
        }
#endif

        vec gout_alias(&gout[0], gout.size(), false, true); // make use of armadillo non-contiguous memory access

        gout_alias.elem(indices_S->n_rows + linspace<uvec>(0,indices_S->n_rows-1,indices_S->n_rows)) = Gt.elem(*indices_S);
        int nRow = cMat->n_rows;
        // Transition P->S only if growth rate is still positive when 1 ind is there:
        gout_alias.elem(2*indices_S->n_rows + linspace<uvec>(0,indices_DP->n_rows-1,indices_DP->n_rows)) = Gt.elem(*indices_DP) -
                                                                                                           *bodymass *
                                                                                                           vec(cMat->diag())(apply_all(*indices_DP,[nRow](int i){return i % nRow;}));
        // set G_ix, (i,x) \in {D} to -1 to avoid tripping root finder
        gout_alias.elem(2*indices_S->n_rows + find(X.elem(*indices_DP) > 0.9*(*bodymass))).fill(-1);
        // set G_ix, (i,x) \in {D,P} to nan if G_ix>0 to avoid tripping root finder
        // these are on the way to P from above
//    gout_alias.elem(2*indices_S->n_rows + find(Gt.elem(*indices_DP) > 0)).fill(datum::nan);

        if (0) {
            X.print("X_root_fun");
            Gt.print("Gt_root_fun");
            cout << "gout_alias " << gout_alias(min((int) gout_alias.size()-1,56)) << endl;
        }
    } else {
        gout[0] = -1; // excess work to avoid redundancy with PSD root finding functionality
    }
    return 0;
}

int CommunityDynamics::set_root_directions(int * direction) const{
    if (g_form_of_dynamics) {
        int nrtfn=number_of_root_functions();
        int i=0;
        for(; i<indices_S->n_rows; i++) { // Poisson clocks for S states, going up
            direction[i] = +1;
        }
        for(; i < 2 * indices_S->n_rows; i++) { // Linear growth rates for S states, going down
            direction[i] = -1;
        }
        for(; i < nrtfn; i++) { // Linear growth rates for P states, going up
            direction[i] = +1;
        }
    }

    return 0;
}

// Without tracking of local growth rates
int CommunityDynamics::dynamics(ODE_vector const & state, ODE_vector & time_derivative) {
    return dynamics(state,time_derivative,derived_quantity_dump);
}

// With tracking of local growth rates
int CommunityDynamics::dynamics(ODE_vector const & state, ODE_vector & time_derivative, ODE_vector & derived_quantities) {

    // summary:
    // numerically approximate (serial/parallel) the solution to the LVMCM (competitive/bipartite)
    // dynamics equations are constructed in a modular way for flexibility
    // interfaces with SUNDIALs ODE solver

    // required members:
    // members of CommunityDynamics object are pointers to matrix objects stored in a object of class Species

    // arguments:
    // state - ODE_vector representing the current state (biomass) of the metacommunity
    // time_derivative - ODE_vector representing the dynamics of the system constructed in a modular way by successive matrix operations

    // output:
    // updates to spp matrices - no return value since memory is shared

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Initialize temporary matrix objects ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const mat X(&state[0], cMat->n_rows, xMat->n_cols, false, true);
    mat dXdt(&time_derivative[0], cMat->n_rows, xMat->n_cols, false, true);
    mat Gt(&derived_quantities[0], cMat->n_rows, xMat->n_cols, false, true);
    mat massEffect(dXdt.n_rows, dXdt.n_cols); // need to use BxD for dB/dt and dPsi/dt therefore need a copy external to dXdt

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Modular construction of model objects ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Def B: state matrix for D,P compartments
    mat B(X.n_rows, X.n_cols, fill::zeros);
    B.elem(*indices_DP) = X.elem(*indices_DP); // fill B using indexing

    compute_intrinsic_growth_rates(B, Gt);
    dXdt = Gt; // keep local growth rates for compartment transitions

    // Set: dXdt = G.B
    dXdt.elem(*indices_DP) %= B.elem(*indices_DP);

    // Set: dXdt = G.B + BD (both consumer and producer dispersal handled simultaneously)
    if (!emMat) {
        massEffect = B * dMat_sp;
        dXdt += massEffect;
    } else {
        mat emMat_N(B.n_rows, B.n_cols);
        for (int i = 0; i < B.n_rows; i++) {
            emMat_N.row(i).fill(emMat->at(i));
        }
        dXdt += massEffect;
        massEffect = (emMat_N % B) * dMat_sp;
    }

    // Set: dPsi/dt = W.H.G / m(G + mu), dPsi_ix/dt = 0 for (i,x) \in {D,P}
    vec dPsi(indices_S->n_rows);
    dPsi = *bodymass_inv * massEffect.elem(*indices_S) % Gt.elem(*indices_S) % pow(Gt.elem(*indices_S) + *mu, -1);

    // Make sure dPsi is always non-negative
    dPsi(find(Gt.elem(*indices_S) < 0)).fill(0);

    // Set: dX_ix/dt = dPsi_ix/dt, for (i,x) \in {S}
    if (indices_S) {
        dXdt.elem(*indices_S) = dPsi;
    }

    // set negative elements to absolute - accounts for small numerical errors but produces large 'oscillations' in pathological case
    const uvec negativeB = find(B < 0);
    dXdt(negativeB) = - B(negativeB);

    // remove nan arising due to infinitesimal biomasses
    const uvec nonfiniteX = find_nonfinite(dXdt);
    if (nonfiniteX.n_rows > 0) {
//        cout << "\nstate[0] " << state[0] << endl;
//        X.print("X_cd");
//        B.print("B_cd");
//        dXdt.print("dXdt");
        cout << "Warning: non finite dXdt detected" << endl;
	return -1;
    }

    // remove negative blow up arising due to infinitesimal immigration rates from numerical zeros
    uvec blowupX = find(abs(dXdt)>1e10);
    if (blowupX.n_rows > 0) {
        cout << "Warning: dXdt blow-up detected" << endl;
        blowupX.print("blowupX");
    }
    dXdt(blowupX).fill(0);

    bool print_mat=false;
    if (print_mat) {
        X.print("X_cd");
//        B.print("B_cd");
        dXdt.print("dXdt_cd");
        dPsi.print("dPsi_cd");
        cout << "bodymass_inv " << *bodymass_inv << endl;
        massEffect.elem(*indices_S).print("massEffect");
        Gt.elem(*indices_S).print("Gt");
        pow(Gt.elem(*indices_S) + *mu, -1).print("Gt inv");
    }

    return 0;
}

#ifdef IDA
int CommunityDynamics::precondition_IDA(ODE_vector const & state,
					ODE_vector const & derived_quantities,
					ODE_vector const & in,
					ODE_vector & out,
					realtype alpha) const{
    
    const mat X(&state[0], cMat->n_rows, xMat->n_cols, false, true);
    const mat Gt(&derived_quantities[0], cMat->n_rows, xMat->n_cols, false, true);
    const mat in_mat(&in[0], cMat->n_rows, xMat->n_cols, false, true);
    mat out_mat(&out[0], cMat->n_rows, xMat->n_cols, false, true);
    mat P_diag(X-Gt); // Assumes interaction matrix has 1 on diagonal!
    P_diag.elem(*indices_S).fill(0);
    P_diag += alpha;
    out_mat = in_mat / P_diag; // Elementwise division [N.B.: Jacobian is (N * S) x (N * S) matrix]
    return 0;
}

int CommunityDynamics::set_inequality_constraints(ODE_vector & constraints) const {
    mat constraint_mat(&constraints[0], cMat->n_rows, xMat->n_cols, false, true);
    constraint_mat.elem(*indices_S).fill(0); // any value
    constraint_mat.elem(*indices_DP).fill(1); // non-negative
    return 0;
};

#endif

inline unsigned int* oneDto2Dind(const int oneDindex, const int nRows) {
    // convert 1D columnwise indices to 2D row-column indices
    static unsigned int res [2];
    res[0] = oneDindex % nRows;
    res[1] = oneDindex / nRows;
    return res;
}

void CommunityDynamics::react_to_roots(ODE_state::root_indices_t &root_indices) {

    bool REPORT_DETECTION = true;

    if( !root_indices.empty() ) {
//        cout << "REACT TO ROOTS TRIGGERED" << endl;
        uvec indices_DP_tmp = *indices_DP; // copy needed to avoid errors in case root_indices.size() > 1
        uvec indices_S_tmp = *indices_S;

        uvec root_indices_vec = conv_to<uvec>::from(abs(conv_to<ivec>::from(root_indices))) - 1; // needed to make use of aramdillo non contiguous memory access
        for (int i=root_indices_vec.n_rows-1; i>=0; i--) {
//	    REPORT(sign(root_indices[i]));
//	    REPORT(abs(root_indices[i])-1);
	    
            if (abs(root_indices[i])-1 < indices_S->n_rows) {
                // detected change in state of Poisson clock, S->D transition
                if (REPORT_DETECTION) {
                    cout << "S->D detected" << endl;
                }
                if (root_indices[i] < 0) {
                    cout << "ERROR: Poisson clock approached zero from above, index " << abs(root_indices[i])-1 << endl;
                }

                // Compute 1 and 2D indices of full model objects (i.e. not subsets stored in indices_Y)
                unsigned int true_1Dindex = (*indices_S)(root_indices_vec(i));
                unsigned int* true_2Dindex = oneDto2Dind(true_1Dindex, xMat->n_rows); //points to static within oneDto2Dind(..)

                if (!g_block_transitions) {
                    // Compute establisment probability of ith species and update X_ix
                    double growthRate_i = (*rMat)(true_1Dindex);
                    vec density_dependence = cMat->row(true_2Dindex[0]) * xMat->col(true_2Dindex[1]);
                    if (density_dependence.n_rows > 1) {
                        cout << "Non-scalar density dependence!" << endl;
                    }
                    growthRate_i -= density_dependence(0);
                    double probEstablish = growthRate_i / (growthRate_i - *mu);
                    (*xMat)(true_1Dindex) = *bodymass / probEstablish;

                    // update indices (temp vectors)
                    indices_S_tmp.shed_row(abs(root_indices[i])-1);
                    indices_DP_tmp.resize(indices_DP_tmp.n_rows+1);
                    indices_DP_tmp(indices_DP_tmp.n_rows-1) = true_1Dindex;

                } else {
                    (*xMat)(true_1Dindex) = -2;
                }

            } else if (abs(root_indices[i])-1 >= indices_S->n_rows && abs(root_indices[i])-1 < 2*indices_S->n_rows) {
                // detected change in growth rate for (i,x) \in {S}, S->P transition
                if (REPORT_DETECTION) {
                    cout << "S->P detected" << endl;
                }

                // Compute 1D index of full model objects (i.e. not subsets stored in indices_Y)
                unsigned int true_1Dindex = (*indices_S)(root_indices_vec(i) - indices_S->n_rows);

                if (root_indices[i] > 0) {
                    cout << root_indices_vec(i) << endl;
                    cout << "\nERROR: Growth rate of S species approached zero from below, true index " << true_1Dindex << endl;
                    REPORT(current_time);
                    continue;
                }

                if (!g_block_transitions) {
                    // Set X_ix = 0 (i.e. delete Poisson clock)
                    (*xMat)(true_1Dindex) = 0;

                    // update indices (temp vectors)
                    indices_S_tmp.shed_row(abs(root_indices[i]) - indices_S->n_rows - 1);
                    indices_DP_tmp.resize(indices_DP_tmp.n_rows + 1);
                    indices_DP_tmp(indices_DP_tmp.n_rows - 1) = true_1Dindex;
                } else {
                    (*xMat)(true_1Dindex) = -2;
                }

            } else if (abs(root_indices[i])-1 >= 2*indices_S->n_rows) {
                // detected change in growth rate for (i,x) \in {P}, P->S/D transition
                if (REPORT_DETECTION) {
                    cout << "P->D/S detected... ";
                }
                if (root_indices[i] < 0) {
                    cout << "ERROR: Growth rate of P species approached zero from above, index " << abs(root_indices[i])-1 << endl;
		            continue;
                }

                // Compute 1 and 2D indices of full model objects (i.e. not subsets stored in indices_Y)
                unsigned int true_1Dindex = (*indices_DP)(root_indices_vec(i) - 2*indices_S->n_rows);
                unsigned int* true_2Dindex = oneDto2Dind(true_1Dindex, xMat->n_rows);

                double p_extant;

                if (!g_block_transitions) {
                    // Compute probability that ith species is extant in xth site
                    // Def B_i: state vector for i \in D,P compartments
                    rowvec B = xMat->row(true_2Dindex[0]);
                    B(find(B < 0)).zeros();
                    vec massEffect_i = B * dMat->col(true_2Dindex[1]);

                    // REPORT(B.row(true_2Dindex[0]));
                    // REPORT(dMat->col(true_2Dindex[1]));
                    if (massEffect_i.n_rows > 1) {
                        cout << "ERROR: Non-scalar mass effect!" << endl;
                    }
                    p_extant = 1 - pow(2, -massEffect_i(0) * (*bodymass_inv) / *mu);
                    // REPORT(massEffect_i(0));
                    // REPORT(-massEffect_i(0) * (*bodymass_inv)  / *mu);
                    // REPORT(p_extant);
                } else {
                    p_extant = -1.0;
                }
                vec z_PD(1);
                z_PD.randu();

                if (z_PD(0) < p_extant) {
                    if (REPORT_DETECTION) {
                        cout << "D selected" << endl;
                    }
                    // Update X_ix
                    (*xMat)(true_1Dindex) = max(*bodymass, min(1.0, abs((*xMat)(true_1Dindex)) / p_extant));
                } else {
                    if (REPORT_DETECTION) {
                        cout << "S selected" << endl;
                    }

                    // Update X_ix
                    vec poissonClock;
                    poissonClock.randu(1); // inefficient but convenient
                    poissonClock = log(poissonClock) - 1; // -1 give buffer against very small negative biomasses in P compartment
                    (*xMat)(true_1Dindex) = poissonClock(0); // Sample Poisson clock; equivalent to z~-1*Exp(1) - 1

                    // update indices (temp vectors)
                    indices_DP_tmp.shed_row(abs(root_indices[i])-2*indices_S->n_rows-1);
                    indices_S_tmp.resize(indices_S_tmp.n_rows+1);
                    indices_S_tmp(indices_S_tmp.n_rows-1) = true_1Dindex;
                }
            }
        }

        if (indices_DP_tmp.n_rows != indices_DP->n_rows) {
            *indices_DP = sort(indices_DP_tmp);
        }

        if (indices_S_tmp.n_rows != indices_S->n_rows) {
            *indices_S = sort(indices_S_tmp);
        }
    }
}

// Local Variables:
// c-file-style: "stroustrup"
// End:
