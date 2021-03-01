////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// The parallelizable Lotka-Volterra Metacommunity assembly Model (pLVMCM) ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Jacob Dinner O'Sullivan -- j.l.dinner@qmul.ac.uk | j.osullivan@zoho.com ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// A metacommunity assembly model parallelised via domain decomposition using MPI ////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
    Copyright (C) 2020  Jacob D. O'Sullivan, Axel G. Rossberg

    This file is part of pLVMCM

    pLVMCM is free software: you can redistribute it and/or modify
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

//////////////////////////////////// PRE-RELEASE VERSION, PLEASE DO NOT DISTRIBUTE /////////////////////////////////////

/*
 * This class contains the members and methods required interfacing with the SUNDIALS numerical ODE solver.
 */

#include "CommunityDynamics.h"

CommunityDynamics::CommunityDynamics(const CommunityDynamics & C) {

    bMat_p = C.bMat_p;
//    bMat_c = C.bMat_c;
    emMat_p = C.emMat_p;
//    emMat_c = C.emMat_c;
    rMat = C.rMat;
    efMat = C.efMat;
    uMat_p = C.uMat_p;
//    uMat_c = C.uMat_c;
    cMat = C.cMat;
    scVec = C.scVec;
    scVec_prime = C.scVec_prime;
//    aMat = C.aMat;
    dMat = C.dMat;
    dMat_n = C.dMat_n;
    dMat_m = C.dMat_m;
    rho = C.rho;
}

CommunityDynamics & CommunityDynamics::operator=(const CommunityDynamics &C) {

    if (this == &C) {
        return *this;
    } else {
        delete bMat_p;
//        delete bMat_c;
        delete emMat_p;
//        delete emMat_c;
        delete rMat;
        delete efMat;
        delete uMat_p;
//        delete uMat_c;
        delete cMat;
        delete scVec;
        delete scVec_prime;
//        delete aMat;
        delete dMat;
        delete dMat_n;
        delete dMat_m;
        delete rho;

        bMat_p = C.bMat_p;
//        bMat_c = C.bMat_c;
        emMat_p = C.emMat_p;
//        emMat_c = C.emMat_c;
        rMat = C.rMat;
        efMat = C.efMat;
        uMat_p = C.uMat_p;
//        uMat_c = C.uMat_c;
        cMat = C.cMat;
        scVec = C.scVec;
        scVec_prime = C.scVec_prime;
//        aMat = C.aMat;
        dMat = C.dMat;
        dMat_n = C.dMat_n;
        dMat_m = C.dMat_m;
        rho = C.rho;

        return *this;
    }
}

int CommunityDynamics::number_of_variables() const {
    // define number of state variables
    return (bMat_p->n_rows*bMat_p->n_cols);// + bMat_c->n_rows*bMat_c->n_cols);
}


void CommunityDynamics::read_state_from(const ODE_vector & state) {
    // converts 1D vector to N dimensional biomassMat
    int k = 0;
    for (int i=0; i<bMat_p->n_cols; i++) {
        for (int j=0; j<bMat_p->n_rows; j++) {
            (*bMat_p)(j,i) = state[k];
            k++;
        }
    }

//    for (int i=0; i<bMat_c->n_cols; i++) {
//        for (int j=0; j<bMat_c->n_rows; j++) {
//            (*bMat_c)(j,i) = state[k];
//            k++;
//        }
//    }
}

void CommunityDynamics::write_state_to(ODE_vector & state) const {
    // converts N dimensional biomassMat to 1D vector
    // note B_p and B_c are concatenated into single ODE_vector state
    int k = 0;
    for (int i=0; i<bMat_p->n_cols; i++) {
        for (int j=0; j<bMat_p->n_rows; j++) {
            state[k] = (*bMat_p)(j,i);
            k++;
        }
    }

//    for (int i=0; i<bMat_c->n_cols; i++) {
//        for (int j=0; j<bMat_c->n_rows; j++) {
//            state[k] = (*bMat_c)(j,i);
//            k++;
//        }
//    }
}

void CommunityDynamics::prepare_for_integration(){
    if (!dMat) {
        if (dMat_n) {
            // parallel dispersal operator
            dMat_n_sp = *dMat_n; // cast as sparse matrix
	    }
    }else{
	// non-parallel dispersal operator
        dMat_sp = *dMat; // cast as sparse matrix
    }
    // No need to set size of prodDisp, consDisp.  Armadillo
    // automatically gets this right.
}

void CommunityDynamics::cleanup_after_integration(){
    // We might want to clean up dMat_n_sp, dMat_sp, prodDisp,
    // consDisp here, but this will happen anyway once the object is
    // removed, and if it is not cleaning them up might be
    // inefficient.  So we don't.
}

int CommunityDynamics::dynamics(ODE_vector const & state, ODE_vector & time_derivative) {

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
        // updates to sppPool matrices - no return value since memory is shared

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Initialize temporary matrix objects ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const mat B_p(&state[0], rMat->n_rows, bMat_p->n_cols, false, true);
    const mat B_c(&state[rMat->n_rows*bMat_p->n_cols], bMat_p->n_rows - rMat->n_rows, bMat_p->n_cols, false, true);

    mat dBdt_p(&time_derivative[0], rMat->n_rows, bMat_p->n_cols, false, true);
    mat dBdt_c(&time_derivative[rMat->n_rows*bMat_p->n_cols], bMat_p->n_rows - rMat->n_rows, bMat_p->n_cols, false, true);

    // growth and dispersal objects used to construct model dynamics in _modular_ way
    mat& prodGrowth(dBdt_p);  // alias, saves one copy operation
    mat& consGrowth(dBdt_c);  // alias, saves one copy operation

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Modular construction of model objects ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // producer internal growth rate
    prodGrowth = *rMat;

    if (testing) {
        if (S_p > 0) {
            prodGrowth.rows(0,S_p-1).zeros(); // resident growth set to zero for invader testing
        }
    }

    if (efMat) { // environmental fluctuations added to growth rate matrix
        prodGrowth += *efMat;
    }

    // producer density dependence
    if (cMat) { // intra-/inter-specific
        if (scVec) { // scaling of local interaction terms - assume cMat hollow
            mat intersp, intrasp;
            intersp = (cMat->submat(0,0,rMat->n_rows-1, rMat->n_rows-1)) * B_p;
            intersp.each_row() %= conv_to<rowvec>::from(*scVec);
            intrasp = B_p; // assume self interaction term = 1
            intrasp.each_row() %= conv_to<rowvec>::from(*scVec_prime);
            prodGrowth -= (intersp + intrasp); // full model
        } else {
            prodGrowth -= (cMat->submat(0,0,rMat->n_rows-1, rMat->n_rows-1)) * B_p;
        }
    } else {
        prodGrowth -= B_p; // intra-specific only
    }

    // trophic interactions -- spatial scaling???
    if (B_c.n_rows != 0) {
        prodGrowth -= (cMat->submat(0,rMat->n_rows,rMat->n_rows-1,cMat->n_rows-1)) * B_c;
        consGrowth = *rho * (cMat->submat(rMat->n_rows,0,cMat->n_rows-1,rMat->n_rows-1) * B_p - 1);
    }

    if (dMat_n) {
        if (B_p.n_rows != 0) {
            if (!emMat_p) {
                prodDisp = B_p * dMat_n_sp + (uMat_p->rows(0, B_p.n_rows-1)); // uMat: time independent _fixed_ unknowns SPARSE
            } else {
                mat emMat_p_N(B_p.n_rows, B_p.n_cols);
                for (int i=0; i<B_p.n_rows; i++) {
                    emMat_p_N.row(i).fill(emMat_p->at(i));
                }
                prodDisp = (emMat_p_N % B_p) * dMat_n_sp + (uMat_p->rows(0, B_p.n_rows-1)); // uMat: time independent _fixed_ unknowns SPARSE
            }
        }

        if (B_c.n_rows != 0) {
            if (!emMat_p) {
                consDisp = B_c * dMat_n_sp;
                mat TEMP = (uMat_p->rows(B_p.n_rows, uMat_p->n_rows-1)); // uMat: time independent _fixed_ unknowns SPARSE
                bool printTEMP=false;
                if (printTEMP) {
                    B_p.print("\nB_p");
                    B_c.print("\nB_c");
                    cout << "\nU_p\n" << uMat_p->rows(0, B_p.n_rows-1) << endl;
                    cout << "\nU_c\n" << uMat_p->rows(B_p.n_rows, uMat_p->n_rows-1) << endl;
                }
                consDisp += TEMP;
            } else {
                mat emMat_c_N(B_c.n_rows, B_c.n_cols);
                for (int i=0; i<B_c.n_rows; i++) {
                    emMat_c_N.row(i).fill(emMat_p->at(B_p.n_rows + i));
                }
                consDisp = (emMat_c_N % B_c) * dMat_n_sp + (uMat_p->rows(B_p.n_rows, uMat_p->n_rows-1)); // uMat: time independent _fixed_ unknowns SPARSE
            }
        }
    }

    if (dMat) {
        if (bMat_p) {
            if (!emMat_p) {
                prodDisp = B_p * dMat_sp;
            } else {
                mat emMat_p_N(B_p.n_rows, B_p.n_cols);
                for (int i=0; i<B_p.n_rows; i++) {
                    emMat_p_N.row(i).fill(emMat_p->at(i));
                }
                prodDisp = (emMat_p_N % B_p) * dMat_sp;
            }
        }
        if (B_c.n_rows != 0) {
            if (!emMat_p) {
                consDisp = B_c * dMat_sp;
            } else {
                mat emMat_c_N(B_c.n_rows, B_c.n_cols);
                for (int i = 0; i < B_c.n_rows; i++) {
                    emMat_c_N.row(i).fill(emMat_p->at(B_p.n_rows + i));
                }
            }
        }
    }

    prodGrowth %= B_p;
    consGrowth %= B_c;

//    dBdt_p = prodGrowth; // don't need this, because the two are just aliases
//    dBdt_c = consGrowth; // don't need this, because the two are just aliases

    if (prodDisp.n_rows > 0) {
        dBdt_p += prodDisp;
    }

    if (consDisp.n_rows > 0) {
        dBdt_c += consDisp;
    }

    bool print_mat=false;
    if (print_mat) {
        B_p.print("\n\n\nBcd1");
        uMat_p->print("\n\n\nUcd");
        B_c.print("\nB_c");
        rMat->print("\nR");
        cMat->print("\nA");
        scVec->print("\nSc");
        prodDisp.print("\nprodDisp");
        consDisp.print("\nconsDisp");
        prodGrowth.print("\nprodGrowth");
        consGrowth.print("\nconsGrowth");
    }

    // set negative elements to absolute - accounts for small numerical errors but produces large 'oscillations' in pathological case
    uvec negativeB_p = find(B_p < 0);
    dBdt_p(negativeB_p) = - B_p(negativeB_p);
    uvec negativeB_c = find(B_c < 0);
    dBdt_c(negativeB_c) = - B_c(negativeB_c);

    // remove nan arising due to infinitesimal biomasses
    uvec nonfiniteB_p = find_nonfinite(dBdt_p);
    dBdt_p(nonfiniteB_p).fill(0);
    uvec nonfiniteB_c = find_nonfinite(dBdt_c);
    dBdt_c(nonfiniteB_c).fill(0);
    if (nonfiniteB_p.n_rows > 0) {
        cout << "Warning: non finite dBdt_c detected" << endl;
    }
    if (nonfiniteB_c.n_rows > 0) {
        cout << "Warning: non finite dBdt_c detected" << endl;
    }

    // remove negative blow up arising due to infinitesimal immigration rates from numerical zeros
    uvec blowupB_p = find(abs(dBdt_p)>1e10);
    dBdt_p(blowupB_p).fill(0);
    uvec blowupB_c = find(abs(dBdt_c)>1e10);
    dBdt_c(blowupB_c).fill(0);
//    if (blowupB_p.n_rows > 0) {
//        cout << "Warning: dBdt_p blow-up detected" << endl;
//    }
//    if (blowupB_c.n_rows > 0) {
//        cout << "Warning: dBdt_c blow-up detected" << endl;
//    }

    return 0;
}

// Local Variables:
// c-file-style: "stroustrup"
// End:
