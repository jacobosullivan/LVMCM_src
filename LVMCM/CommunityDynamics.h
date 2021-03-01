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

#ifndef LVMCM_COMMUNITYDYNAMICS_H
#define LVMCM_COMMUNITYDYNAMICS_H

#include <iostream>
#include "ODE.h"
#include "error.h"
#include <armadillo>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include "Metacommunity.h"

using namespace std;
using namespace arma;

//#define SPARSE_IMAT // switch to sparse format

class CommunityDynamics : public ODE_dynamical_object {
private:
    // Matrices to be pre-allocated before integration:
    sp_mat dMat_n_sp;
    sp_mat dMat_sp;
    mat prodDisp;
    mat consDisp;
public:
// members
    bool testing = false; // switch for invader testing (intrasp, disp switched off)
    int S_p, S_c; // resident species diversity
    // matrix objects
    mat *bMat_p {}; // producer biomass matrix
//    mat *bMat_c {}; // consumer biomass matrix
    mat *emMat_p {}; // producer biomass matrix
//    mat *emMat_c {}; // consumer biomass matrix
    const mat *uMat_p {}; // producer interface state
//    const mat *uMat_c {}; // consumer interface state
#ifdef SPARSE_IMAT
    sp_mat *cMat; // producer interaction matrix (sparse)
#else
    mat *cMat {}; // producer interaction matrix
#endif
    vec *scVec {}; // local interspecific interaction scaling
    vec *scVec_prime {}; // local intraspecific interaction scaling
//    mat *aMat {}; // trophic interaction matrix
    mat *rMat {}; // producer growth rate matrix
    mat *efMat {}; // environmental fluctuation matrix
    mat *dMat {}; // regional dispersal operator
    mat *dMat_n {}; // subdomain dispersal operator
    mat *dMat_m {}; // inter-subdomain dispersal operator
    double *rho {}; // consumer respiration rate

// methods
    virtual int dynamics(ODE_vector const & state, ODE_vector & time_derivative); // CVode dynamical object
    virtual void write_state_to(ODE_vector & state) const; // reads state variables into ODE_vector 'state'
    virtual void read_state_from(const ODE_vector & state); // writes state to shared memory object bMat_p/c
    virtual int number_of_variables() const; // computes number of state variables
    void prepare_for_integration(); // pre-compute some relevant objects
    void cleanup_after_integration();
    // (default) constructor
    CommunityDynamics() {}
    // (default) deconstructor
    ~CommunityDynamics() {}
    // copy constructor
    CommunityDynamics(const CommunityDynamics &C);
    // assignment operator
    CommunityDynamics & operator=(const CommunityDynamics &C);

};

#endif //LVMCM_COMMUNITYDYNAMICS_H

// Local Variables:
// c-file-style: "stroustrup"
// End:
