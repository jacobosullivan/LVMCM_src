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

#ifndef LVMCM_COMMUNITYDYNAMICS_H
#define LVMCM_COMMUNITYDYNAMICS_H

#include <iostream>
#include "ODE.h"
#include "error.h"
#include <armadillo>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include "Metacommunity.h"

//#define IDA // comment for CVode solver
//#define SPARSE_IMAT // switch to sparse format

using namespace std;
using namespace arma;

class CommunityDynamics : public ODE_dynamical_object {
private:
    // Matrices to be pre-allocated before integration:
    sp_mat dMat_sp;
public:
// members
    int S_p, S_c; // resident species diversity

    // matrix objects
    mat *xMat {}; // producer biomass matrix
    mat *emMat {}; // producer biomass matrix
#ifdef SPARSE_IMAT
    sp_mat *cMat; // producer interaction matrix (sparse)
#else
    mat *cMat {}; // producer interaction matrix
#endif
    vec *scVec {}; // local interspecific interaction scaling
    vec *scVec_prime {}; // local intraspecific interaction scaling
    mat *rMat {}; // producer growth rate matrix
    mat *efMat {}; // environmental fluctuation matrix
    mat *dMat {}; // regional dispersal operator
    mat *psdMat_pd {}; // SxN indicator matrix with indices of D,P and S components
    ODE_vector derived_quantity_dump; // write derived quanities here if not needed.
    uvec *indices_DP {}; // indices of species in {D,P}
    uvec *indices_S {}; // indices of species in {S}
    double *rho {}; // consumer respiration rate
    double *bodymass {};
    double *bodymass_inv {}; // inverse body mass
    double *mu {}; // mortality

// methods
    virtual int dynamics(ODE_vector const & state, ODE_vector & time_derivative); // IDA dynamical object
    virtual int dynamics(ODE_vector const & state, ODE_vector & time_derivative, ODE_vector & derived_quantities); // IDA dynamical object
    virtual void write_state_to(ODE_vector & state) const; // reads state variables into ODE_vector 'state'
    virtual void read_state_from(const ODE_vector & state); // writes state to shared memory object xMat
    virtual int number_of_variables() const; // computes number of state variables
private:
    int private_number_of_derived_quantities(); // this can be nonzero even if the next is zero
public:
#ifdef IDA
    virtual int number_of_derived_quantities() {
        return private_number_of_derived_quantities();
    }; // computes number of derived variables
    virtual int precondition_IDA(ODE_vector const & state,
				 ODE_vector const & derived,
				 ODE_vector const & in,
				 ODE_vector & out,
				 realtype alpha) const;
    // Current preconditioner makes code slighly faster but does not
    // help solving issues with "Nonlinear convergence failure rate"
//    virtual bool has_preconditioner(){return true;};
    virtual int set_inequality_constraints(ODE_vector & constraints) const;
#endif
    void prepare_for_integration(); // pre-compute some relevant objects
    void cleanup_after_integration();

    void compute_intrinsic_growth_rates(const mat B, mat & Gt) const;
    
    // root finding machinery:
    virtual int number_of_root_functions() const;
    virtual int set_root_directions(int * direction) const;
    virtual int root_functions(ODE_vector const & state, ODE_vector const & derived, ODE_vector & gout);
    void react_to_roots(ODE_state::root_indices_t & root_indices);

    // (default) constructor - required for classes with pointer members
    CommunityDynamics() {}
    // (default) deconstructor - required for classes with pointer members
    ~CommunityDynamics() {}
};

#endif //LVMCM_COMMUNITYDYNAMICS_H

// Local Variables:
// c-file-style: "stroustrup"
// End:
