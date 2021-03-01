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
 * This class contains the members and methods required for generating and storing the abiotic component of the model
 * including the (domain decomposed) spatial network, and explicitly modelled environmental distribiutions
 */

#ifndef LVMCM_TOPOGRAPHY_H
#define LVMCM_TOPOGRAPHY_H

#include <armadillo>
using namespace std;
using namespace arma;

class Topography {
public:
// members
    unsigned int subdomain = 0; // subdomain identifier
    uvec subdom_nodes; // subdomain node indices

    // parameters
    int no_nodes; // number of nodes (patches) in graph (metacommunity)
    int bisec; // number of bisections for domain decomposition
    int envVar; // number of environmental variables
    double var_e; // variance of the environmental distribution (implicit or explicit)
    double phi; // environmental autocorrelation length
    double T_int; // Intercept of linear temperature gradient
    string xMatFileName; // path to file name for imported network
    string scMatFileName; // path to file name for imported local scaling
    vec skVec; // vector of shape parameters defining environmental sensitivity

    // switches
    bool randGraph = true; // switch between random graph (1) and lattice (0)
    bool gabriel = true; // switch between gabriel (1) and complete graph (0)
    bool consArea_multiplicative; // switch for selecting additive of multiplicative for of perturbation in 'developed' area

    // matrix objects
    mat network; // x,y coords of nodes
    mat distMat; // euclidean distances between nodes
    mat sigEVec; // Eigenvectors of spatial covariance matrix
    mat sigEVal; // Eigenvalues of spatial covariance matrix
    mat adjMat; // spatial adjacency matrix
    uvec fVec; // indicator vector for domain decomposition
    uvec cFVec; // cumilatively count nodes in subdomains
    uvec intIF; // internal interface nodes
    uvec adjIF; // adjacent interface nodes
    mat envMat; // envVarxN matrix encoding the spatial distribution in enviromental variables

    vec consArea_bin; // binary vector for allocating conservation area
    vec scVec; // local interspecific interaction scaling
    vec scVec_prime; // local intraspecific interaction scaling

// methods
    // topo modelling
    void genNetwork(); // generate random network
    void genDistMat(); // generate distance matrix
    void genAdjMat(); // generate adjacency matrix
    void genDomainDecomp(mat netImprtd = {}); // domain decomposition algorithm
    void genEnvironment(); // Samples from envVar GRFs and stores values in environment matrix
    void genTempGrad(); // generate linear temperature gradient T(x) = T_int - sqrt(N)*x

// (default) constructor
    Topography() {}

// (default) deconstructor
    ~Topography() {}
};

#endif //LVMCM_TOPOGRAPHY_H

// Local Variables:
// c-file-style: "stroustrup"
// End:
