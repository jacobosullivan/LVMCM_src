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
 * This class contains the members and methods required for generating and storing the biotic component of the model,
 * the species pool. Methods include invasion (sampling) of new species and extinction (removal) of exluded species.
 */

#ifndef SETC_SPECIES_H
#define SETC_SPECIES_H

#include "Topography.h"
#include <armadillo>

using namespace std;
using namespace arma;

class Species: public Topography {
public:
// members
    Topography topo; // 'topography' nested within 'species'

    // parameters
    double c1; // interspecific competition parameter 1 - (disc: c_ij; cont: alpha)
    double c2; // interspecific competition parameter 2 - (disc: P(c_ij); cont: beta)
    double c3; // interspecific competition parameter 3 - connectance for continuously distributed A_ij
    double rho; // consumer mortality
    double sigma; // standard deviation log-normal attack rate distribution
    double alpha; // base attack rate
    double pProducer = 0.5; // probability of invading a producer species
    double emRate; // emigration rate
    double dispL; // dispersal length
    double thresh = 1e-4; // detection/extinction threshold
    double sigma_t = 0.05; // standard deviation of Ohrstein-Uhlenbeck process
    double omega = 0.5; // parameter controlling the width of the temperature niche
    double delta_g = 2.0; // (0.5x) proprtional difference between range in environmental and environmental optima
    double sigma_r = 0.25; // standard deviation of white noise added to quadratic environmental response function
    // switches
    bool prodComp = true; // select producer competition on/off
    int comp_dist = 0; // select competition distribution: 0 - discrete, 1 - pure beta, 2 - discretized beta
    bool symComp = false; // select symmetric competition
    int dispNorm = 0; // effort (0), degree normalized (1) or passive dispersal (2)

    // matrix objects
    mat bMat_p; // PxN biomass matrix - producers
    mat bavMat_p; // PxN biomass average matrix - producers
    mat bMat_p_src; // PxN biomass matrix - producers, source only
    mat rMat; // PxN r matrix - producers
    mat sMat; // PxN r matrix ignoring temperature dependence - producers
    mat cMat; // PxP competitive overlap matrix - producers
    mat uMat_p; // fixed unknowns - producers
    mat tMat; // environmental tolerances for explicit aboitic modelling
    mat dMat_n; // (sub)domain dispersal matrix
    mat dMat_m; // inter-subdomain dispersal matrix
    mat emMat_p; // species specific emigration rates - producers
    mat ouMat; // Ohrstein-Uhlenback process
    mat efMat; // environmental fluctuations centred on 0 for perturbing R
    vec bias; // species specific bias of environmental fluctuation
    mat gMat; // Sxl matrix of species environmental optima

    // storage objects
    mat trajectories; // matrix that will store the trajectories for analysis
    mat fluctuations; // matrix that will store time dependent growth rate matrices during relaxation
    umat sppRichness; // vector recording species richness as a function of time T

    // counters
    int invasion = 0; // counter used for recording number of invasions
    int invEvent = 0; // number of invasion events
    int inv0 = 0; // record number invasions at import
    int S_p = 0; // producer species richness
    int S_c = 0; // consumer species richness
    int I_p = 0; // counter for multispecies invasions
    int I_c = 0; // counter for multispecies invasions

// methods
    // species modelling
    void genDispMat(); // generate dispersal matrix
    mat genRVec(rowvec zVecExt = {}); // generate a spatially correlated random field of maximum internal growth rates, random variables can be passed from outside
    mat genRVecTemp(); // generate a spatially correlated random field mapped onto the output of a temperature response function
    mat genRVecQuad(); // generate growth rate vector based on quadratic environmental response function
    void updateRVecTemp(); // regenerate internal growth rate matrix after discrete warming event
    void ouProcess(); // updates efMat for modelling temporal abiotic turnover of temporal autocorrelation sigma_t, and bias mu
    mat genRVecERF(); // generate random specific environmental tolerance vector (outcome of implicit Environmental Response Function) and use to generate r_i
    void invade(int trophLev, bool test = true); // add new producer (trophLev=0) or cosumer (1) to specified 'port' patch
    field<uvec> extinct(int wholeDom=1, uvec ind_p = {}, uvec ind_c = {}); // remove extinct species, wholeDom flag to indicate if species testing required (whole domain only)
    void subSet(int domain); // subset whole domain objects into designated subdomain objects

    // (default) constructor
    Species () {}

    // (default) deconstructor
    ~Species () {}
};

#endif //SETC_SPECIES_H

// Local Variables:
// c-file-style: "stroustrup"
// End:
