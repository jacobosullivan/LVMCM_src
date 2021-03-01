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

#include <iostream>
#include <armadillo>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <random>

#include "Species.h"
#include "Topography.h"
#include "LVMCM_rng.h"

using namespace std;
using namespace arma;
using namespace boost;

mat Species::genRVec(rowvec zVecExt) {

    // summary:
        // generate invader spatially autocorrelated growth rate vector
        // modelled using a Gaussian random field generated via spectral decomposition (https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution)
        // environment modelled implicitly 'through the eyes' of species
        // species growth rate distributions are independent

    // arguments:
        // zVecExt - vector of random variables passed from outside, required for modelling time dependent growth rates

    // required members:
        // topo.envVar - number of explicitly modelled environmental variables
        // topo.var_e - variance of the environmental distribution(s)
        // topo.phi - spatial autocorrelation length of the environmental distribution(s)
        // topo.distMat - euclidean distance matrix

    // output:
        // r_i - spatially autocorrelated growth rate vector

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Generate and decompose covariance matrix /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double mu; // mean of growth rate distribution
    if (zVecExt.n_cols == 0) {
        mu = 1.0;
    } else {
        mu = 0.0; // mean set to zero for OU process used to model time dependent growth rates
    }

    vec zVec;
    if (zVecExt.n_cols == 0) {
        zVec.randn(topo.no_nodes); // generate random normal variables
    } else {
        zVec = zVecExt.t();
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Sample from Gaussian random field /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vec r_i;
    r_i = mu + (topo.sigEVec * topo.sigEVal) * zVec; // sample from random field

    if (topo.network.n_rows != topo.no_nodes) { // subset by subdomain
        r_i = r_i.rows(topo.subdom_nodes);
    }

    return(r_i.t()); // return row vector of intrisic growth rates
}

mat Species::genRVecTemp() {

    // summary:
    // each species is assigned a random temperature optimum T_i and autocorrelated vector of intrinsic growth rates S
    // growth rate vector R then defined as S - (T(x) - T_i)^2
    // T(x) = T_int - sqrt(N)*x_1 where x_1 is first component of the cartesian coordinates of the node x

    // required members:
    // topo.T_int - intercept of the temperature distribution - can be 'dialed' to model temperature shifts

    // required methodsL
    // Species::genRVec()

    // output:
    // updates to tMat - matrix of species temperature optima
    // r_i - spatially autocorrelated growth rate vector

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Sample specific temperature optima ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    rowvec T_i;
    T_i.randu(1);
    T_i(0) *= topo.T_int; // randomly sample specific temperature optimum in range 0 <= T_i <= T_int
    rowvec T_iVec(topo.network.n_rows);
    T_iVec.fill(T_i(0));

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Generate spatially autocorrelated noise /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    rowvec s_i(topo.no_nodes);
    s_i = genRVec(); // generate spatially autocorrelated abiotic random field

    if (sMat.n_rows == 0) {
        sMat.set_size(1,topo.no_nodes);
        sMat = s_i;
    } else {
        sMat = join_vert(sMat, s_i);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// Generate growth rate vector ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    rowvec r_i(topo.no_nodes);

    if (omega > 0) { // omega fixed for all species

        double parab_width = pow(2/omega, 2); // calculate width parameter for parabolic temperature response function
        r_i = s_i - parab_width * pow(topo.envMat - T_iVec,2); // r_i generated using quadratic off-set from temperature optimum function
        uvec neg_growth = find(r_i < -1);
        r_i.elem(neg_growth).fill(-1); // set negative growth rates to zero for computation efficiency - result identical

        if (tMat.n_rows == 0) {
            tMat.set_size(1,1);
            tMat.row(0) = T_i;
        } else {
            tMat = join_vert(tMat, T_i);
        }
    } else {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Sample species temperature niche width /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // omega sampled from a discrete uniform distribution
        // values chosen set such that positive domain of the temperature niche covers a fraction p of the landscape

        bool discr_niche_width = false;
        double parab_width; // parabola width parameter
        double omega_i; // proportion of landscape covered by positive domain of parabola

        if (discr_niche_width) {

            // discrete integer random numbers
            typedef boost::random::mt19937 RandomIntGenerator;
            typedef boost::uniform_int<> UnifIntDistribution;
            typedef boost::variate_generator<RandomIntGenerator &, UnifIntDistribution> Generator;


            if (0) { // still testing the number of discrete width classes:
                // generate integers in range 1-3
                Generator getRandomInt(LVMCM_rng::boost_rng, UnifIntDistribution(1, 3));
                int disc_spat_niche = getRandomInt(); // niche width class for species i

                switch (disc_spat_niche) {
                    case 1:
                        parab_width = 100;
                        omega_i = 0.2;
                        break;
                    case 2:
                        parab_width = 25;
                        omega_i = 0.4;
                        break;
                    case 3:
                        parab_width = 11.1;
                        omega_i = 0.6;
                        break;
                }
            } else {
                // generate integers in range 1-5
                Generator getRandomInt(LVMCM_rng::boost_rng, UnifIntDistribution(1, 5));
                int disc_spat_niche = getRandomInt(); // niche width class for species i

                switch (disc_spat_niche) {
                    case 1:
                        parab_width = 400;
                        omega_i = 0.1;
                        break;
                    case 2:
                        parab_width = 100;
                        omega_i = 0.2;
                        break;
                    case 3:
                        parab_width = 44.4;
                        omega_i = 0.3;
                        break;
                    case 4:
                        parab_width = 25;
                        omega_i = 0.4;
                        break;
                    case 5:
                        parab_width = 16;
                        omega_i = 0.5;
                        break;
                }
            }
        } else { // sample omega_i uniformly and calculate parab_width

            vec p_i(1);
            p_i.randu();
            omega_i = p_i(0);
            parab_width = pow(2/omega_i, 2);

        }

        r_i = s_i - parab_width * pow(topo.envMat - T_iVec,2); // r_i generated using quadratic off-set from temperature optimum function
        uvec neg_growth = find(r_i < -1);
        r_i.elem(neg_growth).fill(-1); // set negative growth rates to zero for computation efficiency - result identical

        if (tMat.n_rows == 0) { // species temperature niche width stored in second column of tMat
            tMat.set_size(1,2); // include non-uniform temperature niche width
            tMat(0,0) = T_i(0);
            tMat(0,1) = omega_i;
        } else {
            rowvec T_niche_i(2); // include non-uniform temperature niche width
            T_niche_i(0) = T_i(0);
            T_niche_i(1) = omega_i;
            tMat = join_vert(tMat, T_niche_i);
        }
    }

    return (r_i);
}

mat Species::genRVecQuad() {

    // summary:
    // each species is assigned random environmental optima g_ij
    // growth rate vector R then defined as R_i = R0 - s_k(E_k - g_ik)^2
    // s_k is the kth element of the vector which contains the quadratic
    // shape parameters i.e. the sensitivity to the kth environmental variable
    // In order to centre the distribution on 1, I set R0=0, then subtract mean R_i
    // and add 1 to each new vector

    // required members:
    // topo.sVec - quadratic shape parameters, randomly sampled when environment generated
    // topo.envMat - Nxl matrix of autocorrelated environmental variables

    // required methods:

    // output:
    // updates to gMat - matrix of species temperature optima
    // r_i - spatially autocorrelated growth rate vector

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Sample specific environmental optima ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    rowvec g_ik;
    g_ik.randu(topo.envVar); // gen l random uniform variables
    g_ik *= (2*delta_g); // scale range of g_ik
    g_ik -= delta_g; // translate range of g_ik

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// Generate growth rate vector ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    rowvec r_i(topo.no_nodes);
    r_i.ones(); // R0 set to one for all species
    for (int k=0; k<topo.envVar; k++) {
        rowvec emg = topo.envMat.row(k) - g_ik(k);
        r_i -= topo.skVec(k)*(emg%emg);
    }

    r_i.elem(find(r_i < -1.0)).fill(-1.0); // set all terms < -1 to -1;
    // this speeds up simulation since very large negative coefficients are excluded from dynamics;
    // some weakly negative terms are retained to ensure fast decay outside of enironmental niche.

    if (tMat.n_rows == 0) { // species environmental optima
        tMat.set_size(1,topo.envVar);
        tMat.row(0) = g_ik;
    } else {
        tMat = join_vert(tMat, g_ik);
    }

    return (r_i);
}

void Species::updateRVecTemp() {
    // summary:
    // T_int (intercept of temperature gradient) externally updated to model regional climate warming
    // updates rMat to reflect abiotic change

    topo.genTempGrad();
    mat rMat_cc(rMat.n_rows, rMat.n_cols);
    mat envMat_SN(rMat.n_rows, rMat.n_cols), tMat_SN(rMat.n_rows, rMat.n_cols);
    for (int i=0; i<envMat_SN.n_rows; i++) {
        envMat_SN.row(i) = topo.envMat.row(0);
        tMat_SN.row(i).fill(tMat(i,0));
    }

    if (omega > 0) { // niche widths fixed for all species
        rMat_cc = sMat - omega * pow(envMat_SN - tMat_SN,2); // r_i generated using quadratic off-set from temperature optimum function
    } else {
        mat rMat_t = pow(envMat_SN - tMat_SN,2);
        for (int i=0; i<rMat_t.n_rows; i++) {
            rMat_t.row(i) = tMat(i,1) * rMat_t.row(i); // include non-uniform niche width
        }
        rMat_cc = sMat - rMat_t; // r_i generated using quadratic off-set from temperature optimum function
    }

    uvec neg_growth = find(rMat_cc < 0);
    rMat_cc.elem(neg_growth).fill(-1); // set negative growth rates to zero for computation efficiency - result identical
    rMat = rMat_cc; // over-write growth rate matrix

}

void Species::ouProcess() {

    // summary:
        // simulated abiotic turnover by updating state of matrix of environmental fluctuations
        // efMat updated using spatially autocorrelated ouMat (Ohrstein-Uhlenbeck random process)
        // OU process: Z(t) = (Z(t-1) + sigma_t.z) / sqrt(1 + sigma_t^2)
        // Output of OU random walk passed to genRVec
        // rMat(t) = rMat(0) + efMat(t)
        // if mu != 0 (moving average process), each species responds positively or negatively to abiotic turnover with probability 0.5

    // arguments:
        // sigma_t - temporal autocorrelation (standard deviation of white noise)
        // mu - absolute value of mean of random process, set to zero for stationary OU process

    // output:
        // efMat - current state of environmental fluctuation
        // ouMat process current stat of OU process

    // external function calls:
        // genRVec()

    if (ouMat.n_rows == 0) { // initialize matrices
        ouMat.zeros(bMat_p.n_rows, bMat_p.n_cols);
        efMat.zeros(bMat_p.n_rows, bMat_p.n_cols);
    }

    mat Zmat;
    Zmat.randn(efMat.n_rows, efMat.n_cols); // generate matrix of standard normal random variables
    double mu = 0;
    ouMat = (ouMat + mu + sigma_t * Zmat) / sqrt(1 + sigma_t*sigma_t); // update OU process

    for (int i=0; i<rMat.n_rows; i++) { // use output of OU process to seed random field simulator
        efMat.row(i) = genRVec(ouMat.row(i));
    }
}

mat Species::genRVecERF() {

    // summary:
        // each species is assigned a random environmental response vector which is used to generate r_i by multiplication with the environmental variables
        // tolerances are sample from a spherical distribution so that R is normally distributed

    // required members:
        // topo.envVar - number of explicitly modelled environmental variables
        // topo.envMat - matrix of environmental distribution vectors, dimensions (envVar x no_nodes)

    // output:
        // updates to tMat - matrix of species environmental response functions
        // r_i - spatially autocorrelated growth rate vector

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Sample specific environmental tolerance coefficients ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    rowvec tol;
    tol.randn(topo.envVar);
    rowvec u = pow(tol,2);
    tol = tol / sqrt(sum(u)) * sqrt(topo.envVar);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Generate local internal growth rates ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    rowvec r_i(topo.no_nodes);
    r_i = 1 + (tol * topo.envMat) / sqrt(2 * topo.envVar);

    if (tMat.n_rows == 0) {
        tMat.set_size(1,topo.envVar);
        tMat = tol;
    } else {
        tMat = join_vert(tMat, tol);
    }

    return (r_i);
}

void Species::invade(int trophLev, bool invTest) {

    // summary:
        // generate invader by sampling numerical traits and updating model matrices

    // arguments:
        // port - node to which invader is introduced
        // trophLev - select producer (0) or consumer (1) species
        // invTest - if preparing matrices for invader testing, coupling of invader->residents is suppressed

    // required members:
        // c1, c2 - interspecific competition parameters
        // rho - consumer mortality rate
        // sigma - standard deviation log-normal attack rate distribution
        // alpha - base attack rate
        // pProducer - probability of invading a producer species, set to zero for competitive model
        // emRate - emigration rate
        // dispL - dispersal length
        // prodComp - select coupled (T), uncoupled (F) producer dynamics
        // comp_dist - select discrete (0), pure beta (1), or discretized beta (2) sampling for A_ij

    // output:
        // updates to matrices bMat_p, rMat, cMat, bMat_c, aMat

    double inv = 1e-6; // invasion biomass

    if (trophLev == 0) { // generate producer

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// Add row to bMat_p ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (bMat_p.n_rows == 0) { // replace with insert rows need to set n cols at initialization
            bMat_p.set_size(1, topo.network.n_rows);
            bMat_p.fill(inv);
        } else {
            bMat_p.insert_rows(S_p+I_p,1);
            bMat_p.row(S_p+I_p).fill(inv);
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Add row to rMat ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (topo.envVar == 0) { // implicitly modelled environment
            if (rMat.n_rows == 0) {
                rMat = genRVec();
            } else {
                rMat = join_vert(rMat, genRVec());
            }
        } else if (topo.T_int != -1.0) {
            if (rMat.n_rows == 0) {
                rMat = genRVecTemp();
            } else {
                rMat = join_vert(rMat, genRVecTemp());
            }
        } else { // explicitly modelled environment
            if (topo.skVec(0) > 0.0) {
                if (rMat.n_rows == 0) {
                    rMat = genRVecQuad();
                } else {
                    rMat = join_vert(rMat, genRVecQuad());
                }
            } else {
                if (rMat.n_rows == 0) {
                    rMat = genRVecERF();
                } else {
                    rMat = join_vert(rMat, genRVecERF());
                }
            }
        }

        if (topo.consArea_bin.n_rows > 0) {
            if (topo.consArea_multiplicative) {
                rMat.row(rMat.n_rows - 1) %= topo.consArea_bin.t();
            } else {
                rMat.row(rMat.n_rows - 1) -= topo.consArea_bin.t();
            }
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Add row to emMat ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (emRate < 0) { // sample emigration rate for uniform distribution
            vec emRate_i(1);
            emRate_i.randu();
            if (emMat_p.n_rows == 0) {
                emMat_p = emRate_i;
            } else {
                emMat_p.insert_rows(S_p+I_p,1);
                emMat_p.row(S_p+I_p-1) = emRate_i;
            }
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// Add row/col to cMat //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        cMat.insert_rows(S_p+I_p,1); // row/col added even in case prodComp = false, however matrix operations not performed
        cMat.insert_cols(S_p+I_p,1);

        if (prodComp) { // generate non-zero interspecific competition terms

            if (comp_dist == 0) { // sample from discrete distribution

                boost::random::binomial_distribution<int> distribution(1,c1);
                for (int i=0; i<S_p+I_p; i++) {
                    cMat(S_p+I_p, i) = c2 * distribution(LVMCM_rng::boost_rng);
                    if (symComp) {
                        cMat(i, S_p+I_p) = cMat(S_p+I_p, i);
                    } else {
                        cMat(i, S_p+I_p) = c2 * distribution(LVMCM_rng::boost_rng);
                    }
                }

            } else if (comp_dist == 1) { // sample from continuous (beta) distribution

                typedef boost::random::mt19937 RandomNumberGenerator;
                typedef boost::random::beta_distribution<> BetaDistribution;
                typedef boost::variate_generator<RandomNumberGenerator&, BetaDistribution> Generator;
                BetaDistribution distribution(c1,c2);
                Generator getRandomNumber(LVMCM_rng::boost_rng,distribution);
                for (int i=0; i<S_p+I_p; i++) {
                    cMat(S_p+I_p, i) = distribution(LVMCM_rng::boost_rng);
                    if (symComp) {
                        cMat(i, S_p+I_p) = cMat(S_p+I_p, i);
                    } else {
                        cMat(i, S_p+I_p) = distribution(LVMCM_rng::boost_rng);
                    }
                }

            } else if (comp_dist == 2) { // sample from discretized beta distribution

                boost::random::binomial_distribution<int> binom_distribution(1,c3);
                typedef boost::random::mt19937 RandomNumberGenerator;
                typedef boost::random::beta_distribution<> BetaDistribution;
                typedef boost::variate_generator<RandomNumberGenerator&, BetaDistribution> Generator;
                BetaDistribution beta_distribution(c1,c2);
                Generator getRandomNumber(LVMCM_rng::boost_rng,beta_distribution);

                for (int i=0; i<S_p+I_p; i++) {
                    cMat(S_p+I_p, i) = beta_distribution(LVMCM_rng::boost_rng);
                    cMat(S_p+I_p, i) *= binom_distribution(LVMCM_rng::boost_rng);
                    if (symComp) {
                        cMat(i, S_p+I_p) = cMat(S_p+I_p, i);
                    } else {
                        cMat(i, S_p+I_p) = beta_distribution(LVMCM_rng::boost_rng);
                        cMat(i, S_p+I_p) *= binom_distribution(LVMCM_rng::boost_rng);
                    }
                }
            } else if (comp_dist == 3) { // sample from pure normal distribution
                typedef boost::random::mt19937 RandomNumberGenerator;
                typedef boost::random::normal_distribution<> NormDistribution;
                typedef boost::variate_generator<RandomNumberGenerator&, NormDistribution> Generator;
                NormDistribution normal_distribution(c1,sqrt(c2));
                Generator getRandomNumber(LVMCM_rng::boost_rng,normal_distribution);

                for (int i=0; i<S_p+I_p; i++) {
                    cMat(S_p + I_p, i) = normal_distribution(LVMCM_rng::boost_rng);
                    if (symComp) {
                        cMat(i, S_p + I_p) = cMat(S_p + I_p, i);
                    } else {
                        cMat(i, S_p + I_p) = normal_distribution(LVMCM_rng::boost_rng);
                    }
                }
            } else if (comp_dist == 4) { // sample from discretized normal distribution

                boost::random::binomial_distribution<int> binom_distribution(1,c3);
                typedef boost::random::mt19937 RandomNumberGenerator;
                typedef boost::random::normal_distribution<> NormDistribution;
                typedef boost::variate_generator<RandomNumberGenerator&, NormDistribution> Generator;
                NormDistribution normal_distribution(c1,sqrt(c2));
                Generator getRandomNumber(LVMCM_rng::boost_rng,normal_distribution);

                for (int i=0; i<S_p+I_p; i++) {
                    cMat(S_p + I_p, i) = normal_distribution(LVMCM_rng::boost_rng);
                    cMat(S_p+I_p, i) *= binom_distribution(LVMCM_rng::boost_rng);
                    if (symComp) {
                        cMat(i, S_p + I_p) = cMat(S_p + I_p, i);
                    } else {
                        cMat(i, S_p + I_p) = normal_distribution(LVMCM_rng::boost_rng);
                        cMat(i, S_p+I_p) *= binom_distribution(LVMCM_rng::boost_rng);
                    }
                }
            }
        }

        cMat(S_p+I_p, S_p+I_p) = 1.0;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// Add trophic interactions  ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (pProducer < 1.0) { // trophic interaction coefficients 'below' competition matrix
            if (S_c+I_c > 0) {
                rowvec a_pc = randn<rowvec>(S_c+I_c);
                a_pc *= sigma;
                a_pc = alpha * exp(a_pc);
                cMat.submat(S_p+I_p, S_p+I_p+1, S_p+I_p, cMat.n_cols-1) = a_pc;
                cMat.submat(S_p+I_p+1, S_p+I_p, cMat.n_cols-1, S_p+I_p) = a_pc.t();
            }
        }

        I_p++;

    } else if (trophLev == 1) { // generate consumer

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// Add row to bMat_p /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        bMat_p.insert_rows(bMat_p.n_rows,1);
        bMat_p.row(bMat_p.n_rows-1).fill(inv);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// Add row/col to cMat //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        cMat.insert_rows(cMat.n_rows,1);
        cMat.insert_cols(cMat.n_cols,1);

        rowvec a_pc = randn<rowvec>(S_p+I_p);
        a_pc *= sigma;
        a_pc = alpha * exp(a_pc);

        cMat.submat(cMat.n_rows-1, 0, cMat.n_rows-1, S_p+I_p-1) = a_pc;
        cMat.submat(0, cMat.n_cols-1, S_p+I_p-1, cMat.n_cols-1) = a_pc.t();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Add row to emMat ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (emRate < 0) { // sample emigration rate for uniform distribution
            vec emRate_i(1);
            emRate_i.randu();
            if (emMat_p.n_rows == 0) {
                emMat_p = emRate_i;
            } else {
                emMat_p.insert_rows(bMat_p.n_rows-1,1);
                emMat_p.row(emMat_p.n_rows-1) = emRate_i;
            }
        }
        I_c++;
    }
}

field<uvec> Species::extinct(int wholeDom, uvec ind_p, uvec ind_c) {

    // summary:
        // scan biomass matrices for regionally extinct species (root) and remove corresponding vectors from model matrices (all processes)

    //arguments:
        // wholeDom - indicates function called from root process and initiates scan of biomass matrices
        // ind_p - vector of indices of extinct producers - passed to non-root process for vector removal
        // ind_p - vector of indices of extinct consumers - passed to non-root process for vector removal

    // required members:
        // thresh - detection/extinction threshold

    // output:
        // (root) 2D field of producer/consumer extiction indices
        // updates to matrices bMat_p, rMat, cMat, bMat_c, aMat

    field<uvec> indReturn(2); // return object - 2D uvec of prod/cons extinction indices
    int countS=0;
    uvec ind_ext; // temporary extinction index
    int S_tot = bMat_p.n_rows; // total species richness

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// Check for extinct producers ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (wholeDom) { // check for extinct species at regional scale
        ind_p.set_size(S_p);
        for (int i=0; i<S_p; i++) {
            ind_ext = find(bMat_p.row(i) > thresh); // search for present populations
            if (ind_ext.n_rows == 0) {
                ind_p(countS) = i; // if no populations present, add index to ind_p
                countS++;
            }
        }
        ind_p.resize(countS);
    }
    indReturn(0) = ind_p; // store extinct producer index

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Shed rows/cols model objects ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int i = ind_p.n_rows - 1; i >= 0; i--) {
        bMat_p.shed_row(ind_p(i));
        if (uMat_p.n_rows > 0) {
            uMat_p.shed_row(ind_p(i));
        }
        rMat.shed_row(ind_p(i));
        if (sMat.n_rows > 0) {
            sMat.shed_row(ind_p(i));
        }
        cMat.shed_row(ind_p(i));
        cMat.shed_col(ind_p(i));
        if (tMat.n_rows != 0) {
            tMat.shed_row(ind_p(i));
        }
        if (emMat_p.n_rows != 0) {
            emMat_p.shed_row(ind_p(i));
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Check for extinct consumers ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (wholeDom) { // as above for consumer species
        ind_c.set_size(S_c);
        countS=0;
        for (int i=rMat.n_rows; i<bMat_p.n_rows; i++) {
            ind_ext = find(bMat_p.row(i) > thresh);
            if (ind_ext.n_rows == 0) {
                ind_c(countS) = i;
                countS++;
            }
        }
        ind_c.resize(countS);
    }
    indReturn(1) = ind_c;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Shed rows/cols model objects ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int i = ind_c.n_rows - 1; i >= 0; i--) {
        bMat_p.shed_row(ind_c(i));
        if (uMat_p.n_rows > 0) {
            uMat_p.shed_row(ind_c(i));
        }
        cMat.shed_row(ind_c(i));
        cMat.shed_col(ind_c(i));
        if (emMat_p.n_rows != 0) {
            emMat_p.shed_row(ind_p(i));
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Update temporal species richness vector and S_p/c ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (sppRichness.n_rows == 0) {
        sppRichness.set_size(1,2);
        sppRichness(sppRichness.n_rows-1,0) = rMat.n_rows;
        sppRichness(sppRichness.n_rows-1,1) = bMat_p.n_rows - rMat.n_rows;
    } else { // resize sppRichness vector(s) after each invasion
        sppRichness.resize(sppRichness.n_rows+1,2);
        sppRichness(sppRichness.n_rows-1,0) = rMat.n_rows;
        sppRichness(sppRichness.n_rows-1,1) = bMat_p.n_rows - rMat.n_rows;
    }

    return(indReturn);
}

void Species::genDispMat() {

    // summary:
        // generate the regional dispersal operator
        // dMat(x,y) = (e / sum_y(exp(-d_xy/l)) * exp(-d_xy/l)
        // dMat(x,x) = -e
        // if dMat already exists, regenerate for long distance dispersal perturbation
        // dispNorm: 0 - effort weighted; 1 - degree weighted; 2 - passive

    // required members:
        // dispL - dispersal length
        // emRate - emigration rate
        // topo.distMat
        // topo.adjMat

    // external function calls:
        // genAdjMat()

    // output:
        // dMat_n - regional dispersal operator

    if (dispNorm == 0) {
        cout  << "(effort normalized dispersal model)..." << endl;
    } else if (dispNorm == 1) {
        cout  << "(degree normalized dispersal model)..." << endl;
    } else if (dispNorm == 2) {
        cout  << "(passive dispersal model)..." << endl;
    }

    if (topo.no_nodes == 1) {
        dMat_n.set_size(1,1);
        dMat_n.zeros();
    } else {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Generate exponential dispersal kernal /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (topo.distMat.n_rows == 0) {
            topo.genDistMat();
        }
        dMat_n = exp(-topo.distMat / dispL);
        if (topo.adjMat.n_rows == 0) {
            topo.genAdjMat();
        }
        dMat_n %= topo.adjMat;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// Normalise immigration rate ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // In early versions D_xy were normalized by the edge weights such that a greater proportion of biomass
    // departed toward the nearer adjacent nodes. This weak bias was used to further localize dispersal and
    // reduce the incidence of numerical errors in the estimation of the Jacobian and the regional interact
    // ion matrices. However, since the Jacobian method for estimating C assumes a fixed point, this approa
    // ch is not applicable to steady state dynamics such as seen in the case of a heteroclinic network. In
    // order that the long distance dispersal connectivity perturbation make sense and that biomasses are b
    // alanced, this bias should be removed and immigration terms normalized by degree ONLY. MAKE SURE THIS
    // IS CLEAR IN ALL WRITTEN MATERIAL.

    mat kMat(topo.no_nodes, topo.no_nodes); // matrix encoding degree of each patch
    kMat.zeros();
    for (int i=0; i<topo.no_nodes; i++) { // normalise by weight
        // simple degree normalization
        uvec kVec = find(abs(topo.adjMat.col(i)) == 1); // absolute means that new edges (-1) also included
        for (int j=0; j<topo.no_nodes; j++) {
            if (abs(topo.adjMat(i,j)) == 1) {
                if (dispNorm == 0) {
                    // 'effort' normalized dispersal matrix
                    kMat(i,j) = emRate / sum(dMat_n.col(i));
                } else if (dispNorm == 1) {
                    // degree normalized dispersal matrix
                    kMat(i,j) = emRate / kVec.n_rows;
                } else if (dispNorm == 2) {
                    // 'passive' dispersal
                    kMat(i,j) = emRate;
                }
            }
        }
    }

    if (emRate < 0) { // for non-uniform emRate, set to -1.0, emRate included in dynamics
        kMat = abs(kMat);
    }
    dMat_n(find(topo.adjMat == -1)).fill(1); // required for element-wise multiplication
    dMat_n = kMat % dMat_n; // new edges, if any, assigned e/k (distance independent)
    bool minus_e = false; // include -e term from diagonal of D
    if (minus_e) {
        dMat_n.diag().fill(-1*abs(emRate)); // the dispersal operator includes the (negative) emigration terms on the diagonal
    }
}

void Species::subSet(int domain) {

    // summary:
        // subset model by domain decomposition using restriction operators and column indices
        // generate the matrix of fixed unknowns representingan approximation of the subdomain interface

    // arguments:
        // domain - index of subdomain, i.e. rank of parallel process, from which function called

    // required members:
        // topo.fVec - indicator vector, elements of which denote subdomain allocation of nodes

    // output:
        // updates to matrices bMat_p, rMat, bMat_c, uMat_p, uMat_c, dMat_n, dMat_m
        // define topo.intIF and topo.adjIF which record the interface nodes

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Infer and store intra- and inter- subdomain edges /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uvec index; // storage indices of specified subdomain
    index = find(topo.fVec == domain);

    std::vector<long long unsigned int> internal, adjacent; // storage objects for internal and adjacent edges - push_back functionality
    for (auto it=index.begin(); it!=index.end(); it++) {
        for (int j=0; j<topo.no_nodes; j++) {
            if (topo.adjMat(*it, j) == 1) {
                uvec in = find(index == j);
                if (in.is_empty()) {
                    internal.push_back(*it);
                    adjacent.push_back(j);
                }
            }
        }
    }

    uvec iIF(&internal[0], internal.size()); // convert to armadillo objects
    uvec aIF(&adjacent[0], internal.size());
    topo.intIF = unique(iIF); // dereplicate
    topo.adjIF = unique(aIF);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Extract row/cols corresponding to specified subdomain /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    topo.network = topo.network.rows(index);
    if (topo.envMat.n_rows != 0) {
        topo.envMat = topo.envMat.cols(index);
    }
    if (bMat_p.n_rows != 0) { // imported objects
        uMat_p = bMat_p.cols(topo.adjIF);
        bMat_p = bMat_p.cols(index);
        rMat = rMat.cols(index);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Generate restriction operators for decomposing dispersal operator /////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    umat P_n(index.n_rows, topo.no_nodes), P_m(topo.adjIF.n_rows, topo.no_nodes); // restriction operators
    P_n.zeros();
    P_m.zeros();

    for (int i=0; i<index.n_rows; i++) {
        P_n(i,index(i)) = 1;
    }

    for (int i=0; i<topo.adjIF.n_rows; i++) {
        P_m(i,topo.adjIF(i)) = 1;
    }

    dMat_m = P_m * dMat_n * P_n.t(); // inter-subdomain dispersal operator
    dMat_n = P_n * dMat_n * P_n.t(); // intra-subdomain dispersal operator

    if (uMat_p.n_rows > 0) {
        uMat_p = uMat_p * dMat_m;
    }

    if (topo.distMat.n_rows != 0) {
        topo.distMat = topo.distMat.submat(index,index);
    }

    topo.subdomain = domain;
    topo.subdom_nodes = index;
}
