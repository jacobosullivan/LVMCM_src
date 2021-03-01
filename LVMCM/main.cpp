////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// The parallelizable Lotka-Volterra Metacommunity assembly Model (pLVMCM) ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Jacob Dinner O'Sullivan -- j.osullivan@qmul.ac.uk | j.osullivan@zoho.com ///////////////////////
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

#include <iostream>
#include <stdio.h>
#include <mpi.h>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <csignal>
#include <csetjmp>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "Metacommunity.h"
#include "LVMCM_rng.h"

using namespace std;
using namespace arma;
using std::vector;

// Assembly flags
bool ASSEMBLE = true; // run assembly algorithm or jump directly to generation of analysis object
bool SOURCE_SINK = false; // generate source-sink matrix without assembling
bool SNAPSHOT = false; // store snapshot biomass and growth rate matrices
bool OUTPUT = true; // select write to file
bool WARMING = false; // select warming experiment
bool LONGDISTDISP = false; // select introduction of long distance dispersal
bool CONSAREA = false; // select conservation area experiment
bool TRAJECTORY = false; // select generate and write trajectory object to file
bool FLUCTUATE = false; // select generate and write to file abiotic fluctuation
bool CMAT_REG = false; // select generate and write to file regional competitive overlap matrix
int BETA_T_INVASION = 0; // select generate trajectory matrices after each invasion
int FIX_SEED = 0; // if 0 , all seeds random, if 1 landscape seed fixed, if 2 all seeds fixed
bool NODE_REMOVAL = false; // select node removal experiment
string DD_SOL_PATH = {}; // path to location for dd solution output

// Storage for setjmp/longjmp
jmp_buf jump_buffer;
jmp_buf jump_buffer_continue_assembly;

// Global Metacommunity object and timing variables
Metacommunity meta;
time_t time1, time2;
double mpi_time1, mpi_time2;
int g_seed;

// CVode tolerances
double TolA = 1e-8;
double TolR = 1e-7;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Start of assembly algoritm ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

    time_t start;
    time(&start);
    int sim = 1; // simulation switch

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Store program arguments in variables ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Default Metacommunity parameters
    int a_init = 1;
    string a_bMat = {};
    string a_xMat = {};
    string a_scMat = {};
    int a_invMax = 0;
    int a_iterS = 2;
    int a_deltaT = 100;
    int a_tMax = 500;
    int a_perfRep = 0;
    int a_perfS = 0;
    double a_autoFluct = 0; // if nonzero, gives trajectory recording time
    bool a_harvest = false;
    string a_outputDirectory;
    bool a_continue_assembly = false;
    bool a_write_continue_assembly = false;
    double a_invasionSize = 0.05;
    int a_g_max = 0;
    int a_nodeRemoval = 0;

    // Default Species parameters
    double a_c1 = 0.5;
    double a_c2 = 0.5;
    double a_c3 = 0.0;
    double a_emRate = 0.1;
    double a_dispL = 0.1;
    double a_pProducer = 1.0;
    bool a_prodComp = true;
    bool a_symComp = false;
    double a_alpha = 0;
    double a_sigma = 0;
    double a_sigma_t = 0.1;
    double a_rho = 0;
    int a_comp_dist = 0;
    double a_omega = 0.5;
    int a_dispNorm = 0;

    // Default Topograpy parameters
    int a_no_nodes = 4;
    double a_phi = 1;
    int a_envVar = 0;
    vec a_skVec = {0.0};
    double a_var_e = 0.01;
    bool a_randGraph = true;
    bool a_gabriel = true;
    int a_bisec = 0;
    double a_T_int = -1;
    double a_dTdt = 0.0;
    int a_edges = 1;
    int a_cons_percent = 60;
    double a_cons_phi = 1e-4;
    double a_cons_rep = 10;

    // Default output variables
    double a_parOut = 0;
    string a_experiment = "DEFAULT";
    int a_rep = 0;
    string a_jobID = "NA";

    for (int i = 1; i<argc; i++) { // loop through program arguments an allocate to parameters
        char var1 = argv[i][1];
        char var2 = argv[i][2];
        if (!isalpha(var2)) {
            var2 = '0';
        }
        i++;

        switch (var1) {
            // input parameters
            case 'a' : // set alpha - base attack rate (double)
                a_alpha = atof(argv[i]);
                break;

            case 'b' : // bFile: path to data for importing initialized model (string)
                a_bMat = argv[i];
                a_init = 0;
                break;

            case 'c' : // set c1, c2 - competition parameters (2x double)
                a_c1 = atof(argv[i]);
                i++;
                a_c2 = atof(argv[i]);
                break;

            case 'd' :
                switch (var2) {
                    case '0' : // set emRate, dispL - emigration rate and dispersal length (2x double)
                        a_emRate = atof(argv[i]);
                        i++;
                        a_dispL = atof(argv[i]);
                        break;
                    case 'd' : // set jobID - job ID used for automatic checkpointing (string)
                        DD_SOL_PATH = argv[i];
                        break;
                    case 'n' : // set normalization of dispersal model
                        a_dispNorm = atoi(argv[i]);
                        break;

                }
                break;

            case 'e' : // set envVar - number of explicitly modelled environmental variables (int)
                a_envVar = atoi(argv[i]);
                break;

            case 'f' : // set outputDirectory - location for write to file (string)
                a_outputDirectory = argv[i];
                break;

            case 'g' : // select maximum gamma diversity
                       // CAUTION: if regional limit is less than requested diversity, model will not converge!
                a_g_max = atoi(argv[i]);
                a_invasionSize = 0; // This will ensure gamma is not exceeded but if gamma set high, will slow assembly
                break;

            case 'i' :
                switch (var2) {
                    case '0' : // set invMax - total number of invasions (int)
                        a_invMax = atoi(argv[i]);
                        break;
                    case 'd' : // set jobID - job ID used for automatic checkpointing (string)
                        a_jobID = argv[i];
                        break;
                    case 's' : // set proportion of extra new invaders in each iteration (double)
                        a_invasionSize = atof(argv[i]);
                        break;
                }
                break;

            case 'n' : // set N - number of nodes (int)
                a_no_nodes = atoi(argv[i]);
                break;

            case 'o' : // set parOut, experiment, rep - key parameter, experiment name, replicate number for output filenames (double, string, int)
                a_parOut = atof(argv[i]);
                i++;
                a_experiment = argv[i];
                i++;
                a_rep = atoi(argv[i]);
                break;

            case 'p' :
                switch (var2) {
                    case '0' : // set phi - spatial autocorrelation length of the environment (double)
                        a_phi = atof(argv[i]);
                        break;
                    case 'p' : // set pProducer - probabilty of sampling a producer, bipartite models (double)
                        a_pProducer = atof(argv[i]);
                        break;
                }
                break;

            case 'r' : // set rho - consumer respiration rate (double)
                a_rho = atof(argv[i]);
                break;

            case 's' :
                switch(var2) {
                    case '0' : // set sigma - trophic link distribution parameter (double)
                        a_sigma = atof(argv[i]);
                        break;
                    case 'c' : // scFile: path to data for scaling of local interaction matrix
                        a_scMat = argv[i];
                        break;
                    case 'i' : // set iterS, deltaT - Schwartz iteration and time window (2x int)
                        a_iterS = atoi(argv[i]);
                        i++;
                        a_deltaT = atoi(argv[i]);
                        break;
                    case 'k' : // sk: environmental sensitivity shape parameter
                        a_envVar = atoi(argv[i]);
                        a_skVec.set_size(a_envVar);
                        i++;
                        for (int s=0; s<a_skVec.n_rows; s++) {
                            a_skVec(s) =  atof(argv[i]);
                            if (s < a_skVec.n_rows-1) {
                                i++;
                            }
                        }
                        break;
                    case 't' : // set sigma_t - standard deviation of random environmental fluctuation (double)
                        a_sigma_t = atof(argv[i]);
                        break;
                }
                break;

            case 't' : // set tMax - relaxation time (int)
                switch(var2) {
                    case '0' :
                        a_tMax = atoi(argv[i]);
                        break;
                    case 'n' : // set omega - temperature niche width in units 1/sqrt(N) (double)
                        a_omega = atof(argv[i]);
                        a_T_int = 1.0;
                        a_envVar = 1;
                        break;
                }
                break;

            case 'v' : // set var_e - variance of base environmental distribution (double)
                a_var_e = atof(argv[i]);
                break;

            case 'x' : // xFile: path to data for importing spatial network
                a_xMat = argv[i];
                break;

            // Switches
            case 'C' : // set prodComp - select producer coupling on/off (bool)
                if (!strcmp(argv[i],"F")) {
                    a_prodComp = false;
                }
                break;

            case 'D' : // set comp_dist - select distribution from which A_ij are sampled
            // 0 - discrete, 1 - pure beta, 2 - discretized beta. If 2, an additional argument is passed defining the connectance
                a_comp_dist = atoi(argv[i]);
                if ((a_comp_dist == 2) || (a_comp_dist == 4)) {
                    i++;
                    a_c3 = atof(argv[i]);
                }
                break;

            case 'G' : // set gabriel - select Gabriel/complete graph (bool)
                if (!strcmp(argv[i],"F")) {
                    a_gabriel = false;
                }
                break;

            case 'O' : // set OUTPUT - select write to file (bool)
                if (!strcmp(argv[i],"F")) {
                    OUTPUT = false;
                }
                break;

            case 'R' : // set randGraph - select random spatial network/lattice (bool)
                if (!strcmp(argv[i],"F")) {
                    a_randGraph = false;
                }
                break;

            case 'S' : // set SNAPSHOT - select regular write to file (bool)
                switch(var2) {
                    case '0' :
                        if (!strcmp(argv[i],"T")) {
                            SNAPSHOT = true;
                        }
                        break;
                    case 'C' : // select symmetric competition model
                        if (!strcmp(argv[i],"T")) {
                            a_symComp = true;
                        }
                        break;
                    case 'S' : // gen source sink matrix
                        if (!strcmp(argv[i],"T")) {
                            ASSEMBLE = false;
                            SOURCE_SINK = true;
                        }
                        break;
                }
                break;

            case 'Z' : // set FIX_SEED: 0 - random seed generated randomly; 1 - network/environment fixed; 2 - all sampling fixed
                FIX_SEED = atoi(argv[i]);
                break;

            // Perturbation experiments/analysis objects
            case 'B' : // compute beta t after each invasion
                BETA_T_INVASION = atoi(argv[i]);
                break;
            case 'F' : // set tMax, cons_percent, cons_phi, cons_rep - select conservation area/fragmentation experiment and set
                       // relaxation time, percent landscape conserved, correlation length of binary conservation area, and replicate number
                       // (int, int, double, int)
                a_tMax = atoi(argv[i]);
                i++;
                a_cons_percent = atoi(argv[i]);
                i++;
                a_cons_phi = atof(argv[i]);
                i++;
                if (atoi(argv[i]) > 0) {
                    a_continue_assembly = true;
                    a_invMax = atoi(argv[i]);
                }
                i++;
                a_cons_rep = atoi(argv[i]);
                CONSAREA = true;
                ASSEMBLE = false;
                break;

            case 'H' : // set harvest - estimate regional scale interaction coefficients using harvesting experiment (bool)
                if (!strcmp(argv[i],"T")) {
                    a_harvest = true;
                    ASSEMBLE = false;
                    CMAT_REG = true;
                }
                break;

            case 'K' : // set nodeRemoval
                a_nodeRemoval = atoi(argv[i]);
                ASSEMBLE = false;
                NODE_REMOVAL = true;
                break;

            case 'L' : // set tMax, edges - select long distance dispersal experment and set
                       // relaxation time and number of edges per perturbation event (2x int)
                a_tMax = atoi(argv[i]);
                i++;
                a_edges = atoi(argv[i]);
                LONGDISTDISP = true;
                ASSEMBLE = false;
                break;

            case 'T' : // set autoFluct - store trajectory for studying autonomous fluctuations (bool)
                ASSEMBLE = false;
                a_autoFluct = atof(argv[i]);
                break;

            case 'W' : // set dTdt - select temperature warming experment and set rate of warming (double)
                a_dTdt = atof(argv[i]);
                WARMING = true;
                ASSEMBLE = false;
                break;
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Set return clause parallel algorithm ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Metacommunity meta;
    int rank = 0; // MPI program removed from current version

    if (setjmp(jump_buffer)) { // set return call in longjump clause
        if (rank == 0) {
            time_t finish;
            time(&finish);
            cout << endl << "Simulation stated at: " << ctime(&start);
            cout << "Simulation finished at: " << ctime(&finish);
        }
        return (0);
    }

    if (rank == 0) { // start timer and check for parameter errors
        mpi_time1 = time(nullptr);

        if (!a_randGraph) { // check perfect square in case of lattice
            double srN = sqrt(a_no_nodes);
            if (floor(srN) - srN != 0) {
                cout << "\nError: perfect square N expected" << endl;
                longjmp(jump_buffer, 1);
            }
        }
        if (a_emRate < 0) { // check abs(emRate) = 1, required for non-uniform emRate;
            if (a_emRate != -1.0) {
                cout << "\nError: for non-uniform emigration rate model, emRate = -1.0 required" << endl;
                longjmp(jump_buffer, 1);
            }
        }
        if (a_pProducer < 1) { // check non-zero trophic parameters
            if ((a_alpha == 0) || (a_sigma == 0) || (a_rho == 0)) {
                cout << "\nError: trophic parameters set to zero!" << endl;
                longjmp(jump_buffer, 1);
            }
        }
        if ((a_prodComp == 1) && (a_comp_dist<3)) { // check non-zero competitive parameters
            if ((a_c1 == 0) || (a_c2 == 0)) {
                cout << "\nError: competitive parameters set to zero!" << endl;
                longjmp(jump_buffer, 1);
            }
        }
        double no_sub = a_no_nodes / pow(2, a_bisec);
        if (floor(no_sub) - no_sub != 0) { // check bisec/no_nodes correspond
            cout << "\nError: N / 2^bisec non-integer" << endl;
            longjmp(jump_buffer, 1);
        }

        if (a_deltaT != 0) {
            double tM_dT = a_tMax / a_deltaT;
            if (floor(tM_dT) - tM_dT != 0) { // check tMax/deltaT correspond
                cout << "\nError: tMax / deltaT non-integer" << endl;
                longjmp(jump_buffer, 1);
            }
        }
        if (OUTPUT) {
            if (a_outputDirectory.length() == 0) { // check output directory set
                string homedir = getenv("HOME");
                a_outputDirectory = homedir+"/LVMCM_src/";
                cout << "\nWarning: no output directory given" << endl;
                cout << "Saving data into " << a_outputDirectory << endl << endl;
            }
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// Set random seeds 1  //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (rank == 0) {
        if (FIX_SEED > 0) {
            g_seed = 1;
        } else {
            // generate random seed
            std::random_device rd;
            g_seed = rd();
        }
    }

    cout << "Random seed seen by process " << rank << " " << g_seed << endl;

    LVMCM_rng::boost_rng.seed(g_seed);
    arma_rng::set_seed(g_seed);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Parameterize parallel assembly model ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (rank == 0) { // parameterize root process using program arguments
        printf("\nInitializing LVMCM\n\n");
    }

    { // subdomain initialization scope
        // parameterize subdomains
        if (a_init) { // parameterize all processes from program arguments

            meta = Metacommunity(
                    // Metacommunity parameters
                    a_init,
                    a_bMat,
                    a_xMat,
                    a_scMat,
                    a_invMax,
                    a_iterS,
                    a_deltaT,
                    a_tMax,
                    a_outputDirectory,
                    // Species parameters
                    a_c1,
                    a_c2,
                    a_c3,
                    a_emRate,
                    a_dispL,
                    a_pProducer,
                    a_prodComp,
                    a_symComp,
                    a_alpha,
                    a_sigma,
                    a_sigma_t,
                    a_rho,
                    a_comp_dist,
                    a_omega,
                    a_dispNorm,
                    // Topograpy parameters
                    a_no_nodes,
                    a_phi,
                    a_envVar,
                    a_skVec,
                    a_var_e,
                    a_randGraph,
                    a_gabriel,
                    a_bisec,
                    a_T_int,
                    // output variables
                    a_parOut,
                    a_experiment,
                    a_rep,
                    a_jobID);
            if (rank == 0) {
                cout << "output folder set to " << meta.outputDirectory << endl;
                meta.printParams();
            }
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Broadcast imported metacommunity objects   //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (!a_init) { // broadcast imported metacommunity model objects and parameters to non-root processes
            int S_p, S_c, prodComp, envVar;
            if (rank == 0) {
                meta = Metacommunity(
                        // Metacommunity parameters
                        a_init,
                        a_bMat,
                        a_xMat,
                        a_scMat,
                        a_invMax,
                        a_iterS,
                        a_deltaT,
                        a_tMax,
                        a_outputDirectory,
                        // Species parameters
                        a_c1,
                        a_c2,
                        a_c3,
                        a_emRate,
                        a_dispL,
                        a_pProducer,
                        a_prodComp,
                        a_symComp,
                        a_alpha,
                        a_sigma,
                        a_sigma_t,
                        a_rho,
                        a_comp_dist,
                        a_omega,
                        a_dispNorm,
                        // Topograpy parameters
                        a_no_nodes,
                        a_phi,
                        a_envVar,
                        a_skVec,
                        a_var_e,
                        a_randGraph,
                        a_gabriel,
                        a_bisec,
                        a_T_int,
                        // output variables
                        a_parOut,
                        a_experiment,
                        a_rep,
                        a_jobID);
                if (rank == 0) {
                    cout << "output folder set to " << meta.outputDirectory << endl;
                    meta.printParams();
                }
                S_p = meta.sppPool.rMat.n_rows;
                S_c = meta.sppPool.bMat_p.n_rows - S_p;
                prodComp = meta.sppPool.prodComp;
                envVar = meta.sppPool.topo.envVar;
            }
        }
    } // subdomain initialization scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Generate topo - topography/environment (root) ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    { // domain decomposition scope
        if (rank == 0) {
            if (a_init) {
                meta.sppPool.topo.genDomainDecomp(); // generate and decompose landscape
                meta.sppPool.genDispMat();
            } else {
                meta.sppPool.topo.genDomainDecomp(
                    meta.sppPool.topo.network); // decompose imported topo
            }

            unique(meta.sppPool.topo.fVec.t()).print("\nf.unique");
            if (ASSEMBLE) {
                meta = meta; // store copy of complete domain at root process for outputting and extinction testing
            }
        }
        meta.rank = rank; // store rank for signal handler
    } // domain decomposition scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Simulation /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (rank == 0) {
            printf("\nStarting assembly\n");
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// Set random seeds 2  //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (FIX_SEED == 1) {
        if (rank == 0) {
            // generate random seed for invaders
            std::random_device rd;
            g_seed = rd();
        }
        cout << "Random seed (species only) seen by process " << rank << " " << g_seed << endl;
        LVMCM_rng::boost_rng.seed(g_seed);
        arma_rng::set_seed(g_seed);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Initialize MPI buffers //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    setjmp(jump_buffer_continue_assembly);

    if (ASSEMBLE) {

        { // assembly scope

            // initialise objects
            mat B_p, BStore_p, BStore_c, U_p, Ub_p, Tr, R, S, Cr, Cc; // MPI buffers
            mat B_dd; // for storage of domain decomposed solution if required
            int t, it, T = meta.tMax/meta.deltaT; // number of Schwarz time windows
            int no_invaders_p, no_residents_p, no_extinct_p, no_invaders_c, no_residents_c, no_extinct_c; // counters

            // reset seed for invader sampling (required to ensure all processes synced for invader sampling)
            LVMCM_rng::boost_rng.seed(g_seed);
            arma_rng::set_seed(g_seed);

            do {

                meta.sppPool.invEvent++;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Generate invader (all subdomains in parallel) ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                { // invader testing scope
                    // select number of invaders and sample their trophic levels
                    int S_p0 = meta.sppPool.S_p; // record diversity prior to invader testing
                    int S_c0 = meta.sppPool.S_c;

                    // select number of invaders and sample their trophic levels
                    int no_trophLev[2] = {0}; // no producers, consumers to invade
                    int no_invaders = a_invasionSize * (meta.sppPool.S_p + meta.sppPool.S_c) + 1;

                    if (meta.sppPool.pProducer < 1) { // bipartite
                        // random uniform variables generated for allocating trophic level
                        vec trophLev(no_invaders);
                        int seedProd = 5;
                        if (meta.sppPool.bMat_p.n_rows > seedProd) { // seed metacommunity with at least seedProd producers
                            trophLev.randu();
                        } else {
                            trophLev.zeros();
                        }
                        uvec invaderIndex = find(trophLev <= meta.sppPool.pProducer);
                        no_trophLev[0] = invaderIndex.n_rows; // number of producers to invade
                        no_trophLev[1] = no_invaders - no_trophLev[0]; // number of consumers to invade
                    } else {
                        no_trophLev[0] = no_invaders;
                    }

                    for (int tL = 0; tL < 2; tL++) { // first invade producers, then consumers
                        double spp_tested=0, suc_inv=0;

                        if (no_trophLev[tL] == 0) {
                            continue;
                        } else {

                            do {
                                // sample random invaders and simulate dynamics
                                spp_tested += no_trophLev[tL]*2;
                                uvec posGrowth = meta.invaderSample(tL, no_trophLev[tL]);

                                // select desired number of invaders with positive growth rates
                                suc_inv += posGrowth.n_rows;
                                posGrowth.resize(min((int) posGrowth.n_rows, no_trophLev[tL]));

                                // remove unsucessful/excess invaders
                                mat bInv_max;
                                bInv_max = meta.invaderCleanup(tL, posGrowth);

                                if (posGrowth.n_rows == 0) { // in this case all processes restart invader testing
                                    meta.invaderPopulate(tL, bInv_max);
                                    continue;
                                }

                                // populate interaction coefficients and reset invader biomass
                                meta.invaderPopulate(tL, bInv_max);
                                no_trophLev[tL] -= (int) posGrowth.n_rows;
                            } while (no_trophLev[tL] > 0);
                        }

                        // record numerical invasion probability
                        meta.invasionProb.resize(meta.sppPool.invEvent, 2);
                        if (tL == 0) {
                            meta.invasionProb(meta.sppPool.invEvent - 1, tL) =
                                    suc_inv / spp_tested;
                        } else if (tL == 1) {
                            meta.invasionProb(meta.sppPool.invEvent - 1, tL) =
                                    suc_inv / spp_tested;
                        }
                    }

                    if (rank == 0) { // update cMat for both producer and consumer invasions
                        meta.sppPool.cMat = meta.sppPool.cMat; // straight copy since C is not spatially decomposed
                    }

                    meta.sppPool.invasion += no_invaders;
                } // invader testing scope

                { // dynamics scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Undecomposed numerical solution /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    meta.metaCDynamics(meta.tMax); // simulate metacommunty dynamics

                } // dynamics scope

                { // extinction scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Remove extinct species /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    meta.sppPool.extinct(); // remove extinct species
                    meta.sppPool.S_p = meta.sppPool.rMat.n_rows;
                    meta.sppPool.S_c = meta.sppPool.bMat_p.n_rows - meta.sppPool.rMat.n_rows;
                } // extinction scope

                if (BETA_T_INVASION > 0) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Generate unsaturated trajectory /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    meta.storeTraj = -1 * meta.sppPool.invasion;
                    meta.metaCDynamics(BETA_T_INVASION); // simulate metacommunty dynamics
                    meta.storeTraj = 0;
                }

                if (rank == 0) {
                    mat presAbs;
                    presAbs.zeros(meta.sppPool.bMat_p.n_rows, meta.sppPool.bMat_p.n_cols);
                    double thresh = 1e-4;
                    presAbs.elem(find(meta.sppPool.bMat_p > thresh)).ones();
                    rowvec alpha = sum(presAbs,0);

                    printf("\rInvasions / S_p / S_c = %d / %d / %d         ",
                           meta.sppPool.invasion, (int) meta.sppPool.rMat.n_rows, (int) meta.sppPool.bMat_p.n_rows - (int) meta.sppPool.rMat.n_rows);
                    fflush(stdout);

                    if (a_g_max > 0) { // select regional number of species
                        if (meta.sppPool.bMat_p.n_rows == a_g_max) {
                            sim = 0; // switch simulation off
                        }
                    } else if (meta.sppPool.invasion >= meta.invMax) {
                        sim = 0; // switch simulation off
                    }
                }
            } while (sim);
        } // assembly scope

        if (rank == 0) {
            mpi_time2 = time(nullptr);
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Final book keeping /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (rank == 0) {

            if (OUTPUT) {

                if (meta.sppPool.topo.consArea_bin.n_rows > 0) {
                    ASSEMBLE = false; // write to file handled by perturbation clause in case of continued assemblies
                    a_write_continue_assembly = true;
                }
                if (ASSEMBLE) {
                    meta.metaCDynamics(1e4); // final relaxation
                    meta.sppPool.extinct();
                    meta.sppPool.S_p = meta.sppPool.rMat.n_rows;
                    meta.sppPool.S_c = meta.sppPool.bMat_p.n_rows - meta.sppPool.S_p;
                    time(&time2); // stop timing
                    meta.simTime += time2 - time1; // record assembly time
                    meta.outputData();

                    printf("\n\nDetermining source-sink populations...");
                    meta.genSourceSink();
                    meta.outputData();
                }

                if (SOURCE_SINK) {
                    printf("\n\nGenerating source sink matrix...");
                    meta.metaCDynamics(1e4); // final relaxation
                    meta.sppPool.extinct();
                    meta.sppPool.S_p = meta.sppPool.rMat.n_rows;
                    meta.sppPool.S_c = meta.sppPool.bMat_p.n_rows - meta.sppPool.S_p;
                    meta.outputData();
                    meta.genSourceSink();
                    meta.outputData();
                }

                if (CMAT_REG) {
                    printf("\n\nGenerating regional scale interaction matrix...");
                    meta.genCMatReg();
                    meta.outputData();
                }

                if (a_autoFluct) {
                    meta.metaCDynamics(1e4); // final relaxation
                    meta.sppPool.extinct();
                    meta.sppPool.S_p = meta.sppPool.rMat.n_rows;
                    meta.sppPool.S_c = meta.sppPool.bMat_p.n_rows - meta.sppPool.S_p;
                    meta.outputData();
                    printf("\n\nGenerating static environment trajectory...");
                    meta.storeTraj = 1; // store trajectories in (NxS)xtRelax object
                    meta.metaCDynamics(a_autoFluct);
                }
            }
        }
    longjmp(jump_buffer, 1); // jump to return
}

// Local Variables:
// c-file-style: "stroustrup"
// End:
