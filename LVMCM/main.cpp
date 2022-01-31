////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// The parallelizable Lotka-Volterra Metacommunity assembly Model (LVMCM) ////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Jacob Dinner O'Sullivan -- j.osullivan@qmul.ac.uk | j.osullivan@zoho.com ///////////////////////
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
 * NEED TO UPDATE THE COMPUTION OF U*D_m IN PARALLEL ALGORITHM TO ACCOMODATE THE DISTRIBUTION IN EMIGRATION RATES
 */

#include <iostream>
#include <stdio.h>
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
#include <chrono>
#include <ctime>

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
bool ASSEMBLY_VIDEO = false;

// Storage for setjmp/longjmp
jmp_buf jump_buffer;
jmp_buf jump_buffer_continue_assembly;

// Global Metacommunity object and timing variables
Metacommunity meta;
std::chrono::time_point<std::chrono::system_clock> time1, time2;
std::chrono::duration<double> elapsed_time;
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

    cout << "WARNING: Handling of compartment indices after regional invasion AND extinction not currently compatible with bipartite model" << endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Store program arguments in variables ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Default Metacommunity parameters
    int a_init = 1;
    string a_bMat = {};
    string a_xMat = {};
    string a_scMat = {};
    string a_envMat = {};
    int a_invMax = 0;
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
    double a_c3 = -1;
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
    int a_lattice_width = 0;
    int a_lattice_height = 0;
    double a_phi = 5;
    int a_envVar = 0;
    vec a_skVec = {0.0};
    double a_var_e = 0.01;
    bool a_randGraph = true;
    bool a_gabriel = true;
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
    string a_path = {};

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
                    case 'n' : // set normalization of dispersal model
                        a_dispNorm = atoi(argv[i]);
                        break;
 
                }
                break;

            case 'e' :
                switch (var2) {
                    case '0' : // set envVar - number of explicitly modelled environmental variables (int)
                        a_envVar = atoi(argv[i]);
                        break;
                    case 'm' : // set path to data for environment matrix
                        a_envMat = argv[i];
                        break;
                }
                break;

            case 'f' : // set outputDirectory - location for write to file (string)
                a_outputDirectory = argv[i];
                break;

            case 'g' : // select maximum gamma diversity
                       // CAUTION: if regional limit is less than requested diversity, model will not converge!
                a_g_max = atoi(argv[i]);
                a_invasionSize = 0; // This will ensure gamma is not exceeded but if gamma set high, will slow assembly
                break;

            case 'h' : // set h - height of lattice
                a_randGraph = false;
                a_lattice_height = atoi(argv[i]);
                break;

            case 'i' :
                switch (var2) {
                    case '0' : // set invMax - total number of invasions (int)
                        a_invMax = atoi(argv[i]);
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
                        a_skVec.print("a_skVec_m");
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

            case 'w' : // set w - width of lattice
                a_lattice_width = atoi(argv[i]);
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

            case 'V' : // select save matrices for assembly video
                ASSEMBLY_VIDEO = true;
                a_path = argv[i];

            case 'X' : // set form of dynamics: 0 - LVMCM; 1 - PSD
                g_form_of_dynamics = atoi(argv[i]);
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
                switch(var2) {
                    case '0' :
                        ASSEMBLE = false;
                        a_autoFluct = atof(argv[i]);
                        break;
                    case 'T' : // select symmetric competition model
                        ASSEMBLE = false;
                        a_autoFluct = atof(argv[i]);
                        g_block_transitions = true;
                        break;
                }
                break;

            case 'W' : // set dTdt - select temperature warming experment and set rate of warming (double)
                a_dTdt = atof(argv[i]);
                WARMING = true;
                ASSEMBLE = false;
                break;
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// Set return clause /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        Metacommunity meta;

        if (setjmp(jump_buffer)) { // set return call in longjump clause
            time_t finish;
            time(&finish);
            cout << "\n\nSimulation stated at: " << ctime(&start);
            cout << "Simulation finished at: " << ctime(&finish);
            return (0);
        }

        time1 = std::chrono::system_clock::now();

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

        if (((a_lattice_height != 0) && (a_lattice_width == 0)) || ((a_lattice_height == 0) && (a_lattice_width != 0))) { // check lattice dimensions
            cout <<  a_lattice_height << " " << a_lattice_width << endl;
            cout << "\nError: Missing lattice dimension!" << endl;
            longjmp(jump_buffer, 1);
        }

        if (a_outputDirectory.length() == 0) { // check output directory set
            cout << a_outputDirectory << endl;
            cout << "Warning: no output directory given" << endl;
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// Set random seeds 1  //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (FIX_SEED > 0) {
            g_seed = 1;
        } else {
            // generate random seed
            std::random_device rd;
            g_seed = rd();
        }

        LVMCM_rng::boost_rng.seed(g_seed);
        arma_rng::set_seed(g_seed);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Parameterize assembly model ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        { // initialization scope
            meta = Metacommunity(
                    // Metacommunity parameters
                    a_init,
                    a_bMat,
                    a_xMat,
                    a_scMat,
                    a_invMax,
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
                    a_lattice_height,
                    a_lattice_width,
                    a_phi,
                    a_envVar,
                    a_skVec,
                    a_var_e,
                    a_randGraph,
                    a_gabriel,
                    a_T_int,
                    a_envMat,
                    // output variables
                    a_parOut,
                    a_experiment,
                    a_rep);
            meta.printParams();
        } // initialization scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Generate topo - topography/environment /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        { // abiotic modelling scope
            if (a_init) {
                meta.spp.topo.genLandscape(); // generate landscape
                meta.spp.genDispMat();
            } else {
                meta.spp.topo.genLandscape(
                    meta.spp.topo.network); // import network
            }
        } // abiotic modelling scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Simulation /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (FIX_SEED == 1) {
        // generate random seed for invaders
        std::random_device rd;
        g_seed = rd();

        LVMCM_rng::boost_rng.seed(g_seed);
        arma_rng::set_seed(g_seed);
    }

    setjmp(jump_buffer_continue_assembly);

    if (ASSEMBLE) {
        { // assembly scope
            // initialise objects
            int t;
            int no_invaders_p, no_residents_p, no_extinct_p, no_invaders_c, no_residents_c, no_extinct_c; // counters

            LVMCM_rng::boost_rng.seed(g_seed);
            arma_rng::set_seed(g_seed);

            do {

                meta.spp.invEvent++;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// Generate invader /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                { // invader testing scope
                    // select number of invaders and sample their trophic levels
                    int S_p0 = meta.spp.S_p; // record diversity prior to invader testing
                    int S_c0 = meta.spp.S_c;

                    // select number of invaders and sample their trophic levels
                    int no_trophLev[2] = {0}; // no producers, consumers to invade
                    int no_invaders = a_invasionSize * (meta.spp.S_p + meta.spp.S_c) + 1;

                    if (meta.spp.pProducer < 1) { // bipartite
                        // random uniform variables generated for allocating trophic level
                        vec trophLev(no_invaders);
                        int seedProd = 5; // minimum number of producer prior to introducing consumers
                        if (meta.spp.xMat.n_rows > seedProd) { // seed metacommunity with at least seedProd producers
                            trophLev.randu();
                        } else {
                            trophLev.zeros();
                        }
                        uvec invaderIndex = find(trophLev <= meta.spp.pProducer);
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
                            // sample random invaders and simulate dynamics
                            spp_tested += no_trophLev[tL]*2; // ?
                            meta.invaderSample(tL, no_trophLev[tL]);
                        }
                    }

                    meta.spp.invasion += no_invaders; // record cummulative number of invasions
                } // invader testing scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////// Simulate dynamics ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                { // dynamics scope
                    meta.metaCDynamics(meta.tMax); // simulate metacommunty dynamics
                } // dynamics scope

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Remove extinct species /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                { // extinction scope
                    meta.spp.extinct(); // remove extinct species
                    meta.spp.S_p = meta.spp.rMat.n_rows; // update richness counters
                    meta.spp.S_c = meta.spp.xMat.n_rows - meta.spp.rMat.n_rows; // update richness counters
                } // extinction scope

                if (BETA_T_INVASION > 0) {
                    // compute time series after every invasion to demonstrate emergence of autonomous turnover
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Generate unsaturated trajectory /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    meta.storeTraj = -1 * meta.spp.invasion;
                    meta.metaCDynamics(BETA_T_INVASION); // simulate metacommunty dynamics
                    meta.storeTraj = 0;
                }

                printf("\rInvasions / S_p / S_c = %d / %d / %d         ",
                       meta.spp.invasion, (int) meta.spp.rMat.n_rows, (int) meta.spp.xMat.n_rows - (int) meta.spp.rMat.n_rows);
                fflush(stdout);

                if (ASSEMBLY_VIDEO) { // saves key matrix objects after each invasion to show assembly process
                    meta.saveVideo(a_path);
                }
                
                if (a_g_max > 0) { // select regional number of species
                    if (meta.spp.xMat.n_rows >= a_g_max) { // stop simulation once externally defined regional diversity reached
                        sim = 0; // switch simulation off
                    }
                } else if (meta.spp.invasion >= meta.invMax) { // stop simulation once externally define no invasions reached
                    sim = 0; // switch simulation off
                }

                if (SNAPSHOT) { // write to file after every invasion
                    meta.saveMC();
                }

            } while (sim);
        } // assembly scope
    }

    time2 = std::chrono::system_clock::now();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Final book keeping /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (OUTPUT) { // write to file if requested

        if (ASSEMBLE) { // standard book keeping, assembly algorithm
            meta.metaCDynamics(1e4); // final relaxation
            meta.spp.extinct(); // final relaxation
            meta.spp.S_p = meta.spp.rMat.n_rows; // final relaxation
            meta.spp.S_c = meta.spp.xMat.n_rows - meta.spp.S_p; // final relaxation
            elapsed_time = time2 - time1;
            meta.simTime += elapsed_time.count(); // record assembly time
            meta.saveMC();

            printf("\n\nDetermining source-sink populations...");
            meta.genSourceSink();
            meta.saveMC();
        }

        if (SOURCE_SINK) { // if source sink only requested
            printf("\n\nGenerating source sink matrix...");
            meta.metaCDynamics(1e4); // final relaxation
            meta.spp.extinct();
            meta.spp.S_p = meta.spp.rMat.n_rows;
            meta.spp.S_c = meta.spp.xMat.n_rows - meta.spp.S_p;
            meta.saveMC();
            meta.genSourceSink();
            meta.saveMC();
        }

        if (CMAT_REG) { // if harvest only requested
            printf("\n\nGenerating regional scale interaction matrix...");
            meta.genCMatReg();
            meta.saveMC();
        }

        if (a_autoFluct) { // if time series only requested
            meta.metaCDynamics(1e4); // final relaxation
            meta.spp.extinct();
            meta.spp.S_p = meta.spp.rMat.n_rows;
            meta.spp.S_c = meta.spp.xMat.n_rows - meta.spp.S_p;
            meta.saveMC();
            printf("\n\nGenerating static environment trajectory...");
            meta.storeTraj = 1; // store trajectories in (NxS)xtRelax object
            meta.metaCDynamics(a_autoFluct);
        }

        if (FLUCTUATE) { // if time series, fluctuating environment requested
            printf(" Generating dynamic environment trajectory...");
            meta.storeTraj = 2; // concatenate static/dynamic environment trajectories
            int tRelax = 200; // relaxtion time for trajectories object
            for (int t = 0; t < tRelax; t++) { // simulate temporally fluctuating environment
                meta.envFluct();
            }
        }

        if (WARMING) { // if warming experiment requested
            // Begin with static environment step to demonstrate degree of autonomous fluctuations
            printf("\n\nGenerating static environment trajectory...");
            meta.storeTraj = 1; // store trajectories in (NxS)xtRelax object
            int tRelax = 100; // relaxtion time for trajectories object

            // Dial up intercept of temperature gradient, driving temperature optima to the right
            printf("\n\nSimulating regional warming...");
            meta.storeTraj = 1; // concatenate static/dynamic environment trajectories
            meta.spp.topo.T_int = 1.0; // required for importing models prior to adding T_int to output command
            tRelax = 10000; // relaxtion time for warming experiment object
            int res = 100;
            for (int t = 0; t < (tRelax / res); t++) { // simulate temporally fluctuating environment
                cout << "\nYear " << t << endl;
                meta.warming(a_dTdt, res, t);
            }
        }

        if (LONGDISTDISP) { // if long distance dispersal experiment requested
            printf("Randomly allocating long distance spatial coupling...\n");
            meta.longDistDisp(a_tMax, a_edges);
        }

        if (NODE_REMOVAL) { // if site removal experiment requested
            printf("Removing nodes and simulating dynamics...\n");
            meta.storeTraj=3;
            meta.nodeRemoval(a_tMax, a_nodeRemoval);
        }
    }
    longjmp(jump_buffer, 1); // jump to return
}

// Local Variables:
// c-file-style: "stroustrup"
// End:
