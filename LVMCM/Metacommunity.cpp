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
 * This class contains the members and methods required for simulating metacommunity dynamics,
 * importing and outputing data
 */

#include <iostream>
#include <fstream>
#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random.hpp>
#include <vector>
#include <iomanip>
#include <sstream>
#include "Metacommunity.h"
#include "Topography.h"
#include "Species.h"
#include "ODE.h"
#include "CommunityDynamics.h"
#include "LVMCM_rng.h"
#include <fenv.h>
#include <sys/ioctl.h> // ioctl, TIOCGWINSZ
#include <fcntl.h>     // open
#include <unistd.h>    // close

using namespace std;
using namespace arma;
using namespace boost::numeric::ublas;
using namespace boost::filesystem;

int g_form_of_dynamics = 0; // 0: LVMCM; 1: PSD
int g_block_transitions = 0;

void Metacommunity::metaCDynamics(int T) {

    // summary:
        // numerically solve metacommunity dynamics for trajectory simulation (serial or parallel) or invader testing
        // trajectory approximation uses Sundials CVODE solver (see communityDynamics.h/.cpp)

    // arguments:
        // dispersal - select with/without dispersal term (for invader testing)

    // required memebers:
        // spp - object of class Species, where model state is stored
        // T - relaxation time
        // storeTraj - selects whether trajectory object, dimensions (N*S)x(T), should be stored (large memory cost)

    // external function calls:
        // ODE.h::integrate_until()

    // output:
        // updates to matrices xMat, and bMat_c

//    int res = 10; // set time step resolution for storing trajectories; relaxT must be multiple of res

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Initialize community dynamics machinery ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CommunityDynamics dynamics;  // ODE dynamical object
    dynamics.xMat = &spp.xMat;

    dynamics.indices_DP = &spp.indices_DP;
    dynamics.indices_S = &spp.indices_S;
    dynamics.rMat = &spp.rMat;
    if (spp.cMat.n_rows != 0) {
        dynamics.cMat = &spp.cMat;
    }
    dynamics.rho = &spp.rho;
    dynamics.bodymass = &spp.bodymass;
    dynamics.bodymass_inv = &spp.bodymass_inv;
    dynamics.mu = &spp.mu;
    if (spp.topo.scVec.n_rows != 0) {
        dynamics.scVec = &spp.topo.scVec;
        dynamics.scVec_prime = &spp.topo.scVec_prime;
    }
    dynamics.S_p = spp.S_p;
    dynamics.S_c = spp.S_c;

    // full dynamics - dispersal switched on
    if (spp.emMat.n_rows != 0) {
        dynamics.emMat = &spp.emMat;
    }
    dynamics.dMat = &spp.dMat;

    string Bpath;
    string bMatDir;
    string newFolder;
    if (storeTraj != 0) {
        // write time point to file in SxN biomass matrix
        // file directory reflects the parOut and rep of the current assembly

        std::size_t pos1 = bMatFileName.find_last_of("/");
        bMatDir = bMatFileName.substr(0,pos1);
        pos1 = bMatFileName.find_last_of(")");
        std::size_t pos2 = bMatFileName.find_last_of(".");
        newFolder = bMatFileName.substr(pos1+1, pos2-pos1-1);
        if (!g_block_transitions) {
            Bpath = bMatDir + "/" + newFolder + "/trajectory";
        } else {
            Bpath = bMatDir + "/" + newFolder + "/trajectory_wo_transitions";
        }

        cout << "\nSaving matrices to " << Bpath << endl;

        if (!exists(Bpath)) { // make directory if doesn't currently exist
            create_directories(Bpath);
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Metacommunity relaxation step /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // prints state matrix after each resPrint unit times
    bool PRINT_X_MAT = false;
//    bool PRINT_X_MAT = spp.invasion > 4*spp.rMat.n_rows;
    bool movie_mode = false;
    int resPrint=1; // temporal resolution printing to console
    int resSave=100; // temporal resolution writing to file (must be multiple of resPrint)
    int terminal_lines,terminal_columns;
    int species_per_row,rows_of_species;
    std::stringstream formatP,formatS,formatD;
    int digits_PD=2;
    int width_of_item=7+digits_PD;
    formatP << "\x1B[34m%9.0" << digits_PD << "e\033[0m";
    formatS << "\x1B[31m%9.0" << digits_PD+4 << "f\033[0m";
    formatD << "\x1B[32m%9.0" << digits_PD << "e\033[0m";

    if(PRINT_X_MAT && movie_mode){
	cout << "\033[2J";  // clear screen

	if(!spp.topo.randGraph){
	    // Get current size of terminal:
	    int fd;
	    struct winsize w;
	    fd = open("/dev/tty", O_RDWR);
	    if(fd < 0 || ioctl(fd, TIOCGWINSZ, &w) < 0)
		FATAL_ERROR("Can't get terminal size.");
	    close(fd);
	    terminal_lines = w.ws_row;
	    terminal_columns = w.ws_col;

	    species_per_row =
		    (terminal_columns+1) / 
		    ((width_of_item+1) * spp.topo.lattice_width);
	    rows_of_species =
		    min<int>((terminal_lines - 2) / (spp.topo.lattice_height+1) ,
			     (spp.xMat.n_rows-1) / species_per_row + 1 );
	}
    }
    
    while(dynamics.current_time < T) {

        ODE_state::root_indices_t root_indices;
        { // ODE_state state scope
            ODE_state state(&dynamics);
            for (double tNext = ceil(dynamics.current_time/resPrint+1e-7)*resPrint;
		 dynamics.current_time < T;
		 tNext = ceil(dynamics.current_time/resPrint+1e-7)*resPrint ) {

                state.integrate_to_root(tNext, root_indices);

			// Print only after reaching steady state (it's more fun)
                if (PRINT_X_MAT) {
                    if (dynamics.current_time >= tNext) {
			if (movie_mode) {// movie mode
			    usleep(1000*100); // sleep so many micro seconds
			    cout << "\033[H" ; //  move cursor to top
			    printf("Invasions / S_p / S_c = %d / %d / %d         \n",
				   spp.invasion, (int) spp.rMat.n_rows, (int) spp.xMat.n_rows - (int) spp.rMat.n_rows);
			}
                        mat X_tmp(&state[0], spp.xMat.n_rows, spp.xMat.n_cols, false, true);
                        string status_output = "\nt = " + to_string((int) tNext) + " X = ";
                        cout << status_output;
			if(!(movie_mode && !spp.topo.randGraph)){
			    // one-species-per-line printing
			    for (int i = 0; i < X_tmp.n_rows; i++) {
				cout << endl;
				for (int j = 0; j < X_tmp.n_cols; j++) {
				    int oneDind = i + X_tmp.n_rows * j;

				    if (((uvec) find(spp.indices_DP == oneDind)).n_rows > 0) {
					if (X_tmp(i, j) > spp.bodymass) {
					    printf(formatD.str().c_str(), X_tmp(i, j));
					} else {
					    printf(formatP.str().c_str(), X_tmp(i, j));
					}
				    } else if (((uvec) find(spp.indices_S == oneDind)).n_rows > 0) {
					printf(formatS.str().c_str(), X_tmp(i, j));
				    } else {
					cout << "ERROR: INDEX NOT FOUND!" << endl;
				    }
				    cout << " ";
				}
			    }
			}else{// print each species with it's own "map"
			    for (int R = 0; R < rows_of_species ; R++){ // row of species
				if (R > 0){
				    cout << endl;
				    int row_length=species_per_row * (width_of_item+1) * spp.topo.lattice_width - 1;
				    for (int c=0; c<row_length; c++){
					cout << ((c+1) % ((width_of_item+1) * spp.topo.lattice_width)==0 ? '+' : '-');
				    }
				}
				for (int r = 0; r < spp.topo.lattice_height; r++){ // row of patches
				    cout << endl;
				    for (int C = 0; C < species_per_row; C++){ // col of species
					int i = R*species_per_row+C; //species index
					if(i >= X_tmp.n_rows) continue;
					for (int c=0; c < spp.topo.lattice_width; c++){
					    cout << ( c==0 ? (C==0 ? "" : "|") : " ");						 
					    int j = r*spp.topo.lattice_width+c; // patch index
					    int oneDind = i + X_tmp.n_rows * j;
					    if (((uvec) find(spp.indices_DP == oneDind)).n_rows > 0) {
						if (X_tmp(i, j) > spp.bodymass) {
						    printf(formatD.str().c_str(), X_tmp(i, j));
						} else {
						    printf(formatP.str().c_str(), X_tmp(i, j));
						}
					    } else if (((uvec) find(spp.indices_S == oneDind)).n_rows > 0) {
						printf(formatS.str().c_str(), X_tmp(i, j));
					    } else {
						cout << "ERROR: INDEX NOT FOUND!" << endl;
					    }
					}
				    }
				}
			    }
			}
                    }
                }

                if (dynamics.current_time >= tNext && storeTraj != 0) {
                    if ((int) tNext % resSave == 0) {
                        // write current state to file
                        mat Btofile(&state[0], 1, dynamics.number_of_variables());
                        Btofile.reshape(spp.xMat.n_rows, spp.xMat.n_cols);
                        string Bfile;
                        Bfile = Bpath + "/bMat" + to_string((int) tNext) + ".mat";
                        Btofile.save(Bfile, raw_ascii);
                    }

                }
                if(!root_indices.empty()) {
                    break; // root was found
                }
            }
        } // ODE_state state scope

        dynamics.react_to_roots(root_indices);  // <- outside of ODE_state scope
        root_indices.clear();
    }
}

void Metacommunity::invaderSample(int trophLev, int no_invaders) {
    // summary:
        // introduce new random species and test for positive growth rates in single evaluation of ODE

    // arguments:
        // trophLev - 0: producer; 1: consumer
        // no_invaders - number of successful invaders to introduce

    // required memebers:
        // spp - object of class Species, where model state is stored

    // external function calls:
        // Species::invade()

    // output:
        // vector of positively growing species indices
        // for MPI program, this vector is distributed to all parallel processes and species are selected according to
        // deterministic algorithm

    // testing parameters
    double min_b = 0; // minimum biomass for inclusion in model
    double inv = 1e-6; // invasion biomass
    int S_before = spp.xMat.n_rows; // record previous richness for updating indicator vectors
    uvec invaderIndex, posGrowth;

    if (trophLev == 0) { // sample producer species
        // size of testing pool
        int invExcess_p = 3; // invade excess species to account for difficulty in finding successful invader -- set arbitrarily
        int suc_inv = 0; // count successful invasions
        do {

            for (int i = 0; i < invExcess_p * no_invaders; i++) {
                spp.invade(0); // invade a producer, spp.I_p += 1
            }

            // single evaluation of ODE to test for positive growth, invaders ONLY
            mat bInv;
            mat B;
            B.zeros(spp.xMat.n_rows, spp.xMat.n_cols);
            B.elem(find(spp.xMat > 0)) = spp.xMat.elem(find(spp.xMat > 0)); // fill B using indexing

            // Compute invader growth rate, test for positive growth
            bInv = mat (spp.rMat.rows(spp.S_p, spp.S_p + spp.I_p - 1) -
                        spp.cMat.submat(spp.S_p, 0,
                                            spp.S_p + spp.I_p - 1, spp.cMat.n_cols - 1) * B);

            // locate (and subset if needed) indices of species with positive invasion fitness
            vec bInv_max(bInv.n_rows);
            for (int i = 0; i < bInv.n_rows; i++) { // store maximum local biomass of each invader
                bInv_max(i) = bInv.row(i).max();
            }
            posGrowth = find(bInv_max >= min_b); // index vector of positively growing species
            posGrowth.resize(min((int) posGrowth.n_rows, no_invaders - suc_inv));
            suc_inv += posGrowth.n_rows; // dial up successful invasion counter

            // remove negative invasion fitness and excess species
            umat negGrowth(spp.I_p,1);
            negGrowth.col(0) = linspace<uvec>(spp.S_p, spp.S_p + spp.I_p - 1,
                                              spp.I_p);

            for (int i = posGrowth.n_rows - 1; i>=0; i--) { // negGrowth includes indices of species not in posGrowth
                negGrowth.shed_row(posGrowth(i));
            }

            for (int i = negGrowth.n_rows - 1; i >= 0; i--) {
                spp.xMat.shed_row(negGrowth(i)); // remove prod biomass vec
                if (spp.rMat.n_rows > 0) {
                    spp.rMat.shed_row(negGrowth(i)); // remove prod growth vec
                }
                if (spp.sMat.n_rows > 0) {
                    spp.sMat.shed_row(negGrowth(i)); // remove prod growth vec
                }
                spp.cMat.shed_row(negGrowth(i)); // remove prod comp term
                spp.cMat.shed_col(negGrowth(i));
                if (spp.topo.envVar != 0) {
                    spp.tMat.shed_row(negGrowth(i));  // remove prod env tol vec
                }
                if (spp.emRate < 0) {
                    spp.emMat.shed_row(negGrowth(i));  // remove prod emigration rate
                }
            }

            if (posGrowth.n_rows>0) {
                // populate xMat
                bInv = bInv.rows(posGrowth); // remove negative invasion fitness and excess from testing matrix

                for (int i=0; i<posGrowth.n_rows; i++) {
                    uvec iMax = find(bInv.row(i) == bInv.row(i).max()); // select highest invasion fitness
                    spp.xMat.row(spp.S_p + i).randu(); // all invader populations with positive invasion fitness (\in P) allocated Poisson Clocks
                    spp.xMat.row(spp.S_p + i) = log(spp.xMat.row(spp.S_p + i)) - 1; // -1 gives a buffer against very small negative biomasses in P compartment
                    uvec ind_probabilistic_inv = find(bInv.row(i) <= 0);
                    if (ind_probabilistic_inv.n_rows > 0) {
                        ind_probabilistic_inv = spp.S_p + i + spp.xMat.n_rows * ind_probabilistic_inv;
                        spp.xMat.elem(ind_probabilistic_inv).zeros(); // all invader populations with negative invasion fitness (\in P) set to 0
                    }
                    spp.xMat(spp.S_p + i, iMax(0)) = spp.bodymass; // add biomass
                }
            }

            spp.I_p = 0; // reset invader counter
            spp.S_p = spp.rMat.n_rows;

        } while (suc_inv < no_invaders);

    } else if (trophLev == 1) { // sample consumer species

        // size of testing pool
        int invExcess_c = 3; // invade excess species to account for difficulty in finding successful invader
        int suc_inv = 0; // count successful invasions
        do {

            for (int i = 0; i < invExcess_c * no_invaders; i++) {
                spp.invade(1); // invade a consumer
            }

            // single evaluation of ODE to test for positive growth, invaders ONLY
            mat B;
            B.zeros(spp.xMat.n_rows, spp.xMat.n_cols);
            B.elem(find(spp.xMat > 0)) = spp.xMat.elem(find(spp.xMat > 0)); // fill B using indexing

            // Compute invader growth rate, test for positive growth
            mat bInv = spp.rho *
                       (spp.cMat.submat(spp.S_p + spp.I_p + spp.S_c, 0, spp.cMat.n_rows - 1,
                                            spp.S_p + spp.I_p - 1) *
                        B.rows(0, spp.S_p + spp.I_p - 1) - 1);

            // locate (an subset if needed) indices of species with positive invasion fitness
            vec bInv_max(bInv.n_rows);
            for (int i = 0; i < bInv.n_rows; i++) { // store maximum local biomass of each invader
                bInv_max(i) = bInv.row(i).max();
            }
            posGrowth = find(bInv_max >= min_b);
            posGrowth = shuffle(posGrowth);
            posGrowth.resize(min((int) posGrowth.n_rows, no_invaders - suc_inv));
            suc_inv += posGrowth.n_rows; // dial up successful invasion counter

            // remove negative invasion fitness and excess species
            umat negGrowth(spp.I_c, 1);
            negGrowth.col(0) = linspace<uvec>(spp.S_p + spp.I_p + spp.S_c, spp.xMat.n_rows - 1,
                                              spp.I_c);

            for (int i = posGrowth.n_rows - 1; i >= 0; i--) { // negGrowth includes indices of species not in posGrowth
                negGrowth.shed_row(posGrowth(i));
            }

            for (int i = negGrowth.n_rows - 1; i >= 0; i--) {
                spp.xMat.shed_row(negGrowth(i)); // remove cons biomass vec
                spp.cMat.shed_row(negGrowth(i)); // remove interaction coefficients
                spp.cMat.shed_col(negGrowth(i));
                if (spp.emRate < 0) {
                    spp.emMat.shed_row(negGrowth(i));  // remove prod emigration rate
                }
            }

            if (posGrowth.n_rows > 0) {
                // populate xMat
                bInv = bInv.rows(posGrowth); // remove negative invasion fitness and excess from testing matrix

                for (int i=0; i<posGrowth.n_rows; i++) {
                    uvec iMax = find(bInv.row(i) == bInv.row(i).max()); // select highest invasion fitness
                    spp.xMat.row(spp.S_c + i).randu(); // all invader populations with positive invasion fitness (\in P) allocated Poisson Clocks
                    spp.xMat.row(spp.S_c + i) = log(spp.xMat.row(spp.S_c + i));
                    uvec ind_probabilistic_inv = find(bInv.row(i) <= 0);
                    if (ind_probabilistic_inv.n_rows > 0) {
                        ind_probabilistic_inv = spp.S_c + i + spp.xMat.n_rows * ind_probabilistic_inv;
                        spp.xMat.elem(ind_probabilistic_inv).zeros(); // all invader populations with negative invasion fitness (\in P) set to 0
                    }
                    spp.xMat(spp.S_c + i, iMax(0)) = spp.bodymass; // add biomass
                }
            }

            spp.I_c = 0; // reset invader counter
            spp.S_c = spp.xMat.n_rows - spp.rMat.n_rows;
        } while (suc_inv < no_invaders);
    }

    // Update indicator vectors
    // THIS ASSUMES A COMPETITION MODEL
    // UPDATE FOR BIPARTITE
    // NOTE THAT THE FORMULAE FOR TRANSLATING INDICES ARE WRONG IN THE CASE THAT INVADERS ARE NOT ADDED TO BOTTOM OF MATRIX

    if (g_form_of_dynamics) { // DEF LVMCM|PSD
        if (spp.xMat.n_rows == 1) {
            spp.indices_DP = find(spp.xMat > -0.5); //splitting at -0.5 is ...
            spp.indices_S = find(spp.xMat <= -0.5); //safter than splitting at -1
        } else {
            uvec site_index_DP = floor(spp.indices_DP / S_before);
            uvec site_index_S = floor(spp.indices_S / S_before);

            spp.indices_DP += (spp.xMat.n_rows - S_before) * site_index_DP;
            spp.indices_S += (spp.xMat.n_rows - S_before) * site_index_S;

            uvec invader_indices((spp.xMat.n_rows - S_before) * spp.xMat.n_cols);
            for (int j = 0; j < spp.xMat.n_cols; j++) {
                invader_indices.rows(j * (spp.xMat.n_rows - S_before),
                                     j * (spp.xMat.n_rows - S_before) + (spp.xMat.n_rows - S_before) - 1) =
                        linspace<uvec>((j * spp.xMat.n_rows) + (spp.xMat.n_rows - 1) -
                                       (spp.xMat.n_rows - S_before) + 1,
                                       (j * spp.xMat.n_rows) + (spp.xMat.n_rows - 1),
                                       (spp.xMat.n_rows - S_before));
            }

            uvec invader_DP = find(spp.xMat.rows(S_before, spp.xMat.n_rows - 1) >= 0);
            spp.indices_DP = sort(join_cols(spp.indices_DP, invader_indices.elem(invader_DP)));
            invader_indices.shed_rows(invader_DP);
            spp.indices_S = sort(join_cols(spp.indices_S, invader_indices));
        }
    } else {
        spp.indices_DP = linspace<uvec>(0,spp.xMat.n_elem-1,spp.xMat.n_elem);
    }
}

void Metacommunity::envFluct() {

    // summary:
        // update aboitic turnover and simulates a single CVode timestep

    // required members:
        // spp - object of class Species, where model state is stored

    // external function calls:
        // metaCDynamics()
        // spp.ouProcess()

    // output:
        // updates to spp matrices

    CommunityDynamics dynamics;  // ODE dynamical object
    dynamics.xMat = &spp.xMat;
    dynamics.rMat = &spp.rMat;

    if (spp.cMat.n_rows != 0) {
        dynamics.cMat = &spp.cMat;
    }

    dynamics.dMat = &spp.dMat; // if dynamics.dMat is initialized, single domain relaxation selected
    dynamics.rho = &spp.rho;

    spp.ouProcess(); // update the spatially resolved OU process to be added to the growth rate matrix
    dynamics.efMat = &spp.efMat; // in CommunityDynamics rMat += efMat
    ODE_state state(& dynamics);
    state.integrate_until(1); // simulate single timestep
    mat tmp(&state[0], 1, dynamics.number_of_variables());
    spp.trajectories = join_vert(spp.trajectories, tmp); // store trajectory
    mat RE = spp.rMat + spp.efMat;
    RE.reshape(1,RE.size());
    spp.fluctuations = join_vert(spp.fluctuations, RE); // store current distrition in R (= rMat + efMat)
}

void Metacommunity::warming(double dTdt, int res, int time) {

    // summary:
        // updates temperature gradient and rMat and simulates a single CVode timestep

    // required members:
        // spp - object of class Species, where model state is stored

    // external function calls:
        // metaCDynamics()
        // spp.updateRVecTemp()

    // output:
        // updates to spp matrices

    // write time point to file in SxN biomass matrix
    // file directory reflects the parOut and rep of the current assembly

    std::size_t pos1 = bMatFileName.find_last_of("/");
    string bMatDir = bMatFileName.substr(0,pos1);
    pos1 = bMatFileName.find_last_of(")");
    std::size_t pos2 = bMatFileName.find_last_of(".");
    string newFolder = bMatFileName.substr(pos1+1, pos2-pos1-1);
    string Bpath = bMatDir + "/" + newFolder + "/dTdt=" + to_string(dTdt);

    if (!exists(Bpath)) { // make directory if doesn't currently exist
        create_directories(Bpath);
    }

    CommunityDynamics dynamics;  // ODE dynamical object

    dynamics.xMat = &spp.xMat;

    if (spp.cMat.n_rows != 0) {
        dynamics.cMat = &spp.cMat;
    }
    dynamics.dMat = &spp.dMat; // if dynamics.dMat is initialized, single domain relaxation selected
    dynamics.rho = &spp.rho;

    for (int t=0; t<res; t++) { // update every time step, save every res timesteps
        spp.topo.T_int += dTdt; // apply warming - T_int updated
        spp.updateRVecTemp(); // apply warming - rMat updated
        printf("\rt = %d, T_int = %f", t, spp.topo.T_int);
        fflush(stdout);
        dynamics.rMat = &spp.rMat; // update rMat seen by CommmunityDynamics
        ODE_state state(& dynamics);
        state.integrate_until(1); // simulate single timestep
    }

    string Bfile = Bpath + "/bMat_w" + to_string(time) + ".mat";
    string Rfile = Bpath + "/rMat_w" + to_string(time) + ".mat";
    spp.xMat.save(Bfile, raw_ascii);
    spp.rMat.save(Rfile, raw_ascii);

    if (time == 0) {
        // store parameter file in same directory for convenience
        string filenameP = Bpath + "/pars.mat";
        writePars(filenameP);
    }
}

void Metacommunity::longDistDisp(int tMax, int edges) {

    // summary: ...

    bool stepwise = true; // stepwise addition of new edges

    if (edges < 0) {
        stepwise = false; // select single addition followed by high res storage of trajectories
    }

    // record row minimum of dispersal operator, if non-zero nodes have degree N-1
    colvec row_min = min(abs(spp.dMat),1);

    int non_zero_d;
    { // scoped for deletion of sparse cast D
        sp_mat D_sp(spp.dMat);
        non_zero_d = D_sp.n_nonzero;
    }
    int zero_d = (spp.dMat.n_elem - non_zero_d)/2;
    umat pert_record(zero_d,3); // record the time of edge allocation for analysis

    if (!stepwise) { // for single perturbuation, half of all zeros set to non-zero
        edges = floor(zero_d/1);
        cout << "edges reset to " << edges << endl;
    }

    std::size_t pos1 = bMatFileName.find_last_of("/");
    string bMatDir = bMatFileName.substr(0,pos1);
    pos1 = bMatFileName.find_last_of(")");
    std::size_t pos2 = bMatFileName.find_last_of(".");
    string newFolder = bMatFileName.substr(pos1+1, pos2-pos1-1);
    string Bpath;
    if (stepwise) {
        Bpath = bMatDir + "/" + newFolder + "/dkdt=" + to_string(edges);
    } else {
        Bpath = bMatDir + "/" + newFolder + "/delta_e=" + to_string(edges);
    }

    if (!exists(Bpath)) { // make directory if doesn't currently exist
        create_directories(Bpath);
    }

    unsigned int t_interval=1; // counter number of relaxation events
    unsigned int t_step=0; // counter number of algorithm iterations
    unsigned int edge_count=0; // counter number of edges added per perturbation

    while(!min(row_min)) { // run until graph is complete
        // select random node with degree < N-1
        uvec node = find(row_min == 0);
        node = shuffle(node);

        // select random missing edge of focal node
        uvec edge = find(spp.dMat.row(node(0))==0);
        edge = shuffle(edge);

        // add edges to adjacency matrix (-1 to distinguish between original edges)
        spp.topo.adjMat(node(0), edge(0)) = -1.0; // new edges in adjacency matrix allocated -1
        spp.topo.adjMat(edge(0), node(0)) = -1.0;

        // regenerate dispersal matrix - new normalization due to altered node degree, new edges not distance weighted
        spp.genDispMat();

        // record edge allocation
        pert_record(t_step, 0) = node(0);
        pert_record(t_step, 1) = edge(0);
        pert_record(t_step, 2) = t_interval;
        t_step++;

        // update row min vector
        row_min = min(abs(spp.dMat),1);

        edge_count++;

        cout << "edge count " << edge_count << endl;

        if (edge_count == edges) {
            if (stepwise) {
                printf("%d new edges allocated\n", edges * (t_interval + 1));
                metaCDynamics(tMax);
                string Bfile = Bpath + "/bMat_d" + to_string(t_interval) + ".mat";
                spp.xMat.save(Bfile, raw_ascii);
                string filenameRR = Bpath + "/pert_record.mat";
                pert_record.save(filenameRR, raw_ascii);
                t_interval++;
                edge_count = 0;

                { // switch off connectivity experiment for testing
                    if (t_interval > 100) {
                        row_min.ones();
                    }
                }

            } else {
                cout << "single connectivity perturbation" << endl;

                pert_record.resize(edge_count,pert_record.n_cols);

                int res = 1;

                // create CommunityDynamics instance for high resolution output
                CommunityDynamics dynamics;  // ODE dynamical object

                dynamics.xMat = &spp.xMat;

                if (spp.cMat.n_rows != 0) {
                    dynamics.cMat = &spp.cMat;
                }

                dynamics.rMat = &spp.rMat; // update rMat seen by CommmunityDynamics
                dynamics.dMat = &spp.dMat; // if dynamics.dMat is initialized, single domain relaxation selected
                dynamics.rho = &spp.rho;

                for (int t=0; t<tMax; t++) { // update every time step, save every res timesteps
                    printf("\rt = %d", t);
                    fflush(stdout);
                    ODE_state state(& dynamics);
                    state.integrate_until(t); // simulate single timestep
                    if ((t > 0) && (t % res == 0)) {
                        string Bfile = Bpath + "/bMat_d" + to_string(t) + ".mat";
                        spp.xMat.save(Bfile, raw_ascii);
                    }
                }
                row_min.ones(); // break while loop
            }
        }
    }

    if (edge_count) { // final relaxation unless no_nodes % edges = 0
        metaCDynamics(tMax);
        string Bfile = Bpath + "/bMat_d" + to_string(t_interval) + ".mat";
        spp.xMat.save(Bfile, raw_ascii);
        pert_record.print("\nRecord");
        string filenameRR = Bpath + "/pert_record.mat";
        pert_record.save(filenameRR, raw_ascii);
    }

    pert_record.print("\nRecord");
    string filenameRR = Bpath + "/pert_record.mat";
    pert_record.save(filenameRR, raw_ascii);
}

void Metacommunity::nodeRemoval(int tMax, int no_removals) {

    // summary:
    // remove nodes either randomly (if no_removals == 0) or sequentially
    // simulate dynamics and record biodiversity outcomes

    // required members:
    // spp - object of class Species, where model state is stored

    // external function calls:
    // metaCDynamics()

    // output:
    // stores biomass matrices

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// Generate output file names ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    storeTraj=3;
    if ((spp.topo.scVec.n_rows > 0) & (spp.topo.scVec_prime.n_rows > 0)) {
        cout << "\nLocal interaction matrices scaled" << endl;
    }
    std::size_t pos1 = bMatFileName.find_last_of("/");
    string bMatDir = bMatFileName.substr(0,pos1);
    pos1 = bMatFileName.find_last_of(")");
    std::size_t pos2 = bMatFileName.find_last_of(".");
    string newFolder = bMatFileName.substr(pos1+1, pos2-pos1-1);
    string Bpath;

    Bpath = bMatDir + "/" + newFolder + "/nodeRemoval";

    if (!exists(Bpath)) { // make directory if doesn't currently exist
        create_directories(Bpath);
    }

    uvec node_id;
    node_id = linspace<uvec>(no_removals,spp.topo.no_nodes-1,spp.topo.no_nodes-no_removals);

    bool write_to_file = true;

    if (write_to_file) {

        Species spp_copy = spp;
        for (int x = 0; x < node_id.n_rows; x++) {
            cout << "\nNode " << x << endl;
            spp = spp_copy;
            spp.xMat.shed_col(node_id(x));
            if (spp.rMat.n_rows > 0) {
                spp.rMat.shed_col(node_id(x));
            }
            spp.topo.network.shed_row(node_id(x));
            spp.topo.distMat.shed_row(node_id(x));
            spp.topo.distMat.shed_col(node_id(x));
            spp.dMat.shed_row(node_id(x)); // corresponding edges also removed
            spp.dMat.shed_col(node_id(x)); // corresponding edges also removed
            if (spp.topo.scVec.n_rows > 0) {
                spp.topo.scVec.shed_row(node_id(x));
            }
            if (spp.topo.scVec_prime.n_rows > 0) {
                spp.topo.scVec_prime.shed_row(node_id(x));
            }
            if (spp.topo.envMat.n_cols > 0) {
                spp.topo.envMat.shed_col(node_id(x));
            }
            spp.topo.no_nodes--;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Relax metacommunty and save new state ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            CommunityDynamics dynamics;  // ODE dynamical object

            dynamics.xMat = &spp.xMat;
            dynamics.rMat = &spp.rMat;
            if (spp.cMat.n_rows != 0) {
                dynamics.cMat = &spp.cMat;
            }

            dynamics.rho = &spp.rho;
            if (spp.topo.scVec.n_rows != 0) {
                dynamics.scVec = &spp.topo.scVec;
                dynamics.scVec_prime = &spp.topo.scVec_prime;
            }
            dynamics.S_p = spp.S_p;
            dynamics.S_c = spp.S_c;

            string Bfile = Bpath + "/bMat_nr_init" + to_string(node_id(x)) + ".mat";
            spp.xMat.save(Bfile, raw_ascii);
            if (spp.S_c > 0) {
                string Bfile_c = Bpath + "/bMat_nr_c_init" + to_string(node_id(x)) + ".mat";
                mat B_c = spp.xMat.rows(S_p,spp.xMat.n_rows-1);
                B_c.save(Bfile_c, raw_ascii);
            }
            metaCDynamics(tMax);
            Bfile = Bpath + "/bMat_nr" + to_string(node_id(x)) + ".mat";
            spp.xMat.save(Bfile, raw_ascii);
            if (spp.S_c > 0) {
                string Bfile_c = Bpath + "/bMat_nr_c" + to_string(node_id(x)) + ".mat";
                mat B_c = spp.xMat.rows(S_p,spp.xMat.n_rows-1);
                B_c.save(Bfile_c, raw_ascii);
            }
        }
    }
}

void Metacommunity::genJacobian() {

    // summary:
    // generate the numerical approximation of the Jacobian matrix for computing regional competitive overlap matrix

    // required members:
    // spp - object of class Species, where model state is stored

    // external function calls:
    // ODE.h::numerical_jacobian()

    // output:
    // jacobian - object for storing output of numerical_jacobian()

//    cout << "\nGenerating numerical Jacobian..." << endl;
    CommunityDynamics dynamics;  // ODE dynamical object
    dynamics.xMat = &spp.xMat;
    dynamics.indices_DP = &spp.indices_DP;
    dynamics.indices_S = &spp.indices_S;
    dynamics.rMat = &spp.rMat;
    if (spp.cMat.n_rows != 0) {
        dynamics.cMat = &spp.cMat;
    }
    dynamics.rho = &spp.rho;
    dynamics.bodymass = &spp.bodymass;
    dynamics.bodymass_inv = &spp.bodymass_inv;
    dynamics.mu = &spp.mu;
    if (spp.topo.scVec.n_rows != 0) {
        dynamics.scVec = &spp.topo.scVec;
        dynamics.scVec_prime = &spp.topo.scVec_prime;
    }
    dynamics.S_p = spp.S_p;
    dynamics.S_c = spp.S_c;

    // full dynamics - dispersal switched on
    if (spp.emMat.n_rows != 0) {
        dynamics.emMat = &spp.emMat;
    }
    dynamics.dMat = &spp.dMat;

    {
        ODE_state state(&dynamics);
        for (double t = 0; t < 0; t++) { // accesses back-end CVode initialization machinery
            state.integrate_until(t);
        }
        jacobian.set_size(dynamics.number_of_variables(),
                          dynamics.number_of_variables() );

        // Alias jacobian with a data type that numerical_Jacobian can handle:
        std::vector<double *> Jacobian_pointers(dynamics.number_of_variables());
        for(int i=dynamics.number_of_variables(); i-->0;){
//            REPORT(int((&jacobian(1,i))-(&jacobian(0,i))));
//            FATAL_ERROR("TEST: If the output above equals 1, all is fine, you can remove this and the previous line, otherwise ask Axel");
            Jacobian_pointers[i] = &(jacobian(0,i));
        }

        dynamics.numerical_Jacobian(Jacobian_pointers);
    }
//    cout << "\nFinished" << endl;
}

#define NORMALISE_REG_MAT 1

void Metacommunity::genCMatReg(double h) {

    // summary:
        // generate a numerically approximated regional interaction matrix via a computation harvesting experiment
        // for detail of algorithm see O'Sullivan et al. (2019)

    // arguments:
        // h - harvesting rate

    // required members:
        // jacobian
        // spp - object of class Species, where model state is stored

    // external function calls:
        // genJacobian()

    // output:
        // cMat_reg - spatially unresolved approximation of the metacommunity interaction matrix

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Compute numerical Jacobian and invert /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    genJacobian();
    mat J_inv = inv(jacobian);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Initialize population indices and matrix objects /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int S_tot = spp.xMat.n_rows;
    uvec index = linspace<uvec>(0, spp.topo.no_nodes-1, spp.topo.no_nodes);
    index *= S_tot; // indexes all populations of focal species i=1; indices for species j generated by element-wise addition
    cMat_reg.set_size(S_tot,S_tot); //
    vec H_i; // storage for harvesting vector
    vec dB_jx; // storage for perturbed, vectorized state matrix

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Iteratively harvest each species i and compute dB_jx ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int i=0; i<S_tot; i++) {
//        cout << "Harvesting spp " << i << endl;
        H_i.zeros(S_tot*spp.topo.no_nodes);
        H_i.elem(index + i) = h * spp.xMat.row(i);
        // "local shift in biomasses of species j due to harvesting of species i per unit h"
        dB_jx = -1 * (J_inv * H_i) / h; // -1 due to definition of competitive interaction coefficent as (+)c_ij
        mat dB_jx_Mat(&dB_jx(0), S_tot, spp.topo.no_nodes); // state vector in matrix form (for row sum operation)
        cMat_reg.col(i) = sum(dB_jx_Mat, 1); // store dB_j = sum_x(dB_jx)
    }

    jacobian.reset();
    cMat_reg = inv(cMat_reg); // invert

#if NORMALISE_REG_MAT
    mat norm, sgn;
    norm.zeros(cMat_reg.n_rows, cMat_reg.n_cols);
    sgn.zeros(cMat_reg.n_rows, cMat_reg.n_cols);
    norm.diag() = 1/sqrt(abs(cMat_reg.diag())); // normalization
    sgn.diag() = cMat_reg.diag() / abs(cMat_reg.diag());
    cMat_reg = sgn * norm * cMat_reg * norm; // normalize regional competitive overlap matrix
#endif
}

void Metacommunity::genSourceSink(int tFullRelax) {

    // summary:
        // infers which populations are dependent upon immigration for local detectability by switching off dispersal and relaxing to equilibrium

    // arguments:
        // tFullRelax - (long) relaxation time

    // required members:
        // spp - object of class Species, where model state is stored

    // external function calls:
        // metaCDynamics()

    // output:
        // matrices of dimensions SxN with 1 indicating source, -1 sink populations

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Set off-diagonal elements of dMat to zero and relax ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    mat DStore = spp.dMat;
    mat BStore_p = spp.xMat;
    mat Src_p, Snk_p;
    spp.dMat.zeros();
    spp.dMat.diag() = DStore.diag();
    metaCDynamics(tFullRelax); // relax

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Generate discrete source (1) sink (-1) matrices //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Src_p.zeros(spp.xMat.n_rows, spp.xMat.n_cols);
    Src_p(find(spp.xMat > spp.thresh)).ones(); // 1 allocated to detected populations after dispersal switched off

    spp.xMat = BStore_p; // reset metacommunity objects
    spp.dMat = DStore;

    Snk_p.zeros(spp.xMat.n_rows,spp.xMat.n_cols);

    Snk_p(find(spp.xMat > spp.thresh)).ones(); // 1 allocated to detected populations before dispersal switched off
    Snk_p = Src_p - Snk_p; // -1 allocated to sink populations, 1 allocated to pops excluded by immigration
    Snk_p(find(Snk_p == 1)).zeros(); // weak source populations removed

//    spp.xMat_src = Src_p + Snk_p; // final matrix source (1) sink (-1)
    spp.xMat_src = spp.xMat;
}

void Metacommunity::printParams() {

    // summary:
        // print model parameterization to console

    printf("\nModel parameters:\n");
    printf("\ntMax %d", tMax);
    printf("\nparOut %f", parOut);
    cout << "\nexperiment " << experiment;
    printf("\nrep %d", rep);
    printf("\nform_of_dynamics %d", g_form_of_dynamics);
    printf("\nsimTime %f", simTime);
    cout << "\ndate " << date;
    printf("\n\nc1 %f", spp.c1);
    printf("\nc2 %f", spp.c2 );
    printf("\nemRate %f", spp.emRate);
    printf("\ndispL %f", spp.dispL);
    printf("\npProducer %f", spp.pProducer) ;
    printf("\nalpha %e", spp.alpha);
    printf("\nsigma %f", spp.sigma);
    printf("\nrho %f", spp.rho);
    printf("\ncomp_dist %d", spp.comp_dist);
    printf("\nomega %f", spp.omega);
    printf("\n\nno_nodes %d", spp.topo.no_nodes);
    printf("\nphi %f", spp.topo.phi);
    if (spp.topo.skVec.n_rows > 0) {
        printf("\nsk %f", spp.topo.skVec(0));
    }
    printf("\nenvVar %d", spp.topo.envVar);
    printf("\nT_int %f", spp.topo.T_int);
    printf("\nvar_e %f", spp.topo.var_e);
    printf("\nrandGraph %d", spp.topo.randGraph);
    printf("\ngabriel %d", spp.topo.gabriel);
    printf("\ninvMax %d", invMax);
    printf("\ndispNorm %d\n\n", spp.dispNorm);
}

void Metacommunity::writePars(string filenameP) {

    boost::filesystem::ofstream params;
    params.open(filenameP);
    params << "date " << date << endl;
    params << "experiment " << experiment << endl;
    params << "form_of_dynamics " << g_form_of_dynamics << endl;
    params <<  "invMax " << invMax << endl;
    params << "parOut " << parOut << endl;
    params << "rep " << rep << endl;
    params << "simTime " << simTime << endl;
    params << "no_nodes " << spp.topo.no_nodes << endl;
    params << "randGraph " << spp.topo.randGraph << endl;
    params << "phi " << spp.topo.phi << endl;
    params << "envVar " << spp.topo.envVar << endl;
    params << "c1 " << spp.c1 << endl;
    params << "c2 " << spp.c2 << endl;
    if (spp.c3 > 0) {
        params << "c3 " << spp.c3 << endl;
    }
    params << "emRate " << spp.emRate << endl;
    params << "dispL " << spp.dispL << endl;
    params << "invasion " << spp.invasion << endl;
    params << "tMax " << tMax << endl;
    params << "pProducer " << spp.pProducer << endl;
    params << "prodComp " << spp.prodComp << endl;
    params << "var_e " << spp.topo.var_e << endl;
    params << "alpha " << spp.alpha << endl;
    params << "sigma " << spp.sigma << endl;
    params << "rho " << spp.rho << endl;
    params << "comp_dist " << spp.comp_dist << endl;
    if (spp.topo.T_int > 0) {
        params << "T_int " << spp.topo.T_int << endl;
    }
    params << "omega " << spp.omega << endl;
    if (spp.topo.skVec.n_rows > 0) {
        params << "sk " << spp.topo.skVec(0) << endl;
    } else {
        params << "sk " << "NULL" << endl;
    }
    params << "delta_g " << spp.delta_g << endl;
    params << "sigma_r " << spp.sigma_r << endl;
    params << "symComp " << spp.symComp << endl;
    params.close();

}

void Metacommunity::saveVideo(string path) {

    genCMatReg();
    genSourceSink();
    string Cfile = path + "cMat_reg_" + to_string(spp.invasion) + ".mat";
    cMat_reg.save(Cfile,raw_ascii);
    string Rfile = path + "rMat_" + to_string(spp.invasion) + ".mat";
    spp.rMat.save(Rfile,raw_ascii);
    string Xfile = path + "bMat_" + to_string(spp.invasion) + ".mat";
    spp.xMat.save(Xfile,raw_ascii);
    string XSfile = path + "bMat_src_" + to_string(spp.invasion) + ".mat";
    spp.xMat_src.save(XSfile,raw_ascii);
    cMat_reg.reset();
    spp.xMat_src.reset();
    if (store_params) {
        string Pfile = path + "params.mat";
        writePars(Pfile);
        string Nfile = path + "network.mat";
        spp.topo.network.save(Nfile, raw_ascii);
        string Efile = path + "envMat.mat";
        spp.topo.envMat.save(Efile, raw_ascii);
        string Dfile = path + "dMat.mat";
        spp.dMat.save(Dfile, raw_ascii);
        store_params = false;
    }
}

void Metacommunity::saveMC() {

    // summary:
        // generate file names and save model matrices to file

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Generate strings for file names //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    cout << "\nOutputting data... ";

    string filenameA, filenameB, filenameBC, filenameBs, filenameBCs, filenameC,
           filenameD, filenameE, filenameEM, filenameEMC, filenameEf, filenameI,
           filenameIP, filenameIPp, filenameIPc, filenameF, filenameN, filenameP,
           filenameR, filenameS, filenameSR, filenameT, filenameTr;

    if (bMatFileName.length() == 0) { // generate path from scratch
        ostringstream nameA, nameB, nameBC, nameBs, nameBCs, nameC, nameD,
        nameE, nameEM, nameEMC, nameEf, nameF, nameI, nameIP, nameIPp,
        nameIPc, nameN, nameP, nameR, nameS, nameSR, nameT, nameTr;

        ostringstream p;
        p << "SimulationData/N=" << spp.topo.no_nodes << "/" << experiment
          << "_experiment/" << date << "/";
        outputDirectory += p.str();

        if (!exists(outputDirectory)) { // make directory if doesn't currently exist
            create_directories(outputDirectory);
        }

        nameR << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "rMat" << rep << ".mat";
        filenameR = nameR.str();
        nameSR << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "sMat" << rep << ".mat";
        filenameSR = nameSR.str();
        nameA << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "aMat" << rep << ".mat";
        filenameA = nameA.str();
        nameB << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "bMat" << rep << ".mat";
        filenameB = nameB.str();
        nameBC << outputDirectory << date << "_" << experiment << "("
               << invMax << ")" << parOut << "bMat_c" << rep << ".mat";
        filenameBC = nameBC.str();
        nameP << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "params" << rep << ".mat";
        filenameP = nameP.str();
        nameN << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "network" << rep << ".mat";
        filenameN = nameN.str();
        nameTr << outputDirectory << date << "_" << experiment << "("
               << invMax << ")" << parOut << "trajec" << rep << ".mat";
        filenameTr = nameTr.str();
        nameEf << outputDirectory << date << "_" << experiment << "("
               << invMax << ")" << parOut << "envFluct" << rep << ".mat";
        filenameEf = nameEf.str();
        nameS << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "S" << rep << ".mat";
        filenameS = nameS.str();
        nameD << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "dMat" << rep << ".mat";
        filenameD = nameD.str();
        nameEM << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "emMat" << rep << ".mat";
        filenameEM = nameEM.str();
        nameEMC << outputDirectory << date << "_" << experiment << "("
               << invMax << ")" << parOut << "emMat_c" << rep << ".mat";
        filenameEMC = nameEMC.str();
        nameI << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "cMat" << rep << ".mat";
        filenameI = nameI.str();
        nameC << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "cMat_reg" << rep << ".mat";
        filenameC = nameC.str();
        nameE << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "envMat" << rep << ".mat";
        filenameE = nameE.str();
        nameT << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "tMat" << rep << ".mat";
        filenameT = nameT.str();
        filenameF = nameF.str();
        nameIP << outputDirectory << date << "_" << experiment << "("
                << invMax << ")" << parOut << "invProb" << rep << ".mat";
        filenameIP = nameIP.str();
        nameIPp << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "invProb_p" << rep << ".mat";
        filenameIPp = nameIPp.str();
        nameIPc << outputDirectory << date << "_" << experiment << "("
                << invMax << ")" << parOut << "invProb_c" << rep << ".mat";
        filenameIPc = nameIPc.str();
        nameBs << outputDirectory << date << "_" << experiment << "("
                << invMax << ")" << parOut << "bMat_src" << rep << ".mat";
        filenameBs = nameBs.str();
        nameBCs << outputDirectory << date << "_" << experiment << "("
                << invMax << ")" << parOut << "bMat_c_src" << rep << ".mat";
        filenameBCs = nameBCs.str();
        bMatFileName = filenameB;
        cout << "Biomass matrix file name: ";
        cout << bMatFileName << endl;
    } else { // use find and replace to generate paths from bMatFileName

        size_t b = bMatFileName.find("bMat");
        filenameBC = bMatFileName;
        filenameBC.replace(b, string("bMat").length(), "bMat_c");
        filenameD = bMatFileName;
        filenameD.replace(b, string("bMat").length(), "dMat");
        filenameEM = bMatFileName;
        filenameEM.replace(b, string("bMat").length(), "emMat");
        filenameEMC = bMatFileName;
        filenameEMC.replace(b, string("bMat").length(), "emMat_c");
        filenameI = bMatFileName;
        filenameI.replace(b, string("bMat").length(), "cMat");
        filenameC = bMatFileName;
        filenameC.replace(b, string("bMat").length(), "cMat_reg");
        filenameA = bMatFileName;
        filenameA.replace(b, string("bMat").length(), "aMat");
        filenameN = bMatFileName;
        filenameN.replace(b, string("bMat").length(), "network");
        filenameP = bMatFileName;
        filenameP.replace(b, string("bMat").length(), "params");
        filenameR = bMatFileName;
        filenameR.replace(b, string("bMat").length(), "rMat");
        filenameSR = bMatFileName;
        filenameSR.replace(b, string("bMat").length(), "sMat");
        filenameS = bMatFileName;
        filenameS.replace(b, string("bMat").length(), "S");
        filenameT = bMatFileName;
        filenameT.replace(b, string("bMat").length(), "tMat");
        filenameE = bMatFileName;
        filenameE.replace(b, string("bMat").length(), "envMat");
        filenameTr = bMatFileName;
        filenameTr.replace(b, string("bMat").length(), "trajec");
        filenameEf = bMatFileName;
        filenameEf.replace(b, string("bMat").length(), "envFluct");
        filenameIP = bMatFileName;
        filenameIP.replace(b, string("bMat").length(), "invProb");
        filenameIPp = bMatFileName;
        filenameIPp.replace(b, string("bMat").length(), "invProb_p");
        filenameIPc = bMatFileName;
        filenameIPc.replace(b, string("bMat").length(), "invProb_c");
        filenameBs = bMatFileName;
        filenameBs.replace(b, string("bMat").length(), "bMat_src");
        filenameBCs = bMatFileName;
        filenameBCs.replace(b, string("bMat").length(), "bMat_c_src");
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Store matrix objects ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Note, if specific experiment run, only relevant objects are output
    if (spp.trajectories.n_rows != 0) { // output trajectories ONLY
        spp.trajectories.save(filenameTr, raw_ascii);
        if (spp.efMat.n_rows != 0) {
            mat R = spp.rMat;
            R.reshape(1, R.size());
            rowvec N(R.n_cols);
            N.fill(datum::nan);
            R = join_vert(R, N);
            spp.fluctuations = join_vert(R, spp.fluctuations);
            spp.fluctuations.save(filenameEf, raw_ascii);
        }
    } else if (cMat_reg.n_rows != 0) { // output regional competitive overlap matrix ONLY
        cMat_reg.save(filenameC, raw_ascii);
    } else if (spp.xMat_src.n_rows != 0) { // output source-sink matrices ONLY
        mat B_p_src = spp.xMat_src.rows(0,spp.rMat.n_rows-1);
        B_p_src.save(filenameBs, raw_ascii);
        if (spp.xMat.n_rows > spp.rMat.n_rows) {
            mat B_c_src = spp.xMat_src.rows(spp.rMat.n_rows, spp.xMat.n_rows-1);
            B_c_src.save(filenameBCs, raw_ascii);
        }
    } else { // output standard model objects
        spp.topo.network.save(filenameN, raw_ascii);
        spp.rMat.save(filenameR, raw_ascii);
        if (spp.sMat.n_rows != 0) {
            spp.sMat.save(filenameSR, raw_ascii);
        }
        if (bMatFileName.length() == 0) {
            mat B_p = spp.xMat.rows(0,spp.S_p-1);
            B_p.save(filenameB, raw_ascii);
        } else {
            mat B_p = spp.xMat.rows(0,spp.S_p-1);
            B_p.save(bMatFileName, raw_ascii);
        }
        if (spp.emMat.n_rows != 0) {
            mat Em_p = spp.emMat.rows(0,spp.rMat.n_rows-1);
            Em_p.save(filenameEM, raw_ascii);
            if (spp.xMat.n_rows > spp.rMat.n_rows) {
                mat Em_c = spp.emMat.rows(spp.rMat.n_rows, spp.xMat.n_rows-1);
                Em_c.save(filenameEMC, raw_ascii);
            }
        }

        if (spp.xMat.n_rows > spp.rMat.n_rows) {
            mat B_c = spp.xMat.rows(spp.rMat.n_rows, spp.xMat.n_rows-1);
            B_c.save(filenameBC, raw_ascii);
        }

        if (spp.topo.envMat.n_rows != 0) {
            spp.topo.envMat.save(filenameE, raw_ascii);
            spp.tMat.save(filenameT, raw_ascii);
        }
        spp.sppRichness.save(filenameS, raw_ascii);
        spp.dMat.save(filenameD, raw_ascii);
        if (spp.S_p > 0) {
            mat C = spp.cMat.submat(0,0,spp.S_p-1,spp.S_p-1);
            C.save(filenameI, raw_ascii);
        }
        if (spp.xMat.n_rows > spp.rMat.n_rows) {
            mat A = spp.cMat.submat(spp.S_p,0,spp.cMat.n_rows-1,spp.S_p-1);
            A.save(filenameA, raw_ascii);
        }
        writePars(filenameP);
    }
}

bool fexists(const std::string& filename) { // check if file/directory exists
    std::ifstream ifile(filename.c_str());
    return (bool)ifile;
}

void Metacommunity::loadMC(string bMatFile) {

    // summary:
        // import metacommunity model

    // arguments:
        // bMatFile - path/file name of produce biomass matrix from which the rest of the file names are generated

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Generate strings for file names //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bMatFileName = bMatFile;

    cout << "Importing data, file name " << bMatFileName << endl;

    string filenameA, filenameR, filenameSR, filenameBC, filenameC, filenameF, filenameP, filenameN,
            filenameT, filenameTr, filenameS, filenameD, filenameEM, filenameEMC,
            filenameI, filenameE, filenameIPp, filenameIPc;
    size_t b = bMatFileName.find("bMat");
    filenameBC = bMatFileName;
    filenameBC.replace(b, string("bMat").length(), "bMat_c");
    filenameD = bMatFileName;
    filenameD.replace(b, string("bMat").length(), "dMat");
    filenameEM = bMatFileName;
    filenameEM.replace(b, string("bMat").length(), "emMat");
    filenameEMC = bMatFileName;
    filenameEMC.replace(b, string("bMat").length(), "emMat_c");
    filenameI = bMatFileName;
    filenameI.replace(b, string("bMat").length(), "cMat");
    filenameA = bMatFileName;
    filenameA.replace(b, string("bMat").length(), "aMat");
    filenameN = bMatFileName;
    filenameN.replace(b, string("bMat").length(), "network");
    filenameP = bMatFileName;
    filenameP.replace(b, string("bMat").length(), "params");
    filenameR = bMatFileName;
    filenameR.replace(b, string("bMat").length(), "rMat");
    filenameSR = bMatFileName;
    filenameSR.replace(b, string("bMat").length(), "sMat");
    filenameS = bMatFileName;
    filenameS.replace(b, string("bMat").length(), "S");
    filenameT = bMatFileName;
    filenameT.replace(b, string("bMat").length(), "tMat");
    filenameE = bMatFileName;
    filenameE.replace(b, string("bMat").length(), "envMat");
    filenameTr = bMatFileName;
    filenameTr.replace(b, string("bMat").length(), "trajec");
    filenameIPp = bMatFileName;
    filenameIPp.replace(b, string("bMat").length(), "invProb_p");
    filenameIPc = bMatFileName;
    filenameIPc.replace(b, string("bMat").length(), "invProb_c");

    cout << "Biomass matrix file name:\n" << bMatFileName << endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Extract data stored in param file ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int number_of_params = 0;
    std::string paramLine;
    boost::filesystem::ifstream countParams;
    countParams.open(filenameP);

    while (std::getline(countParams, paramLine))
        ++number_of_params;

    boost::filesystem::ifstream param(filenameP);

    string parName[number_of_params];
    string pars[number_of_params];
    if (param.is_open())
    {
        for (int i=0; i<number_of_params; i++) {
            param >> parName[i] >> pars[i];
        }
    }

    for(int i=0; i<number_of_params; i++) {
        if (parName[i] == "date") {
            date = pars[i];
            continue;
        } else if (parName[i] == "experiment") {
            experiment = pars[i];
            continue;
        } else if (parName[i] == "form_of_dynamics") {
            g_form_of_dynamics = stoi(pars[i]);
            continue;
        } else if (parName[i] == "invMax") {
            invMax = stoi(pars[i]);
            continue;
        } else if (parName[i] == "parOut") {
            parOut = stof(pars[i]);
            continue;
        } else if (parName[i] == "rep") {
            rep = stoi(pars[i]);
            continue;
        } else if (parName[i] == "simTime") {
            simTime += stof(pars[i]);
            continue;
        } else if (parName[i] == "no_nodes"){
            spp.topo.no_nodes = stoi(pars[i]);
            continue;
        } else if (parName[i] == "randGraph") {
            spp.topo.randGraph = stoi(pars[i]);
            continue;
        } else if (parName[i] == "phi") {
            spp.topo.phi = stof(pars[i]);
            continue;
        } else if (parName[i] == "envVar") {
            spp.topo.envVar = stoi(pars[i]);
            continue;
        } else if (parName[i] == "c1") {
            spp.c1 = stof(pars[i]);
            continue;
        } else if (parName[i] == "c2") {
            spp.c2 = stof(pars[i]);
            continue;
        } else if (parName[i] == "c3") {
            spp.c3 = stof(pars[i]);
            continue;
        } else if (parName[i] == "emRate") {
            spp.emRate = stof(pars[i]);
            continue;
        } else if (parName[i] == "dispL") {
            spp.dispL = stof(pars[i]);
            continue;
        } else if (parName[i] == "invasion") {
            spp.invasion = stoi(pars[i]);
            spp.invEvent = stoi(pars[i]);
            spp.inv0 = spp.invasion;
            continue;
        } else if (parName[i] == "tMax") {
            tMax = stoi(pars[i]);
            continue;
        } else if (parName[i] == "pProducer") {
            spp.pProducer = stof(pars[i]);
            continue;
        } else if (parName[i] == "pProducer") {
            if (stoi(pars[i])) {
                spp.prodComp = true;
            } else {
                spp.prodComp = false;
            }
            continue;
        } else if (parName[i] == "var_e") {
            spp.topo.var_e = stof(pars[i]);
            continue;
        } else if (parName[i] == "alpha") {
            spp.alpha = stof(pars[i]);
            continue;
        } else if (parName[i] == "sigma") {
            spp.sigma = stof(pars[i]);
            continue;
        } else if (parName[i] == "rho") {
            spp.rho = stof(pars[i]);
            continue;
        } else if (parName[i] == "comp_dist") {
            spp.comp_dist = stoi(pars[i]);
            continue;
        } else if (parName[i] == "T_int") {
            spp.topo.T_int = stof(pars[i]);
            continue;
        } else if (parName[i] == "omega") {
            spp.omega = stof(pars[i]);
            continue;
        } else if (parName[i] == "sk") {
            spp.topo.skVec.set_size(1);
            if ((pars[i] == "import") || (pars[i] == "NULL")){
                spp.topo.skVec(0) = 0.0;
            } else {
                spp.topo.skVec(0) = stof(pars[i]);
            }
            continue;
        } else if (parName[i] == "delta_g") {
            spp.delta_g = stof(pars[i]);
            continue;
        } else if (parName[i] == "sigma_r") {
            spp.sigma_r = stof(pars[i]);
            continue;
        } else if (parName[i] == "symComp") {
            if (stoi(pars[i])) {
                spp.symComp = true;
            } else {
                spp.symComp = false;
            }
            continue;
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Load matrix objects ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    spp.xMat.load(bMatFileName);
    spp.S_p = spp.xMat.n_rows;

    spp.dMat.load(filenameD);
    if(fexists(filenameI)) {
        spp.cMat.load(filenameI);
    }
    if(fexists(filenameEM)) {
        spp.emMat.load(filenameEM);
    }
    spp.topo.network.load(filenameN);
    if (fexists(filenameT)) {
        spp.tMat.load(filenameT);
    }

    if (fexists(filenameE)) {
        spp.topo.envMat.load(filenameE);
    }
    spp.rMat.load(filenameR);
    if (fexists(filenameSR)) {
        spp.sMat.load(filenameSR);
    }
    spp.sppRichness.load(filenameS);

    cout << "\nImported network rows 0-4 = " << endl;
    cout << spp.topo.network.rows(0,min((int) spp.topo.network.n_rows-1, 4));
    spp.topo.genDistMat();
    printf("\nSpecies richness, S_p = %d, S_c = %d\n", (int) spp.rMat.n_rows, (int) spp.xMat.n_rows - (int) spp.rMat.n_rows);

    if (g_form_of_dynamics) {
        cout << "\nAllocating compartment indices..." << endl;
        spp.indices_DP = find(spp.xMat > -0.5); //splitting at -0.5 is ...
        spp.indices_S = find(spp.xMat <= -0.5); //safter than splitting at -1
    } else {
        cout << "\nLVMCM, all indices allocated to DP..." << endl;
        spp.indices_DP = linspace<uvec>(0, spp.xMat.n_elem-1, spp.xMat.n_elem);
    }
}

// Local Variables:
// c-file-style: "stroustrup"
// End:
