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

#include "Metacommunity.h"
#include "Topography.h"
#include "Species.h"
#include "ODE.h"
#include "CommunityDynamics.h"
#include "LVMCM_rng.h"

using namespace std;
using namespace arma;
using namespace boost::numeric::ublas;
using namespace boost::filesystem;

class integrating_dynamical_object : public ODE_dynamical_object {
    ODE_dynamical_object & main;
    vec integrals;
    integrating_dynamical_object();
public:
    int dynamics(ODE_vector const & state, 
		 ODE_vector & time_derivative) {
        int retval = main.dynamics(state,time_derivative);
        int n = main.number_of_variables();
        for(int i = 0; i < n; ++i){
            time_derivative[n+i]=state[i];
        }
        return(retval);
    }
    void write_state_to(ODE_vector & state) const {
        main.write_state_to(state);
        int n = main.number_of_variables();
        for(int i = 0; i < n; ++i){
            state[n+i]=integrals[i];
        }
    }
    void read_state_from(const ODE_vector & state) {
        main.read_state_from(state);
        int n = main.number_of_variables();
        for(int i = 0; i < n; ++i){
            integrals[i]=state[n+i];
        }
    }
    int number_of_variables() const {
        return(2*main.number_of_variables());
        }
        vec get_integrals() const {return integrals;}
        integrating_dynamical_object(ODE_dynamical_object & other) :
        main(other) , integrals(other.number_of_variables(),fill::zeros) {
        }
    virtual void prepare_for_integration(){
	    main.prepare_for_integration();
    }; 
    virtual void cleanup_after_integration(){
	    main.cleanup_after_integration();
    };
};

void Metacommunity::metaCDynamics(int relaxT, bool testing) {

    // summary:
        // numerically solve metacommunity dynamics for trajectory simulation (serial or parallel) or invader testing
        // trajectory approximation uses Sundials CVODE solver (see communityDynamics.h/.cpp)

    // arguments:
        // relaxT - relaxtion time (if != tMax)
        // dispersal - select with/without dispersal term (for invader testing)

    // required memebers:
        // sppPool - object of class Species, where model state is stored
        // tMax - relaxation time
        // storeTraj - selects whether trajectory object, dimensions (N*S)x(tMax), should be stored (large memory cost)
        // if parallel approximation used:
            // iterS - number Schwarz iterations
            // deltaT - size of Schwarz timewindow

    // external function calls:
        // ODE.h::integrate_until()

    // output:
        // updates to matrices bMat_p, and bMat_c

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Initialize community dynamics machinery ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CommunityDynamics dynamics;  // ODE dynamical object

    mat B_p, B_c, R, D; // storage for matrix subsets in case of parallel invader testing only

    dynamics.bMat_p = &sppPool.bMat_p;
    dynamics.rMat = &sppPool.rMat;
    if (sppPool.cMat.n_rows != 0) {
        dynamics.cMat = &sppPool.cMat;
    }
    dynamics.rho = &sppPool.rho;
    if (sppPool.topo.scVec.n_rows != 0) {
        dynamics.scVec = &sppPool.topo.scVec;
        dynamics.scVec_prime = &sppPool.topo.scVec_prime;
    }
    dynamics.S_p = sppPool.S_p;
    dynamics.S_c = sppPool.S_c;

    if (testing) {
        dynamics.testing = true; // invader testing - dispersal switched off for efficiency
    } else {
        // full dynamics - dispersal switched on
        if (sppPool.emMat_p.n_rows != 0) {
            dynamics.emMat_p = &sppPool.emMat_p;
        }
        if (sppPool.topo.network.n_rows == sppPool.topo.no_nodes) {
            // whole domain scale - dMat only set
            dynamics.dMat = &sppPool.dMat_n; // if dynamics.dMat is initialized, single domain relaxation selected
        } else {
            // parallel LV relaxtion, dMat_n/m and uMat set
            dynamics.uMat_p = &sppPool.uMat_p;
            dynamics.dMat_n = &sppPool.dMat_n;
            dynamics.dMat_m = &sppPool.dMat_m;
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Metacommunity relaxation step /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // This extends dynamics to compute also time integrals over
    // biomasses, used to compute averages.
    integrating_dynamical_object i_dynamics(dynamics);

    int res = 10; // set time step resolution for storing trajectories; relaxT must be multiple of res
    if (storeTraj == 0) { // standard relaxation
        {
            ODE_state state(& i_dynamics);
            state.integrate_until(relaxT);
        }

        if (sppPool.topo.network.n_rows < sppPool.topo.no_nodes) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Compute temporal averages, e.g. for domain decomposed solution //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//            sppPool.bavMat_p =
//            reshape(i_dynamics.get_integrals(),
//                sppPool.bMat_p.n_rows, sppPool.bMat_p.n_cols) /
//            relaxT;
	    // average resource biomasses
//	    sppPool.bavMat_c =
//		reshape(i_dynamics.get_integrals().
//			tail( dynamics.number_of_variables() - sppPool.bMat_p.n_rows * sppPool.bMat_p.n_cols ),
//			sppPool.bMat_c.n_rows, sppPool.bMat_c.n_cols) /
//		relaxT; // average consumer biomasses
	    }

    } else if (storeTraj == 1) { // relax, storing trajectory object
        {
            // write time point to file in SxN biomass matrix
            // file directory reflects the parOut and rep of the current assembly

            std::size_t pos1 = bMatFileName.find_last_of("/");
            string bMatDir = bMatFileName.substr(0,pos1);
            pos1 = bMatFileName.find_last_of(")");
            std::size_t pos2 = bMatFileName.find_last_of(".");
            string newFolder = bMatFileName.substr(pos1+1, pos2-pos1-1);
            string Bpath = bMatDir + "/" + newFolder + "/trajectory";
            cout << "\nSaving matrices to " << Bpath << endl;

            if (!exists(Bpath)) { // make directory if doesn't currently exist
                create_directories(Bpath);
            }

            ODE_state state(& i_dynamics);
            for (int t = 0; t < relaxT/res; t++) {
                cout << t << endl;
                state.integrate_until(t*res);
                mat Btofile(&state[0], 1, dynamics.number_of_variables());
                Btofile.reshape(sppPool.S_p, sppPool.topo.no_nodes);
                string Bfile = Bpath + "/bMat" + to_string(t) + ".mat";
                Btofile.save(Bfile, raw_ascii);
            }

            // store parameter file in same directory for convenience
            string filenameP = Bpath + "/pars.mat";
            boost::filesystem::ofstream params;
            params.open(filenameP);
            params << "date " << date << endl;
            params << "jobID " << jobID << endl;
            params << "experiment " << experiment << endl;
            params <<  "invMax " << invMax << endl;
            params << "parOut " << parOut << endl;
            params << "rep " << rep << endl;
            params << "simTime " << simTime << endl;
            params << "no_nodes " << sppPool.topo.no_nodes << endl;
            params << "bisec " << sppPool.topo.bisec << endl;
            params << "phi " << sppPool.topo.phi << endl;
            params << "envVar " << sppPool.topo.envVar << endl;
            params << "c1 " << sppPool.c1 << endl;
            params << "c2 " << sppPool.c2 << endl;
            if (sppPool.c3 > 0) {
                params << "c3 " << sppPool.c3 << endl;
            }
            params << "emRate " << sppPool.emRate << endl;
            params << "dispL " << sppPool.dispL << endl;
            params << "invasion " << sppPool.invasion << endl;
            params << "iterS " << iterS << endl;
            params << "deltaT " << deltaT << endl;
            params << "tMax " << tMax << endl;
            params << "pProducer " << sppPool.pProducer << endl;
            params << "prodComp " << sppPool.prodComp << endl;
            params << "var_e " << sppPool.topo.var_e << endl;
            params << "alpha " << sppPool.alpha << endl;
            params << "sigma " << sppPool.sigma << endl;
            params << "rho " << sppPool.rho << endl;
            params << "comp_dist " << sppPool.comp_dist << endl;
            params << "T_int " << sppPool.topo.T_int << endl;
            params.close();
        }
    } else if (storeTraj == 2) { // concatenate previous trajectory object adding zero vectors for new invaders

        // This is used if a continuous trajectory *including invasions* is desired
        // If so this will need to be updated

        if (sppPool.trajectories.n_rows == 0) {
            ODE_state state(& i_dynamics);
            sppPool.trajectories.set_size(relaxT/res, dynamics.number_of_variables());
            for (int t = 0; t < relaxT; t++) {
                state.integrate_until(t*res);
                mat tmp(&state[0], 1, dynamics.number_of_variables());
                sppPool.trajectories.row(t) = tmp;
            }
        } else {
            ODE_state state(& i_dynamics);
            mat trajTemp; trajTemp.set_size(relaxT/res, dynamics.number_of_variables());
            for (int t = 0; t < relaxT/res; t++) {
                state.integrate_until(t*res);
                mat tmp(&state[0], 1, dynamics.number_of_variables());
                trajTemp.row(t) = tmp;
            }
            if (trajTemp.n_cols != sppPool.trajectories.n_cols) {
                mat zeroBiomass(sppPool.trajectories.n_rows, (trajTemp.n_cols - sppPool.trajectories.n_cols));
                zeroBiomass.fill(datum::nan);
                sppPool.trajectories = join_horiz(sppPool.trajectories, zeroBiomass);
                sppPool.trajectories = join_vert(sppPool.trajectories, trajTemp);
            }
        }
    } else if (storeTraj == 3) { // relax, printing after each res timesteps

        cout << "\nSimulating dynamics" << endl;

        ODE_state state(& i_dynamics);
        for (int t = 0; t < relaxT/res; t++) {
            cout << t << endl;
            state.integrate_until(t*res);
        }
    } else if (storeTraj < 0) {
        {
            res=100; // lower res to reduce storage costs
            if (bMatFileName.length() == 0) { // generate path from scratch
                ostringstream p, nameBF;
                p << "SimulationData/N=" << sppPool.topo.no_nodes << "/" << experiment
                  << "_experiment/" << date << "/";
                outputDirectory += p.str();

                if (!exists(outputDirectory)) { // make directory if doesn't currently exist
                    create_directories(outputDirectory);
                }

                nameBF << outputDirectory << date << "_" << experiment << "("
                       << invMax << ")" << parOut << "bMat" << rep << ".mat";
                string filenameBF = nameBF.str();
                bMatFileName = filenameBF;
            }

            // write time point to file in SxN biomass matrix
            // file directory reflects the parOut and rep of the current assembly
            std::size_t pos1 = bMatFileName.find_last_of("/");
            string bMatDir = bMatFileName.substr(0,pos1);
            pos1 = bMatFileName.find_last_of(")");
            std::size_t pos2 = bMatFileName.find_last_of(".");
            string newFolder = bMatFileName.substr(pos1+1, pos2-pos1-1);
            string Bpath = bMatDir + "/" + newFolder + "/trajectory_" + std::to_string(abs(storeTraj));
//            cout << "\nSaving matrices to " << Bpath << endl;

            if (!exists(Bpath)) { // make directory if doesn't currently exist
                create_directories(Bpath);
            }

            ODE_state state(& i_dynamics);
            for (int t = 0; t < relaxT/res; t++) {
//                cout << t << endl;
                state.integrate_until(t*res);
                mat Btofile(&state[0], 1, dynamics.number_of_variables());
                Btofile.reshape(sppPool.S_p, sppPool.topo.no_nodes);
                string Bfile = Bpath + "/bMat" + to_string(t) + ".mat";
                Btofile.save(Bfile, raw_ascii);
            }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Set off-diagonal elements of dMat_n to zero and relax ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            mat DStore = sppPool.dMat_n;
            mat BStore_p = sppPool.bMat_p;
            mat Src_p, Snk_p;
            sppPool.dMat_n.zeros();
            sppPool.dMat_n.diag() = DStore.diag();
            storeTraj=0;
            metaCDynamics(1000); // relax

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Generate discrete source (1) sink (-1) matrices //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            Src_p.zeros(sppPool.bMat_p.n_rows, sppPool.bMat_p.n_cols);
            Src_p(find(sppPool.bMat_p > sppPool.thresh)).ones(); // 1 allocated to detected populations after dispersal switched off

            sppPool.bMat_p = BStore_p; // reset metacommunity objects
            sppPool.dMat_n = DStore;

            Snk_p.zeros(sppPool.bMat_p.n_rows,sppPool.bMat_p.n_cols);

            Snk_p(find(sppPool.bMat_p > sppPool.thresh)).ones(); // 1 allocated to detected populations before dispersal switched off
            Snk_p = Src_p - Snk_p; // -1 allocated to sink populations, 1 allocated to pops excluded by immigration
            Snk_p(find(Snk_p == 1)).zeros(); // weak source populations removed

            string Bfile = Bpath + "/bMat_src.mat";
            Src_p.save(Bfile, raw_ascii);
        }
    }
}

uvec Metacommunity::invaderSample(int trophLev, int no_invaders) {
    // summary:
        // introduce new random species and test for positive growth rates in single evaluation of ODE

    // arguments:
        // trophLev - 0: producer; 1: consumer
        // no_invaders - number of successful invaders to introduce

    // required memebers:
        // sppPool - object of class Species, where model state is stored

    // external function calls:
        // Species::invade()

    // output:
        // vector of positively growing species indices
        // for MPI program, this vector is distributed to all parallel processes and species are selected according to
        // deterministic algorithm

    // testing parameters
    double min_b = 0; // minimum biomass for inclusion in model
    double inv = 1e-6; // invasion biomass
    uvec invaderIndex, posGrowth;

    if (trophLev == 0) { // sample producer species
        // size of testing pool
        int invExcess_p = 3; // invade excess species to account for difficulty in finding successful invader -- set arbitrarily

        for (int i = 0; i < invExcess_p * no_invaders; i++) {
            sppPool.invade(0); // invade a producer
        }

        // single evaluation of ODE to test for positive growth, invaders ONLY
        mat bInv;
        if (sppPool.prodComp) {
            bInv = (sppPool.rMat.rows(sppPool.S_p, sppPool.S_p + sppPool.I_p - 1) -
                    sppPool.cMat.submat(sppPool.S_p, 0,
                                        sppPool.S_p + sppPool.I_p - 1, sppPool.cMat.n_cols - 1) *
                    sppPool.bMat_p) %
                    sppPool.bMat_p.rows(sppPool.S_p, sppPool.S_p + sppPool.I_p - 1) +
                    sppPool.bMat_p.rows(sppPool.S_p, sppPool.S_p + sppPool.I_p - 1) * sppPool.dMat_n;
        } else {
            bInv = (sppPool.rMat.rows(sppPool.S_p, sppPool.S_p + sppPool.I_p - 1) - sppPool.bMat_p.rows(sppPool.S_p, sppPool.S_p + sppPool.I_p - 1)) %
                   sppPool.bMat_p.rows(sppPool.S_p, sppPool.S_p + sppPool.I_p - 1) +
                   sppPool.bMat_p.rows(sppPool.S_p, sppPool.S_p + sppPool.I_p - 1) * sppPool.dMat_n;
        }

        vec bInv_max(bInv.n_rows);
        for (int i = 0; i < bInv.n_rows; i++) { // store maximum local biomass of each invader
            bInv_max(i) = bInv.row(i).max();
        }

        sppPool.bMat_p.rows(sppPool.S_p, sppPool.S_p + sppPool.I_p - 1) = bInv; // save into bMat for assessment of maxima

        posGrowth = find(bInv_max >= min_b); // index vector of positively growing species

    } else if (trophLev == 1) { // sample consumer species

        // size of testing pool
        int invExcess_c = 3; // invade excess species to account for difficulty in finding successful invader

        for (int i = 0; i < invExcess_c * no_invaders; i++) {
            sppPool.invade(1); // invade a consumer
        }

        mat bInv = sppPool.rho * (sppPool.cMat.submat(sppPool.S_p + sppPool.I_p + sppPool.S_c, 0, sppPool.cMat.n_rows-1, sppPool.S_p + sppPool.I_p - 1) *
                                  sppPool.bMat_p.rows(0, sppPool.S_p + sppPool.I_p - 1) - 1) %
                   sppPool.bMat_p.rows(sppPool.S_p + sppPool.I_p + sppPool.S_c, sppPool.bMat_p.n_rows-1) +
                   sppPool.bMat_p.rows(sppPool.S_p + sppPool.I_p + sppPool.S_c, sppPool.bMat_p.n_rows-1) * sppPool.dMat_n;

        sppPool.bMat_p.rows(sppPool.S_p + sppPool.I_p + sppPool.S_c, sppPool.bMat_p.n_rows-1) = bInv; // save into bMat for assessment of maxima
        vec bInv_max(bInv.n_rows);
        for (int i = 0; i < bInv.n_rows; i++) { // store maximum local biomass of each invader
            bInv_max(i) = bInv.row(i).max();
        }
        posGrowth = find(bInv_max >= min_b);
    }
    return(posGrowth); // return index of successful invaders
}

mat Metacommunity::invaderCleanup(int trophLev, uvec posGrowth) {
    // summary:
        // remove unsuccessful or excess invaders from model

    // arguments:
        // trophLev - 0: producer; 1: consumer
        // posGrowth - indices species to retain, output of invaderSample
        // this is gathered from all MPI processes in main.cpp

    // required memebers:
        // sppPool - object of class Species, where model state is stored

    // external function calls:
        // Species::invade()

    // output:
        // vector of maximum biomasses of invaders used to select 'port patch'

    mat bInv_max; // record max b_ix for subdomain selection

    if (trophLev == 0) { // clean up producers

        umat negGrowth(sppPool.I_p,1);
        negGrowth.col(0) = linspace<uvec>(sppPool.S_p, sppPool.S_p + sppPool.I_p - 1,
                                          sppPool.I_p);

        for (int i = posGrowth.n_rows - 1; i>=0; i--) {
            negGrowth.shed_row(posGrowth(i));
        }

        for (int i = negGrowth.n_rows - 1; i >= 0; i--) {
            sppPool.bMat_p.shed_row(negGrowth(i)); // remove prod biomass vec
            if (sppPool.rMat.n_rows > 0) {
                sppPool.rMat.shed_row(negGrowth(i)); // remove prod growth vec
            }
            if (sppPool.sMat.n_rows > 0) {
                sppPool.sMat.shed_row(negGrowth(i)); // remove prod growth vec
            }
            sppPool.cMat.shed_row(negGrowth(i)); // remove prod comp term
            sppPool.cMat.shed_col(negGrowth(i));
            if (sppPool.topo.envVar != 0) {
                sppPool.tMat.shed_row(negGrowth(i));  // remove prod env tol vec
            }
            if (sppPool.emRate < 0) {
                sppPool.emMat_p.shed_row(negGrowth(i));  // remove prod emigration rate
            }
        }

        if (posGrowth.n_rows > 0) {
            bInv_max.set_size(posGrowth.n_rows,1);

            for (int i = 0; i < posGrowth.n_rows; i++) { // store maximum local biomass of each invader
                bInv_max.row(i) = sppPool.bMat_p.row(sppPool.S_p + i).max();
            }
        }

    } else if (trophLev == 1) { // clean up consumers

        umat negGrowth(sppPool.I_c,1);
        negGrowth.col(0) = linspace<uvec>(sppPool.S_p + sppPool.I_p + sppPool.S_c, sppPool.bMat_p.n_rows-1,
                                          sppPool.I_c);

        for (int i = posGrowth.n_rows - 1; i>=0; i--) {
            negGrowth.shed_row(posGrowth(i));
        }

        for (int i = negGrowth.n_rows - 1; i >= 0; i--) {
            sppPool.bMat_p.shed_row(negGrowth(i)); // remove cons biomass vec
            sppPool.cMat.shed_row(negGrowth(i)); // remove interaction coefficients
            sppPool.cMat.shed_col(negGrowth(i));
            if (sppPool.emRate < 0) {
                sppPool.emMat_p.shed_row(negGrowth(i));  // remove prod emigration rate
            }
        }

        if (posGrowth.n_rows > 0) {
            bInv_max.set_size(posGrowth.n_rows,1);

            for (int i = 0; i < posGrowth.n_rows; i++) { // store maximum local biomass of each invader
                bInv_max.row(i) = sppPool.bMat_p.row(sppPool.S_p + sppPool.I_p + i).max();
            }
        }
    }
    return(bInv_max);
}

void Metacommunity::invaderPopulate(int trophLev, mat bInv_max) {
    // summary:
        // reset biomass vectors for selected invaders
        // species introduced at low biomass into patch in which growth rate is maximum during testing

    // arguments:
        // trophLev - 0: producer; 1: consumer
        // bInv_max - maximum biomasses of invaders
        // this is gathered from all MPI processes in main.cpp

    // required memebers:
        // sppPool - object of class Species, where model state is stored

    double inv = 1e-6;

    ucolvec subDom_max, bMax_index(1);

    if (bInv_max.n_rows > 0) {
        // check which subdomain species are growing fastest
        subDom_max = index_max(bInv_max,1);
    }

    if (trophLev == 0) { // populate producer species
        if (bInv_max.n_rows > 0) {
            // reset b_i and invade into favoured node
            for (int i=0; i<subDom_max.n_rows; i++) {
                if (subDom_max(i) == sppPool.topo.subdomain) {
                    bMax_index = index_max(sppPool.bMat_p.row(sppPool.S_p+i));
                    sppPool.bMat_p.row(sppPool.S_p+i).zeros();
                    sppPool.bMat_p(sppPool.S_p+i, bMax_index(0)) = inv;
                } else {
                    sppPool.bMat_p.row(sppPool.S_p+i).zeros();
                }
            }

            if (sppPool.topo.bisec > 0) {
                if (sppPool.uMat_p.n_rows == 0) {
                    sppPool.uMat_p.set_size(sppPool.bMat_p.n_rows, sppPool.bMat_p.n_cols);
                    sppPool.uMat_p.zeros();
                } else {
                    sppPool.uMat_p.insert_rows(sppPool.S_p,sppPool.rMat.n_rows-sppPool.S_p);
                }
            }
        }
        sppPool.I_p = 0; // reset invader counter
        sppPool.S_p = sppPool.rMat.n_rows;
    } else if (trophLev == 1) { // populate consumer species
        if (bInv_max.n_rows > 0) {
            // reset b_i and invade into favoured node
            for (int i=0; i<subDom_max.n_rows; i++) {
                if (subDom_max(i) == sppPool.topo.subdomain) {
                    bMax_index = index_max(sppPool.bMat_p.row(sppPool.S_p+sppPool.I_p+i));
                    sppPool.bMat_p.row(sppPool.S_p+sppPool.I_p+i).zeros();
                    sppPool.bMat_p(sppPool.S_p+sppPool.I_p+i, bMax_index(0)) = inv;
                } else {
                    sppPool.bMat_p.row(sppPool.S_p+sppPool.I_p+i).zeros();
                }
            }

            if (sppPool.topo.bisec > 0) {
                if (sppPool.uMat_p.n_rows == 0) {
                    sppPool.uMat_p.set_size(sppPool.bMat_p.n_rows, sppPool.bMat_p.n_cols);
                    sppPool.uMat_p.zeros();
                } else {
                    sppPool.uMat_p.insert_rows(sppPool.uMat_p.n_rows, sppPool.bMat_p.n_rows-sppPool.uMat_p.n_rows);
                }
            }
        }
        sppPool.I_c = 0; // reset invader counter
        sppPool.S_c = sppPool.bMat_p.n_rows - sppPool.rMat.n_rows;
    }
}

void Metacommunity::envFluct() {

    // summary:
        // update aboitic turnover and simulates a single CVode timestep

    // required members:
        // sppPool - object of class Species, where model state is stored

    // external function calls:
        // metaCDynamics()
        // sppPool.ouProcess()

    // output:
        // updates to sppPool matrices

    CommunityDynamics dynamics;  // ODE dynamical object
    dynamics.bMat_p = &sppPool.bMat_p;
    dynamics.rMat = &sppPool.rMat;

    if (sppPool.cMat.n_rows != 0) {
        dynamics.cMat = &sppPool.cMat;
    }

    dynamics.dMat = &sppPool.dMat_n; // if dynamics.dMat is initialized, single domain relaxation selected
    dynamics.rho = &sppPool.rho;

    sppPool.ouProcess(); // update the spatially resolved OU process to be added to the growth rate matrix
    dynamics.efMat = &sppPool.efMat; // in CommunityDynamics rMat += efMat
    ODE_state state(& dynamics);
    state.integrate_until(1); // simulate single timestep
    mat tmp(&state[0], 1, dynamics.number_of_variables());
    sppPool.trajectories = join_vert(sppPool.trajectories, tmp); // store trajectory
    mat RE = sppPool.rMat + sppPool.efMat;
    RE.reshape(1,RE.size());
    sppPool.fluctuations = join_vert(sppPool.fluctuations, RE); // store current distrition in R (= rMat + efMat)
}

void Metacommunity::warming(double dTdt, int res, int time) {

    // summary:
        // updates temperature gradient and rMat and simulates a single CVode timestep

    // required members:
        // sppPool - object of class Species, where model state is stored

    // external function calls:
        // metaCDynamics()
        // sppPool.updateRVecTemp()

    // output:
        // updates to sppPool matrices

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

    dynamics.bMat_p = &sppPool.bMat_p;

    if (sppPool.cMat.n_rows != 0) {
        dynamics.cMat = &sppPool.cMat;
    }
    dynamics.dMat = &sppPool.dMat_n; // if dynamics.dMat is initialized, single domain relaxation selected
    dynamics.rho = &sppPool.rho;

    for (int t=0; t<res; t++) { // update every time step, save every res timesteps
        sppPool.topo.T_int += dTdt; // apply warming - T_int updated
        sppPool.updateRVecTemp(); // apply warming - rMat updated
        printf("\rt = %d, T_int = %f", t, sppPool.topo.T_int);
        fflush(stdout);
        dynamics.rMat = &sppPool.rMat; // update rMat seen by CommmunityDynamics
        ODE_state state(& dynamics);
        state.integrate_until(1); // simulate single timestep
    }

    string Bfile = Bpath + "/bMat_w" + to_string(time) + ".mat";
    string Rfile = Bpath + "/rMat_w" + to_string(time) + ".mat";
    sppPool.bMat_p.save(Bfile, raw_ascii);
    sppPool.rMat.save(Rfile, raw_ascii);

    if (time == 0) {
        // store parameter file in same directory for convenience
        string filenameP = Bpath + "/pars.mat";
        boost::filesystem::ofstream params;
        params.open(filenameP);
        params << "date " << date << endl;
        params << "jobID " << jobID << endl;
        params << "experiment " << experiment << endl;
        params <<  "invMax " << invMax << endl;
        params << "parOut " << parOut << endl;
        params << "rep " << rep << endl;
        params << "simTime " << simTime << endl;
        params << "no_nodes " << sppPool.topo.no_nodes << endl;
        params << "bisec " << sppPool.topo.bisec << endl;
        params << "phi " << sppPool.topo.phi << endl;
        params << "envVar " << sppPool.topo.envVar << endl;
        params << "c1 " << sppPool.c1 << endl;
        params << "c2 " << sppPool.c2 << endl;
        if (sppPool.c3 > 0) {
            params << "c3 " << sppPool.c3 << endl;
        }
        params << "emRate " << sppPool.emRate << endl;
        params << "dispL " << sppPool.dispL << endl;
        params << "invasion " << sppPool.invasion << endl;
        params << "iterS " << iterS << endl;
        params << "deltaT " << deltaT << endl;
        params << "tMax " << tMax << endl;
        params << "pProducer " << sppPool.pProducer << endl;
        params << "prodComp " << sppPool.prodComp << endl;
        params << "var_e " << sppPool.topo.var_e << endl;
        params << "alpha " << sppPool.alpha << endl;
        params << "sigma " << sppPool.sigma << endl;
        params << "rho " << sppPool.rho << endl;
        params << "comp_dist " << sppPool.comp_dist << endl;
        params << "T_int " << sppPool.topo.T_int << endl;
        params << "dTdt " << dTdt << endl;
        params.close();
    }
}

void Metacommunity::longDistDisp(int tMax, int edges) {

    // summary: ...

    bool stepwise = true; // stepwise addition of new edges

    if (edges < 0) {
        stepwise = false; // select single addition followed by high res storage of trajectories
    }

    // record row minimum of dispersal operator, if non-zero nodes have degree N-1
    colvec row_min = min(abs(sppPool.dMat_n),1);

    int non_zero_d;
    { // scoped for deletion of sparse cast D
        sp_mat D_sp(sppPool.dMat_n);
        non_zero_d = D_sp.n_nonzero;
    }
    int zero_d = (sppPool.dMat_n.n_elem - non_zero_d)/2;
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
        uvec edge = find(sppPool.dMat_n.row(node(0))==0);
        edge = shuffle(edge);

        // add edges to adjacency matrix (-1 to distinguish between original edges)
        sppPool.topo.adjMat(node(0), edge(0)) = -1.0; // new edges in adjacency matrix allocated -1
        sppPool.topo.adjMat(edge(0), node(0)) = -1.0;

        // regenerate dispersal matrix - new normalization due to altered node degree, new edges not distance weighted
        sppPool.genDispMat();

        // record edge allocation
        pert_record(t_step, 0) = node(0);
        pert_record(t_step, 1) = edge(0);
        pert_record(t_step, 2) = t_interval;
        t_step++;

        // update row min vector
        row_min = min(abs(sppPool.dMat_n),1);

        edge_count++;

        cout << "edge count " << edge_count << endl;

        if (edge_count == edges) {
            if (stepwise) {
                printf("%d new edges allocated\n", edges * (t_interval + 1));
                metaCDynamics(tMax);
                string Bfile = Bpath + "/bMat_d" + to_string(t_interval) + ".mat";
                sppPool.bMat_p.save(Bfile, raw_ascii);
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

                int res = 10;

                // create CommunityDynamics instance for high resolution output
                CommunityDynamics dynamics;  // ODE dynamical object

                dynamics.bMat_p = &sppPool.bMat_p;

                if (sppPool.cMat.n_rows != 0) {
                    dynamics.cMat = &sppPool.cMat;
                }

                dynamics.rMat = &sppPool.rMat; // update rMat seen by CommmunityDynamics
                dynamics.dMat = &sppPool.dMat_n; // if dynamics.dMat is initialized, single domain relaxation selected
                dynamics.rho = &sppPool.rho;

                for (int t=0; t<tMax; t++) { // update every time step, save every res timesteps
                    printf("\rt = %d", t);
                    fflush(stdout);
                    ODE_state state(& dynamics);
                    state.integrate_until(t); // simulate single timestep
                    if ((t > 0) && (t % res == 0)) {
                        string Bfile = Bpath + "/bMat_d" + to_string(t) + ".mat";
                        sppPool.bMat_p.save(Bfile, raw_ascii);
                    }
                }
                row_min.ones(); // break while loop
            }
        }
    }

    if (edge_count) { // final relaxation unless no_nodes % edges = 0
        metaCDynamics(tMax);
        string Bfile = Bpath + "/bMat_d" + to_string(t_interval) + ".mat";
        sppPool.bMat_p.save(Bfile, raw_ascii);
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
    // sppPool - object of class Species, where model state is stored

    // external function calls:
    // metaCDynamics()

    // output:
    // stores biomass matrices

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// Generate output file names ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    storeTraj=3;
    if ((sppPool.topo.scVec.n_rows > 0) & (sppPool.topo.scVec_prime.n_rows > 0)) {
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
    node_id = linspace<uvec>(no_removals,sppPool.topo.no_nodes-1,sppPool.topo.no_nodes-no_removals);

    bool write_to_file = true;

    if (write_to_file) {

        Species sppPool_copy = sppPool;
        for (int x = 0; x < node_id.n_rows; x++) {
            cout << "\nNode " << x << endl;
            sppPool = sppPool_copy;
            sppPool.bMat_p.shed_col(node_id(x));
            if (sppPool.rMat.n_rows > 0) {
                sppPool.rMat.shed_col(node_id(x));
            }
            sppPool.topo.network.shed_row(node_id(x));
            sppPool.topo.distMat.shed_row(node_id(x));
            sppPool.topo.distMat.shed_col(node_id(x));
            sppPool.dMat_n.shed_row(node_id(x)); // corresponding edges also removed
            sppPool.dMat_n.shed_col(node_id(x)); // corresponding edges also removed
            if (sppPool.topo.scVec.n_rows > 0) {
                sppPool.topo.scVec.shed_row(node_id(x));
            }
            if (sppPool.topo.scVec_prime.n_rows > 0) {
                sppPool.topo.scVec_prime.shed_row(node_id(x));
            }
            if (sppPool.topo.envMat.n_cols > 0) {
                sppPool.topo.envMat.shed_col(node_id(x));
            }
            sppPool.topo.no_nodes--;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Relax metacommunty and save new state ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            CommunityDynamics dynamics;  // ODE dynamical object

            dynamics.bMat_p = &sppPool.bMat_p;
            dynamics.rMat = &sppPool.rMat;
            if (sppPool.cMat.n_rows != 0) {
                dynamics.cMat = &sppPool.cMat;
            }

            dynamics.rho = &sppPool.rho;
            if (sppPool.topo.scVec.n_rows != 0) {
                dynamics.scVec = &sppPool.topo.scVec;
                dynamics.scVec_prime = &sppPool.topo.scVec_prime;
            }
            dynamics.S_p = sppPool.S_p;
            dynamics.S_c = sppPool.S_c;

            string Bfile = Bpath + "/bMat_nr_init" + to_string(node_id(x)) + ".mat";
            sppPool.bMat_p.save(Bfile, raw_ascii);
            if (sppPool.S_c > 0) {
                string Bfile_c = Bpath + "/bMat_nr_c_init" + to_string(node_id(x)) + ".mat";
                mat B_c = sppPool.bMat_p.rows(S_p,sppPool.bMat_p.n_rows-1);
                B_c.save(Bfile_c, raw_ascii);
            }
            metaCDynamics(tMax);
            Bfile = Bpath + "/bMat_nr" + to_string(node_id(x)) + ".mat";
            sppPool.bMat_p.save(Bfile, raw_ascii);
            if (sppPool.S_c > 0) {
                string Bfile_c = Bpath + "/bMat_nr_c" + to_string(node_id(x)) + ".mat";
                mat B_c = sppPool.bMat_p.rows(S_p,sppPool.bMat_p.n_rows-1);
                B_c.save(Bfile_c, raw_ascii);
            }
        }
    }
}

void Metacommunity::genJacobian() {

    // summary:
        // generate the numerical approximation of the Jacobian matrix for computing regional competitive overlap matrix

    // required members:
        // sppPool - object of class Species, where model state is stored

    // external function calls:
        // ODE.h::numerical_jacobian()

    // output:
        // jacobian - object for storing output of numerical_jacobian()

    CommunityDynamics dynamics;

    dynamics.bMat_p = &sppPool.bMat_p;
    dynamics.rMat = &sppPool.rMat;
    if (sppPool.cMat.n_rows != 0) {
        dynamics.cMat = &sppPool.cMat;
    }
    dynamics.dMat = &sppPool.dMat_n; // if dynamics.dMat is initialized, single domain relaxation selected
    dynamics.rho = &sppPool.rho;

    {
        ODE_state state(&dynamics);
        for (double t = 0; t < 0; t++) { // accesses back-end CVode initialization machinery
            state.integrate_until(t);
        }
        jacobian.set_size(sppPool.topo.no_nodes * sppPool.bMat_p.n_rows,
                          sppPool.topo.no_nodes * sppPool.bMat_p.n_rows);
        dynamics.numerical_Jacobian(jacobian);
    }
}


void Metacommunity::genCMatReg(double h) {

    // summary:
        // generate a numerically approximated regional interaction matrix via a computation harvesting experiment
        // for detail of algorithm see O'Sullivan et al. (2019)

    // arguments:
        // h - harvesting rate

    // required members:
        // jacobian
        // sppPool - object of class Species, where model state is stored

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

    int S_tot = sppPool.bMat_p.n_rows;
    uvec index = linspace<uvec>(0, sppPool.topo.no_nodes-1, sppPool.topo.no_nodes);
    index *= S_tot; // indexes all populations of focal species i=1; indices for species j generated by element-wise addition
    cMat_reg.set_size(S_tot,S_tot); //
    vec H_i; // storage for harvesting vector
    vec dB_jx; // storage for perturbed, vectorized state matrix

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Iteratively harvest each species i and compute dB_jx ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int i=0; i<S_tot; i++) {
        H_i.zeros(S_tot*sppPool.topo.no_nodes);
        H_i.elem(index + i) = h * sppPool.bMat_p.row(i);
        // "local shift in biomasses of species j due to harvesting of species i per unit h"
        dB_jx = -1 * (J_inv * H_i) / h; // -1 due to definition of competitive interaction coefficent as (+)c_ij
        mat dB_jx_Mat(&dB_jx(0), S_tot, sppPool.topo.no_nodes); // state vector in matrix form (for row sum operation)
        cMat_reg.col(i) = sum(dB_jx_Mat, 1); // store dB_j = sum_x(dB_jx)
    }


    jacobian.reset();
    cMat_reg = inv(cMat_reg); // invert
//    mat norm, sgn;
//    norm.zeros(cMat_reg.n_rows, cMat_reg.n_cols);
//    sgn.zeros(cMat_reg.n_rows, cMat_reg.n_cols);
//    norm.diag() = 1/sqrt(abs(cMat_reg.diag())); // normalization
//    sgn.diag() = cMat_reg.diag() / abs(cMat_reg.diag());
//    cMat_reg = sgn * norm * cMat_reg * norm; // normalize regional competitive overlap matrix
}

void Metacommunity::genSourceSink(int tFullRelax) {

    // summary:
        // infers which populations are dependent upon immigration for local detectability by switching off dispersal and relaxing to equilibrium

    // arguments:
        // tFullRelax - (long) relaxation time

    // required members:
        // sppPool - object of class Species, where model state is stored

    // external function calls:
        // metaCDynamics()

    // output:
        // matrices of dimensions SxN with 1 indicating source, -1 sink populations

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Set off-diagonal elements of dMat_n to zero and relax ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    mat DStore = sppPool.dMat_n;
    mat BStore_p = sppPool.bMat_p;
    mat Src_p, Snk_p;
    sppPool.dMat_n.zeros();
    sppPool.dMat_n.diag() = DStore.diag();
    metaCDynamics(tFullRelax); // relax

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Generate discrete source (1) sink (-1) matrices //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Src_p.zeros(sppPool.bMat_p.n_rows, sppPool.bMat_p.n_cols);
    Src_p(find(sppPool.bMat_p > sppPool.thresh)).ones(); // 1 allocated to detected populations after dispersal switched off

    sppPool.bMat_p = BStore_p; // reset metacommunity objects
    sppPool.dMat_n = DStore;

    Snk_p.zeros(sppPool.bMat_p.n_rows,sppPool.bMat_p.n_cols);

    Snk_p(find(sppPool.bMat_p > sppPool.thresh)).ones(); // 1 allocated to detected populations before dispersal switched off
    Snk_p = Src_p - Snk_p; // -1 allocated to sink populations, 1 allocated to pops excluded by immigration
    Snk_p(find(Snk_p == 1)).zeros(); // weak source populations removed

    sppPool.bMat_p_src = Src_p + Snk_p; // final matrix source (1) sink (-1)
}

void Metacommunity::performanceTest(bool init, int repPer, int S_p, int S_c) {
    // Fill elements of matrices B and R from lognormal and normal distributions for assembly free performance testing

    if (sppPool.pProducer == 1) {
        S_c = 0;
    }

    if (init) {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Randomly populate matrix objects without L-V pruning ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        cout << "Generating random community of " << S_p + S_c << " species, rep " << repPer << endl;

        for (int i = 0; i < S_p; i++) {
            sppPool.invade(0, false); // invade producers - gerate complete cMat and rMat
        }

        if (sppPool.pProducer < 1.0) {
            for (int i = 0; i < S_c; i++) {
                sppPool.invade(1, false); // invade consumer
            }
        }

        sppPool.bMat_p.randn(); // populate bMat_p from normalized log normal distribution
        sppPool.bMat_p = exp(sppPool.bMat_p);
        sppPool.bMat_p = sppPool.bMat_p / sppPool.bMat_p.max();

    } else {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Custom output command ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        string filenameBPer, filenameBCPer, filenamePPer, filenameDPer, filenameNPer, filenameNOPer;
        ostringstream nameBPer, nameBCPer, namePPer, nameDPer, nameNPer, nameNOPer;
        string directory;

        ostringstream p;
        p << "SimulationData/Performance_experiment/N=" << sppPool.topo.no_nodes << "/" << date << "/";
        outputDirectory += p.str();

        if (!exists(outputDirectory)) { // make directory if doesn't currently exist
            create_directories(outputDirectory);
        }

        nameBPer << outputDirectory << sppPool.topo.bisec << "_bMatP" << repPer << "_" << sppPool.bMat_p.n_rows
                 << ".mat";
        nameBCPer << outputDirectory << sppPool.topo.bisec << "_bMat_cP" << repPer << "_"
                  << sppPool.bMat_p.n_rows << ".mat";
        namePPer << outputDirectory << sppPool.topo.bisec << "_paramsP" << repPer << "_" << sppPool.bMat_p.n_rows
                 << ".mat";
        filenameBPer = nameBPer.str();
        filenameBCPer = nameBCPer.str();
        filenamePPer = namePPer.str();

        cout << "Saving to " << filenameBPer << endl;

        if (repPer == 0) { // for first repeat of performance test only, store biomass matrices
            if (sppPool.rMat.n_rows != 0) {
                mat B_p = sppPool.bMat_p.rows(0,sppPool.rMat.n_rows-1);
                B_p.save(filenameBPer, raw_ascii);
            }
            if (sppPool.bMat_p.n_rows > sppPool.rMat.n_rows) {
                mat B_c = sppPool.bMat_p.rows(sppPool.rMat.n_rows, sppPool.bMat_p.n_rows-1);
                B_c.save(filenameBCPer, raw_ascii);
            }
        }

        // generate parameter file
        boost::filesystem::ofstream params;
        params.open(filenamePPer);
        params << "date " << date << endl;
        params << "experiment " << experiment << endl;
        params << "simTime " << simTime << endl;
        params << "no_nodes " << sppPool.topo.no_nodes << endl;
        params << "bisec " << sppPool.topo.bisec << endl;
        params << "phi " << sppPool.topo.phi << endl;
        params << "envVar " << sppPool.topo.envVar << endl;
        params << "c1 " << sppPool.c1 << endl;
        params << "c2 " << sppPool.c2 << endl;
        params << "emRate " << sppPool.emRate << endl;
        params << "dispL " << sppPool.dispL << endl;
        params << "invasion " << sppPool.invasion << endl;
        params << "iterS " << iterS << endl;
        params << "deltaT " << deltaT << endl;
        params << "tMax " << tMax << endl;
        params << "pProducer " << sppPool.pProducer << endl;
        params << "prodComp " << sppPool.prodComp << endl;
        params << "var_e " << sppPool.topo.var_e << endl;
        params << "alpha " << sppPool.alpha << endl;
        params << "sigma " << sppPool.sigma << endl;
        params << "rho " << sppPool.rho << endl;
        params << "comp_dist " << sppPool.comp_dist << endl;
        params.close();
    }
}

void Metacommunity::writePars(int repPer) {
    string filenameBPer, filenameBCPer, filenamePPer, filenameDPer, filenameNPer, filenameNOPer;
    ostringstream nameBPer, nameBCPer, namePPer, nameDPer, nameNPer, nameNOPer;
    string directory;

    ostringstream p;
    p << "SimulationData/Assembly_experiment/N=" << sppPool.topo.no_nodes << "/" << date << "/";
    outputDirectory += p.str();

    if (!exists(outputDirectory)) { // make directory if doesn't currently exist
        create_directories(outputDirectory);
    }

    namePPer << outputDirectory << sppPool.topo.bisec << "_paramsP" << repPer << "_" << sppPool.bMat_p.n_rows
             << ".mat";
    filenamePPer = namePPer.str();

    cout << "Saving to " << filenamePPer << endl;

    // generate parameter file
    boost::filesystem::ofstream params;
    params.open(filenamePPer);
    params << "date " << date << endl;
    params << "experiment " << experiment << endl;
    params << "simTime " << simTime << endl;
    params << "no_nodes " << sppPool.topo.no_nodes << endl;
    params << "bisec " << sppPool.topo.bisec << endl;
    params << "phi " << sppPool.topo.phi << endl;
    params << "envVar " << sppPool.topo.envVar << endl;
    params << "c1 " << sppPool.c1 << endl;
    params << "c2 " << sppPool.c2 << endl;
    params << "emRate " << sppPool.emRate << endl;
    params << "dispL " << sppPool.dispL << endl;
    params << "invasion " << sppPool.invasion << endl;
    params << "iterS " << iterS << endl;
    params << "deltaT " << deltaT << endl;
    params << "tMax " << tMax << endl;
    params << "pProducer " << sppPool.pProducer << endl;
    params << "prodComp " << sppPool.prodComp << endl;
    params << "var_e " << sppPool.topo.var_e << endl;
    params << "alpha " << sppPool.alpha << endl;
    params << "sigma " << sppPool.sigma << endl;
    params << "rho " << sppPool.rho << endl;
    params << "comp_dist " << sppPool.comp_dist << endl;
    params.close();
}

void Metacommunity::printParams() {

    // summary:
        // print model parameterization to console

    printf("\nModel parameters:\n");
    printf("\niterS %d", iterS);
    printf("\ndeltaT %d", deltaT);
    printf("\ntMax %d", tMax);
    printf("\nparOut %f", parOut);
    cout << "\nexperiment " << experiment;
    printf("\nrep %d", rep);
    printf("\nsimTime %f", simTime);
    cout << "\njobID " << jobID;
    cout << "\ndate " << date;
    printf("\n\nc1 %f", sppPool.c1);
    printf("\nc2 %f", sppPool.c2 );
    printf("\nemRate %f", sppPool.emRate);
    printf("\ndispL %f", sppPool.dispL);
    printf("\npProducer %f", sppPool.pProducer) ;
    printf("\nalpha %e", sppPool.alpha);
    printf("\nsigma %f", sppPool.sigma);
    printf("\nrho %f", sppPool.rho);
    printf("\ncomp_dist %d", sppPool.comp_dist);
    printf("\nomega %f", sppPool.omega);
    printf("\n\nno_nodes %d", sppPool.topo.no_nodes);
    printf("\nphi %f", sppPool.topo.phi);
    printf("\nenvVar %d", sppPool.topo.envVar);
    printf("\nT_int %f", sppPool.topo.T_int);
    printf("\nvar_e %f", sppPool.topo.var_e);
    printf("\nrandGraph %d", sppPool.topo.randGraph);
    printf("\ngabriel %d", sppPool.topo.gabriel);
    printf("\nbisec %d", sppPool.topo.bisec);
    printf("\n\ninvMax %d", invMax);
}

void Metacommunity::snapShot() {
    // summary:
        // store biomass and growth rate matrices for continuous analysis of metacommunity during assembly

    string Bpath;

    if (bMatFileName.length() == 0) { // generate path from scratch

	    string filenameBF;
        ostringstream nameBF;

        ostringstream p, q;
        p << "SimulationData/N=" << sppPool.topo.no_nodes << "/" << experiment
          << "_experiment/" << date << "/";
        outputDirectory += p.str();

        nameBF << outputDirectory << date << "_" << experiment << "("
               << invMax << ")" << parOut << "bMat" << rep << ".mat";
        filenameBF = nameBF.str();
        bMatFileName = filenameBF; // set bMatFileName for follow up snapshots

        q << parOut << "bMat" << rep << "/assembly/";
        outputDirectory += q.str();

	    Bpath = outputDirectory;

    } else { // use find and replace to generate paths from bMatFileName

        std::size_t pos1 = bMatFileName.find_last_of("/");
        string bMatDir = bMatFileName.substr(0,pos1);
        pos1 = bMatFileName.find_last_of(")");
        std::size_t pos2 = bMatFileName.find_last_of(".");
        string newFolder = bMatFileName.substr(pos1+1, pos2-pos1-1);
        Bpath = bMatDir + "/" + newFolder + "/assembly/";

    }

    if (!exists(Bpath)) { // make directory if doesn't currently exist
	create_directories(Bpath);
    }

    double averaging_time = 0; // this should be a command line option
    if(averaging_time > 0){ // hack to output temporal averages

	if(sppPool.topo.network.n_rows < sppPool.topo.no_nodes){
	    cerr << "\nAveraging in parallel mode not implemented";
	    exit(1);
	}
	cout << "\nComputing averages... ";
	cout.flush();
	bool old_compute_averages = compute_averages;
	compute_averages = true;
	metaCDynamics(averaging_time);
	compute_averages = old_compute_averages;
	
	string filenameBav = Bpath + "bavMat_" + to_string(sppPool.invEvent) + ".mat";
	if (sppPool.bavMat_p.n_rows > 0) {
	    sppPool.bavMat_p.save(filenameBav, raw_ascii);
	}
	cout << "done.";
	cout.flush();
    }

    string filenameB = Bpath + "bMat_" + to_string(sppPool.invEvent) + ".mat";
    string filenameR = Bpath + "rMat_" + to_string(sppPool.invEvent) + ".mat";
    string filenameBC = Bpath + "bMat_c_" + to_string(sppPool.invEvent) + ".mat";

    if (sppPool.rMat.n_rows != 0) {
        mat B_p = sppPool.bMat_p.rows(0,sppPool.rMat.n_rows-1);
        B_p.save(filenameB, raw_ascii);
    }
    if (sppPool.rMat.n_rows > 0) {
        sppPool.rMat.save(filenameR, raw_ascii);
    }
    if (sppPool.bMat_p.n_rows > sppPool.rMat.n_rows) {
        mat B_c = sppPool.bMat_p.rows(sppPool.rMat.n_rows, sppPool.bMat_p.n_rows-1);
        B_c.save(filenameBC, raw_ascii);
    }
}

void Metacommunity::cleanup() {

    // summary:
        // Reset model objects to pre-invasion state in case signal received during invader tesing

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Remove untested producer species for signal handler output //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    mat BStore_p, BStore_c, RStore, UStore_p, UStore_c, AStore_r, AStore_c, TStore, CStore_r, CStore_c;

    if (sppPool.cMat.n_rows > sppPool.S_p) {
        CStore_r = sppPool.cMat.rows(sppPool.S_p, sppPool.cMat.n_rows - 1);
        sppPool.cMat.resize(sppPool.S_p, sppPool.bMat_p.n_rows);
    }
    if (sppPool.cMat.n_cols > sppPool.S_p) {
        CStore_c = sppPool.cMat.cols(sppPool.S_p, sppPool.cMat.n_cols-1);
        sppPool.cMat.resize(sppPool.S_p, sppPool.S_p);
    }
    if (sppPool.bMat_p.n_rows > sppPool.S_p) {
        BStore_p = sppPool.bMat_p.rows(sppPool.S_p, sppPool.bMat_p.n_rows - 1);
        sppPool.bMat_p.resize(sppPool.S_p, sppPool.topo.no_nodes);
    }
    if (sppPool.rMat.n_rows > sppPool.S_p) {
        RStore = sppPool.rMat.rows(sppPool.S_p, sppPool.rMat.n_rows - 1);
        sppPool.rMat.resize(sppPool.S_p, sppPool.topo.no_nodes);
    }
    if (sppPool.tMat.n_rows > sppPool.S_p) {
        TStore = sppPool.tMat.rows(sppPool.S_p, sppPool.tMat.n_rows-1);
        sppPool.tMat.resize(sppPool.S_p, sppPool.topo.envVar);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Output data ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    outputData();

}

void Metacommunity::outputData() {

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
        p << "SimulationData/N=" << sppPool.topo.no_nodes << "/" << experiment
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
              << invMax << ")" << parOut << "dMat_n" << rep << ".mat";
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
        nameF << outputDirectory << date << "_" << experiment << "("
              << invMax << ")" << parOut << "fVec" << rep << ".mat";
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
        filenameD.replace(b, string("bMat").length(), "dMat_n");
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
        filenameF = bMatFileName;
        filenameF.replace(b, string("bMat").length(), "fVec");
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
/////////////////////////////////// Create/edit file recording current state ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (jobID != "NA") {
        cout << "\nOutputting current state... ";
        ostringstream st;
        st << jobID << ".txt"; // name of file in which state recorded
        string state = st.str();
        std::ifstream stateFile;
        stateFile.open(state);

        if (!stateFile.is_open()) {
            cout << "generating state file " << state << endl;
            std::ofstream stateFileNew;
            stateFileNew.open(state);
            stateFileNew << sppPool.invasion << endl;
            stateFileNew << invMax << endl;
            stateFileNew << bMatFileName << endl;
            stateFileNew << "1" << endl;
            stateFileNew.close();
        } else {
            cout << "updating state file " << state << endl;
            std::vector<string> state_vec;
            string state_line;
            while ( stateFile.good() ) {
                getline (stateFile,state_line);
                state_vec.push_back(state_line);
            }
            stateFile.close();
            state_vec[0] = to_string(sppPool.invasion);
            std::ofstream stateFileNew;
            stateFileNew.open(state);
            if (stateFileNew.is_open()) {
                for (int i =0; i<state_vec.size()-1; i++) {
                    stateFileNew << state_vec[i] << endl;
                }
                stateFileNew.close();
            }
            stateFileNew.close();
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Store matrix objects ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Note, if specific experiment run, only relevant objects are output
    if (sppPool.trajectories.n_rows != 0) { // output trajectories ONLY
        sppPool.trajectories.save(filenameTr, raw_ascii);
        if (sppPool.efMat.n_rows != 0) {
            mat R = sppPool.rMat;
            R.reshape(1, R.size());
            rowvec N(R.n_cols);
            N.fill(datum::nan);
            R = join_vert(R, N);
            sppPool.fluctuations = join_vert(R, sppPool.fluctuations);
            sppPool.fluctuations.save(filenameEf, raw_ascii);
        }
    } else if (cMat_reg.n_rows != 0) { // output regional competitive overlap matrix ONLY
        cMat_reg.save(filenameC, raw_ascii);
    } else if (sppPool.bMat_p_src.n_rows != 0) { // output source-sink matrices ONLY
        mat B_p_src = sppPool.bMat_p_src.rows(0,sppPool.rMat.n_rows-1);
        B_p_src.save(filenameBs, raw_ascii);
        if (sppPool.bMat_p.n_rows > sppPool.rMat.n_rows) {
            mat B_c_src = sppPool.bMat_p_src.rows(sppPool.rMat.n_rows, sppPool.bMat_p.n_rows-1);
            B_c_src.save(filenameBCs, raw_ascii);
        }
    } else { // output standard model objects
        sppPool.topo.network.save(filenameN, raw_ascii);
        sppPool.rMat.save(filenameR, raw_ascii);
        if (sppPool.sMat.n_rows != 0) {
            sppPool.sMat.save(filenameSR, raw_ascii);
        }
        if (bMatFileName.length() == 0) {
            mat B_p = sppPool.bMat_p.rows(0,sppPool.S_p-1);
            B_p.save(filenameB, raw_ascii);
        } else {
            mat B_p = sppPool.bMat_p.rows(0,sppPool.S_p-1);
            B_p.save(bMatFileName, raw_ascii);
        }
        if (sppPool.emMat_p.n_rows != 0) {
            mat Em_p = sppPool.emMat_p.rows(0,sppPool.rMat.n_rows-1);
            Em_p.save(filenameEM, raw_ascii);
            if (sppPool.bMat_p.n_rows > sppPool.rMat.n_rows) {
                mat Em_c = sppPool.emMat_p.rows(sppPool.rMat.n_rows, sppPool.bMat_p.n_rows-1);
                Em_c.save(filenameEMC, raw_ascii);
            }
        }

        if (sppPool.bMat_p.n_rows > sppPool.rMat.n_rows) {
            mat B_c = sppPool.bMat_p.rows(sppPool.rMat.n_rows, sppPool.bMat_p.n_rows-1);
            B_c.save(filenameBC, raw_ascii);
        }
        if (sppPool.topo.envMat.n_rows != 0) {
            sppPool.topo.envMat.save(filenameE, raw_ascii);
            sppPool.tMat.save(filenameT, raw_ascii);
        }
        sppPool.sppRichness.save(filenameS, raw_ascii);
        sppPool.dMat_n.save(filenameD, raw_ascii);
        if (sppPool.S_p > 0) {
            mat C = sppPool.cMat.submat(0,0,sppPool.S_p-1,sppPool.S_p-1);
            C.save(filenameI, raw_ascii);
        }
        if (sppPool.bMat_p.n_rows > sppPool.rMat.n_rows) {
            mat A = sppPool.cMat.submat(sppPool.S_p,0,sppPool.cMat.n_rows-1,sppPool.S_p-1);
            A.save(filenameA, raw_ascii);
        }

        boost::filesystem::ofstream params;
        params.open(filenameP);
        params << "date " << date << endl;
        params << "jobID " << jobID << endl;
        params << "experiment " << experiment << endl;
        params <<  "invMax " << invMax << endl;
        params << "parOut " << parOut << endl;
        params << "rep " << rep << endl;
        params << "simTime " << simTime << endl;
        params << "no_nodes " << sppPool.topo.no_nodes << endl;
        params << "bisec " << sppPool.topo.bisec << endl;
        params << "randGraph " << sppPool.topo.randGraph << endl;
        params << "phi " << sppPool.topo.phi << endl;
        params << "envVar " << sppPool.topo.envVar << endl;
        params << "c1 " << sppPool.c1 << endl;
        params << "c2 " << sppPool.c2 << endl;
        if (sppPool.c3 > 0) {
            params << "c3 " << sppPool.c3 << endl;
        }
        params << "emRate " << sppPool.emRate << endl;
        params << "dispL " << sppPool.dispL << endl;
        params << "invasion " << sppPool.invasion << endl;
        params << "iterS " << iterS << endl;
        params << "deltaT " << deltaT << endl;
        params << "tMax " << tMax << endl;
        params << "pProducer " << sppPool.pProducer << endl;
        params << "prodComp " << sppPool.prodComp << endl;
        params << "var_e " << sppPool.topo.var_e << endl;
        params << "alpha " << sppPool.alpha << endl;
        params << "sigma " << sppPool.sigma << endl;
        params << "rho " << sppPool.rho << endl;
        params << "comp_dist " << sppPool.comp_dist << endl;
        params << "T_int " << sppPool.topo.T_int << endl;
        params << "omega " << sppPool.omega << endl;
        if (sppPool.topo.skVec.n_rows > 0) {
            params << "sk " << sppPool.topo.skVec(0) << endl;
        } else {
            params << "sk " << "NULL" << endl;
        }
        params << "delta_g " << sppPool.delta_g << endl;
        params << "sigma_r " << sppPool.sigma_r << endl;
        params << "symComp " << sppPool.symComp << endl;
        params.close();
    }
}

bool fexists(const std::string& filename) { // check if file/directory exists
    std::ifstream ifile(filename.c_str());
    return (bool)ifile;
}

void Metacommunity::importData(string bMatFile) {

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
    filenameD.replace(b, string("bMat").length(), "dMat_n");
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
    filenameF = bMatFileName;
    filenameF.replace(b, string("bMat").length(), "fVec");

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
        } else if (parName[i] == "jobID") {
            jobID = pars[i];
            continue;
        } else if (parName[i] == "experiment") {
            experiment = pars[i];
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
            sppPool.topo.no_nodes = stoi(pars[i]);
            continue;
        } else if (parName[i] == "bisec") {
            sppPool.topo.bisec = stoi(pars[i]);
            continue;
        } else if (parName[i] == "randGraph") {
            sppPool.topo.randGraph = stoi(pars[i]);
            continue;
        } else if (parName[i] == "phi") {
            sppPool.topo.phi = stof(pars[i]);
            continue;
        } else if (parName[i] == "envVar") {
            sppPool.topo.envVar = stoi(pars[i]);
            continue;
        } else if (parName[i] == "c1") {
            sppPool.c1 = stof(pars[i]);
            continue;
        } else if (parName[i] == "c2") {
            sppPool.c2 = stof(pars[i]);
            continue;
        } else if (parName[i] == "c3") {
            sppPool.c3 = stof(pars[i]);
            continue;
        } else if (parName[i] == "emRate") {
            sppPool.emRate = stof(pars[i]);
            continue;
        } else if (parName[i] == "dispL") {
            sppPool.dispL = stof(pars[i]);
            continue;
        } else if (parName[i] == "invasion") {
            sppPool.invasion = stoi(pars[i]);
            sppPool.invEvent = stoi(pars[i]);
            sppPool.inv0 = sppPool.invasion;
            continue;
        } else if (parName[i] == "iterS") {
            iterS = stoi(pars[i]);
            continue;
        } else if (parName[i] == "deltaT") {
            deltaT = stoi(pars[i]);
            continue;
        } else if (parName[i] == "tMax") {
            tMax = stoi(pars[i]);
            continue;
        } else if (parName[i] == "pProducer") {
            sppPool.pProducer = stof(pars[i]);
            continue;
        } else if (parName[i] == "pProducer") {
            if (stoi(pars[i])) {
                sppPool.prodComp = true;
            } else {
                sppPool.prodComp = false;
            }
            continue;
        } else if (parName[i] == "var_e") {
            sppPool.topo.var_e = stof(pars[i]);
            continue;
        } else if (parName[i] == "alpha") {
            sppPool.alpha = stof(pars[i]);
            continue;
        } else if (parName[i] == "sigma") {
            sppPool.sigma = stof(pars[i]);
            continue;
        } else if (parName[i] == "rho") {
            sppPool.rho = stof(pars[i]);
            continue;
        } else if (parName[i] == "comp_dist") {
            sppPool.comp_dist = stoi(pars[i]);
            continue;
        } else if (parName[i] == "T_int") {
            sppPool.topo.T_int = stof(pars[i]);
            continue;
        } else if (parName[i] == "omega") {
            sppPool.omega = stof(pars[i]);
            continue;
        } else if (parName[i] == "sk") {
            sppPool.topo.skVec.set_size(1);
            if ((pars[i] == "import") || (pars[i] == "NULL")){
                sppPool.topo.skVec(0) = 0.0;
            } else {
                sppPool.topo.skVec(0) = stof(pars[i]);
            }
            continue;
        } else if (parName[i] == "delta_g") {
            sppPool.delta_g = stof(pars[i]);
            continue;
        } else if (parName[i] == "sigma_r") {
            sppPool.sigma_r = stof(pars[i]);
            continue;
        } else if (parName[i] == "symComp") {
            if (stoi(pars[i])) {
                sppPool.symComp = true;
            } else {
                sppPool.symComp = false;
            }
            continue;
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Load matrix objects ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sppPool.bMat_p.load(bMatFileName);
    sppPool.S_p = sppPool.bMat_p.n_rows;

    sppPool.dMat_n.load(filenameD);
    if(fexists(filenameI)) {
        sppPool.cMat.load(filenameI);
    }
    if(fexists(filenameEM)) {
        sppPool.emMat_p.load(filenameEM);
    }
    sppPool.topo.network.load(filenameN);
    if (fexists(filenameT)) {
        sppPool.tMat.load(filenameT);
    }

    if (fexists(filenameE)) {
        sppPool.topo.envMat.load(filenameE);
    }
    if (fexists(filenameF)) {
        sppPool.topo.fVec.load(filenameF);
    }
    sppPool.rMat.load(filenameR);
    if (fexists(filenameSR)) {
        sppPool.sMat.load(filenameSR);
    }
    sppPool.sppRichness.load(filenameS);
    if (fexists(filenameIPp)) {
        invasionProb_p.load(filenameIPp);
    }
    if (fexists(filenameIPc)) {
        invasionProb_c.load(filenameIPc);
    }

    cout << "\nImported network rows 0-4 = " << endl;
    cout << sppPool.topo.network.rows(0,min((int) sppPool.topo.network.n_rows-1, 4));
    sppPool.topo.genDistMat();
    printf("\nSpecies richness, S_p = %d, S_c = %d\n", (int) sppPool.rMat.n_rows, (int) sppPool.bMat_p.n_rows - (int) sppPool.rMat.n_rows);
}

// Local Variables:
// c-file-style: "stroustrup"
// End:
