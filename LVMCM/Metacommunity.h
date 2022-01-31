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

#ifndef LVMCM_METACOMMUNITY_H
#define LVMCM_METACOMMUNITY_H

#include <iostream>
#include <armadillo>

#include "Topography.h"
#include "Species.h"
#include "ODE.h"
#include "error.h"
#include "CommunityDynamics.h"

extern int g_form_of_dynamics; // 0: LVMCM; 1: PSD
extern int g_block_transitions; // 0: S->P/D; P->D transitions permited; 1: blocked

using namespace std;
using namespace arma;

class Metacommunity: public Species {
public:
    // members
    Species spp; // biotic community sppP
    
    // assembly parameters
    int invMax; // total number invasions
    int tMax; // relaxation time
    int S_p = 0; // record species richness for cleanup
    int S_c = 0; // record species richness for cleanup

    // matrix objects - check each of these is required
    mat jacobian; // stores numerical approximation of Jacobian matrix
    mat cMat_reg; // stores effective interaction matrix computed numerically using harvesting experiment

    // data handling objects
    string bMatFileName;
    double parOut;
    string experiment;
    int rep;
    string date;
    double simTime = 0;
    string outputDirectory;

    // switches
    int storeTraj = 0; // 0, turn off; 1, record, don't concatenate; 2, concatenate trajectories for multiple invasions
    bool compute_averages = false; // when true, unconditionally populate bavMat_p and bavMat_c after standard relaxation
    bool store_params = true; // when true params and abiotic distributions saved by assembly video function

    // methods
    void metaCDynamics(int T); // ODE solver - NxS coupled ODEs - simulation metacommunity dynamics
    void invaderSample(int trophLev, int no_invaders); // introduce and test new species
    void envFluct(); // simulate dynamics in context of temporal abiotic turnover
    void warming(double dTdt, int res, int time); // simulate dynamics in context of regional warming
    void longDistDisp(int tMax, int edges); // randomly add fast, dispersal between (potentially) distant nodes
    void consArea(int tMax, int percentLandscape, double phi, int rep, bool write_to_file = false); // model conservation areas using binary random field and simulate dynamics
    void nodeRemoval(int tMax, int no_removals = 0); //
    void genJacobian(); // generates numerical approximation of Jacobian
    void genCMatReg(double h = 0.001); // simulates harvesting experiment
    void genSourceSink(int tFullRelax = 1000); // switch off dispersal to assign source-sink populations
    void writePars(string filenameP); // write parameters to file
    void saveVideo(string path); // store objects during assembly for generating assembly videos

    // book keeping functions - store, output and import data
    void printParams();
    void saveMC();
    void loadMC(string bMat);

    // (default) constructor
    Metacommunity () {}

    // (default) deconstructor
    ~Metacommunity () {}

    // intialization constructor
    Metacommunity (
            // Metacommunity parameters
            bool a_init,
            string a_bMat,
            string a_xMat,
            string a_scMat,
            int a_invMax,
            int a_tMax,
            string a_outputDirectory,
            // Species parameters
            double a_c1,
            double a_c2,
            double a_c3,
            double a_emRate,
            double a_dispL,
            double a_pProducer,
            bool a_prodComp,
            bool a_symComp,
            double a_alpha,
            double a_sigma,
            double a_sigma_t,
            double a_rho,
            int a_comp_dist,
            double a_omega,
            double a_dispNorm,
            // Topograpy parameters
            int a_no_nodes,
            int a_lattice_height,
            int a_lattice_width,
            double a_phi,
            int a_envVar,
            vec a_skVec,
            double a_var_e,
            bool a_randGraph,
            bool a_gabriel,
            double a_T_int,
            string a_envMat,
            // output variables
            double a_parOut,
            string a_experiment,
            int a_rep)
    {
        // Parameterize
        if (a_init) { // initialize new metacommunity model
            invMax = a_invMax;
            tMax = a_tMax;
            parOut = a_parOut;
            experiment = a_experiment;
            rep = a_rep;
            time_t t = time(0);
            struct tm * now = localtime( & t );
            ostringstream dateTemp;
            dateTemp << (now->tm_year + 1900) << '-'
                     << (now->tm_mon + 1) << '-'
                     <<  now->tm_mday;
            date = dateTemp.str();
            outputDirectory = a_outputDirectory;
            spp.c1 = a_c1;
            spp.c2 = a_c2;
            spp.c3 = a_c3;
            spp.emRate = a_emRate;
            spp.dispL = a_dispL;
            spp.pProducer = a_pProducer;
            spp.prodComp = a_prodComp;
            spp.symComp = a_symComp;
            spp.alpha = a_alpha;
            spp.sigma = a_sigma;
            spp.sigma_t = a_sigma_t;
            spp.rho = a_rho;
            spp.comp_dist = a_comp_dist;
            spp.omega = a_omega;
            spp.dispNorm = a_dispNorm;
            spp.topo.xMatFileName = a_xMat;
            spp.topo.scMatFileName = a_scMat;
            spp.topo.no_nodes = a_no_nodes;
            spp.topo.lattice_height = a_lattice_height;
            spp.topo.lattice_width = a_lattice_width;
            spp.topo.phi = a_phi;
            spp.topo.envVar = a_envVar;
            if (a_skVec(0) > 0) {
                spp.topo.skVec = a_skVec;
            }
            spp.topo.var_e = a_var_e;
            spp.topo.randGraph = a_randGraph;
            spp.topo.gabriel = a_gabriel;
            spp.topo.envMatFileName = a_envMat;
            spp.topo.network.reset();
            spp.topo.distMat.reset();
            spp.topo.adjMat.reset();
        } else { // import metacommunity model
            loadMC(a_bMat);
            spp.topo.scMatFileName = a_scMat;
            outputDirectory = a_outputDirectory;
            if (a_invMax != 0) {
                invMax += a_invMax;
                cout << "\ninvMax updated to " << invMax << endl;
            }
            if (a_tMax != tMax) {
                tMax = a_tMax;
                cout << "\ntMax updated to " << tMax << endl;
            }
        }
    }
};

#endif //LVMCM_METACOMMUNITY_H

// Local Variables:
// c-file-style: "stroustrup"
// End:
