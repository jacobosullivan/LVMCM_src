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

#ifndef LVMCM_METACOMMUNITY_H
#define LVMCM_METACOMMUNITY_H

#include <iostream>
#include <armadillo>

#include "Topography.h"
#include "Species.h"
#include "ODE.h"
#include "error.h"
#include "CommunityDynamics.h"

using namespace std;
using namespace arma;

class Metacommunity: public Species {
public:
    // members
    Species sppPool; // biotic community sppP
    
    // assembly parameters
    int invMax; // total number invasions
    int iterS; // number Schwarz iterations
    int deltaT; // size of Schwarz timewindow
    int tMax; // relaxation time
    int S_p = 0; // record species richness for cleanup
    int S_c = 0; // record species richness for cleanup

    // matrix objects - check each of these is required
    mat jacobian; // stores numerical approximation of Jacobian matrix
    mat cMat_reg; // stores effective interaction matrix computed numerically using harvesting experiment
    mat invasionProb_p; // record fraction of successful invaders - producers --- obselete
    mat invasionProb_c; // record fraction of successful invaders - consumers --- obselete
    mat invasionProb;

    // data handling objects
    string bMatFileName;
    string jobID = "NA";
    double parOut;
    string experiment;
    int rep;
    string date;
    double simTime = 0;
    int rank = 0;
    string outputDirectory;

    // switches
    int storeTraj = 0; // 0, turn off; 1, record, don't concatenate; 2, concatenate trajectories for multiple invasions
    bool compute_averages = false; // when true, unconditionally populate bavMat_p and bavMat_c after standard relaxation

    // methods
    void metaCDynamics(int T, bool testing = false); // ODE solver - NxS coupled ODEs - simulation metacommunity dynamics
    uvec invaderSample(int trophLev, int no_invaders); // introduce and test new species
    mat invaderCleanup(int trophLev, uvec posGrowth); // remove unsuccessful invaders
    void invaderPopulate(int trophLev, mat bInv_max); // reset biomasses of invaders
    void envFluct(); // simulate dynamics in context of temporal abiotic turnover
    void warming(double dTdt, int res, int time); // simulate dynamics in context of regional warming
    void longDistDisp(int tMax, int edges); // randomly add fast, dispersal between (potentially) distant nodes
    void consArea(int tMax, int percentLandscape, double phi, int rep, bool write_to_file = false); // model conservation areas using binary random field and simulate dynamics
    void nodeRemoval(int tMax, int no_removals = 0); //
    void genJacobian(); // generates numerical approximation of Jacobian
    void genCMatReg(double h = 0.001); // simulates harvesting experiment
    void genSourceSink(int tFullRelax = 1000); // switch off dispersal to assign source-sink populations
    void performanceTest(bool init, int repPer, int S_p = 0, int S_c = 0); // peformance testing of dynamics
    void writePars(int repPer); // write parameters to file for performance testing

    // book keeping functions - store, output and import data
    void printParams();
    void outputData();
    void snapShot();
    void cleanup();
    void importData(string bMat);

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
            int a_iterS,
            int a_deltaT,
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
            double a_phi,
            int a_envVar,
            vec a_skVec,
            double a_var_e,
            bool a_randGraph,
            bool a_gabriel,
            int a_bisec,
            double a_T_int,
            // output variables
            double a_parOut,
            string a_experiment,
            int a_rep,
            string a_jobID
        )
    {
        // Parameterize
        if (a_init) { // initialize new metacommunity model
            invMax = a_invMax;
            iterS = a_iterS;
            deltaT = a_deltaT;
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
            jobID = a_jobID;
            sppPool.c1 = a_c1;
            sppPool.c2 = a_c2;
            sppPool.c3 = a_c3;
            sppPool.emRate = a_emRate;
            sppPool.dispL = a_dispL;
            sppPool.pProducer = a_pProducer;
            sppPool.prodComp = a_prodComp;
            sppPool.symComp = a_symComp;
            sppPool.alpha = a_alpha;
            sppPool.sigma = a_sigma;
            sppPool.sigma_t = a_sigma_t;
            sppPool.rho = a_rho;
            sppPool.comp_dist = a_comp_dist;
            sppPool.omega = a_omega;
            sppPool.dispNorm = a_dispNorm;
            sppPool.topo.xMatFileName = a_xMat;
            sppPool.topo.scMatFileName = a_scMat;
            sppPool.topo.no_nodes = a_no_nodes;
            sppPool.topo.phi = a_phi;
            sppPool.topo.envVar = a_envVar;
            if (a_skVec(0) > 0) {
                sppPool.topo.skVec = a_skVec;
            }
            sppPool.topo.var_e = a_var_e;
            sppPool.topo.randGraph = a_randGraph;
            sppPool.topo.gabriel = a_gabriel;
            sppPool.topo.bisec = a_bisec;
            sppPool.topo.T_int = a_T_int;
            sppPool.topo.network.reset();
            sppPool.topo.distMat.reset();
            sppPool.topo.adjMat.reset();


        } else { // import metacommunity model
            importData(a_bMat);
            sppPool.topo.scMatFileName = a_scMat;
            outputDirectory = a_outputDirectory;
            if (a_invMax != 0) {
                invMax += a_invMax;
                cout << "\ninvMax updated to " << invMax << endl;
            }
            if (a_bisec != 0) {
                sppPool.topo.bisec = a_bisec;
                cout << "\nbisec updated to " << sppPool.topo.bisec << endl;
            }
            if (a_iterS != 0) {
                iterS = a_iterS;
                cout << "\niterS updated to " << iterS << endl;
            }
            if (a_deltaT != 0) {
                deltaT = a_deltaT;
                cout << "\ndeltaT updated to " << deltaT << endl;
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
