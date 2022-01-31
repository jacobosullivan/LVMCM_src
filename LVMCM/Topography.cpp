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
 * This class contains the members and methods required for generating and storing the abiotic component of the model
 * including the (domain decomposed) spatial network, and explicitly modelled environmental distribiutions
 */

#include "Topography.h"
#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

void Topography::genNetwork() {

    // summary:
    // randomly sample patch coordinates from uniform distribution, store in network
    // coordinates are ranked in order of norm (vector distance from origin)

    // required members:
    // no_nodes - number of nodes in spatial network
    // randGraph - select random planar graph (T) or regular lattice (F)

    // output:
    // network - cartesian coordinates of spatial network

    // external function calls:
    // genDistMat()
    // genAdjMat()

    if (xMatFileName.length() == 0) { // generate network

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// Sample coordinates of patches /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (network.n_rows == 0) {
            if (randGraph) { // generate random graph
                vec xcoord = sqrt(no_nodes) * randu(no_nodes);
                vec ycoord = sqrt(no_nodes) * randu(no_nodes);
                network = join_horiz(xcoord, ycoord);
            } else { // generate regular lattice
                no_nodes = lattice_height*lattice_width;
                vec xcoord = linspace(0, lattice_width - 1, lattice_width);
                vec ycoord = linspace(0, lattice_height - 1, lattice_height);
                network.set_size(no_nodes, 2);
                int i, j, k = 0;
                for (j = 0; j < ycoord.n_rows; j++) { // grid expand algorithm
                    for (i = 0; i < xcoord.n_rows; i++) {
                        network(k, 0) = xcoord(i);
                        network(k, 1) = ycoord(j);
                        k++;
                    }
                }
            }

            if (adjMat.n_rows == 0) {
                genAdjMat(); // generate adjacency matrix
            }
            if (distMat.n_rows == 0) {
                genDistMat(); // generate distance matrix
            }
            if (adjMat.n_rows == 0) {
                genAdjMat(); // generate adjacency matrix
            }
        }
    } else { // import network
        cout << "\nImporting network " << xMatFileName << endl;
        mat X;
        X.load(xMatFileName);
        network = X;
        no_nodes = X.n_rows;
        genDistMat();
        genAdjMat();
    }
}

void Topography::genDistMat() {

    // summary:
    // generate euclidean distance matrix

    // required members:
    // network - cartesian coordinates of spatial network

    // output:
    // distMat - euclidean distance matrix

    mat xMat = repmat(network.col(0), 1, network.n_rows);
    mat yMat = repmat(network.col(1), 1, network.n_rows);
    distMat.set_size(network.n_rows, network.n_rows);
    distMat = sqrt(pow(xMat - trans(xMat),2) + pow(yMat - trans(yMat),2));

    mat SIGMA, UMat, sqrt_LAMBDA;  // objects for storing eigenvalues/vectors
    vec lambda;
    SIGMA = var_e * exp(-1 * distMat / phi); // covariance matrix

    eig_sym(lambda, UMat, SIGMA); // eigen decomposition
    lambda(find(lambda < 0)).zeros(); // remove numerical errors, should be no eigenvalues < 0
    sqrt_LAMBDA.zeros(distMat.n_rows, distMat.n_rows);
    sqrt_LAMBDA.diag() = sqrt(lambda); // diagonal matrix of square root eigenvalues

    sigEVec = UMat;
    sigEVal = sqrt_LAMBDA;
}

void Topography::genAdjMat() {

    // summary:
        // generate the adjacency matrix of the spatial network

    // required members:
        // gabriel - select edge assignement using Gabriel algorithm (https://en.wikipedia.org/wiki/Gabriel_graph) (T) or complete graph (F)

    // output:
        // adjMat - unweighted graph adjacency matrix

    adjMat.zeros(no_nodes,no_nodes); // initialize adjMat

    if (gabriel) {

        if (randGraph) {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// Gabriel algorithm  ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            auto distance = [] (double x1, double y1, double x2, double y2){
                return ((x1-x2)*(x1-x2)) + ((y1-y2)*(y1-y2));
            };

            int n_edges=0; // counter
            double mx, my, rad; // mean x, y, radius
            vec g1(3*no_nodes); // vector of nodes
            vec g2(3*no_nodes); // vector of nodes

            for(int i=0;i<no_nodes;i++) {
                // check, for each pair of nodes, if there is a third node within the circle whose diameter is defined by the line segment joining them
                for(int j=i+1;j<no_nodes;j++) {
                    mx=(network.col(0)[i] + network.col(0)[j])/2;
                    my=(network.col(1)[i] + network.col(1)[j])/2;
                    rad=distance(mx,my,network.col(0)[i], network.col(1)[i]);
                    int l;
                    for (l=0;l<no_nodes;l++) {
                        if((l!=i)&&(l!=j)&&(distance(mx,my,network.col(0)[l],network.col(1)[l])<rad)) {
                            break;
                        }
                    }
                    if (l==no_nodes) {
                        g1(n_edges)=i; g2(n_edges++)=j;
                    }
                }
            }

            for (int i=0; i<n_edges; i++) {
                adjMat(g1(i),g2(i)) = adjMat(g2(i),g1(i)) = 1.0;
            }
        } else { // lattice

//            unsigned int no_nodes_sqrt = sqrt(no_nodes);
            adjMat.set_size(no_nodes, no_nodes);
            for (int y=0; y<lattice_height; y++) {
                for (int x=0; x<lattice_width; x++) {
                    int i = y*lattice_width + x;
                    if (x > 0) {
                        adjMat(i-1,i) = 1;
                        adjMat(i,i-1) = 1;
                    }
                    if (y > 0) {
                        adjMat(i-lattice_width,i) = 1;
                        adjMat(i,i-lattice_width) = 1;
                    }
                }
            }

//            for (int k=1; k<=no_nodes; k++) {
//                int up = k - no_nodes_sqrt;
//                int down = k + no_nodes_sqrt;
//                int left = k - 1;
//                int right = k + 1;
//                if (up > 0) {
//                    adjMat(k-1, up-1) = 1.0;
//                }
//                if (down <= no_nodes_sqrt) {
//                    adjMat(k-1, down-1) = 1.0;
//                }
//                if (left % no_nodes_sqrt != 0) {
//                    adjMat(k-1, left-1) = 1.0;
//                }
//                if (k % no_nodes_sqrt != 0) {
//                    adjMat(k-1, right-1) = 1.0;
//                }
//            }

        }

    } else {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// Complete graph  ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        adjMat.ones();
        adjMat.diag().zeros(); // only diagonal elements set to zero

    }
}

void Topography::genEnvironment() {

    // summary:
        // generate an explicit environmental distribution of 1 or more spatially autocorrelated variables
        // modelled using a Gaussian random field generated via spectral decomposition (https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution)

    // required members:
        // envVar - number of explicitly modelled environmental variables
        // var_e - variance of the environmental distribution(s)
        // phi - spatial autocorrelation length of the environmental distribution(s)
        // distMat - euclidean distance matrix

    // output:
        // envMat - matrix of environmental distribution vectors, dimensions (envVar x no_nodes)

    if (envMatFileName.length() != 0) {
        cout << "\nImporting environment matrix " << envMatFileName << endl;
        envMat.load(envMatFileName);
    } else {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Sample from envVar random fields /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        vec zVec, eRow;
        envMat.set_size(envVar, no_nodes); // initialize enviroment object
        for (int i = 0; i < envVar; i++) {
            zVec.randn(no_nodes); // standard normal variables
            eRow = (sigEVec * sigEVal) * zVec; // sample from random field
            eRow = eRow - mean(eRow); // make standard normal
            eRow = eRow/stddev(eRow);
            envMat.row(i) = eRow.t(); // store in environment object
        }
    }
    rangeEnv.set_size(envVar);
    minEnv.set_size(envVar);

    for (int i = 0; i < envVar; i++) {
        rangeEnv(i) = envMat.row(i).max() - envMat.row(i).min();
        minEnv(i) = envMat.row(i).min();
    }
}

void Topography::genTempGrad() {
    // summary: ...
    envMat.reset();
    rowvec T_x(no_nodes);
    double T_grad = 1 / sqrt((double) no_nodes);
    T_x = T_int - T_grad * network.col(0).t();
    envMat = T_x; // define linear temperature gradient in single dimension
}

void Topography::genLandscape(mat netImprtd) {

    // summary:
        // generate a spatial network and domain decompose via recursive spectral bisection (https://en.wikipedia.org/wiki/Graph_partition)
        // minimum cost cut approach for MPI efficiency, subdomains fully connected and balanced

    // arguments:
        // netImprtd - imported network for simulation extension

    // required members:
        // network - cartesian coordinates of spatial network
        // distMat - euclidean distance matrix
        // adjMat - unweighted graph adjacency matrix

    // output:
        // abiotic landscape objects

    // external function calls:
        // genNetwork()
        // genDistMat()
        // genAdjMat()
        // genEnvironment()

    if (netImprtd.n_rows == 0) {
        genNetwork();
    }
    genDistMat();
    genAdjMat();
    if ((envVar != 0) && (envMatFileName.length() == 0) && envMat.n_rows == 0) { // generate enviroment
        cout << "\nGenerating abiotic environment" << endl;
        genEnvironment();
    }
    if (scMatFileName.length() != 0) { // import scaling matrix
        cout << "\nImporting scaling matrix " << scMatFileName << endl;
        mat Sc;
        Sc.load(scMatFileName);
        scVec_prime = Sc.col(0);
        scVec = Sc.col(1);
    }
}

// Local Variables:
// c-file-style: "stroustrup"
// End:
