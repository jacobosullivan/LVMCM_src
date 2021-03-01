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
                vec xcoord;
                xcoord = linspace(0, sqrt(no_nodes) - 1, sqrt(no_nodes));
                network.set_size(no_nodes, 2);
                int i, j, k = 0;
                for (i = 0; i < xcoord.n_rows; i++) { // grid expand algorithm
                    for (j = 0; j < xcoord.n_rows; j++) {
                        network(k, 0) = xcoord(i);
                        network(k, 1) = xcoord(j);
                        k++;
                    }
                }
            }

            if (adjMat.n_rows == 0) {
                genAdjMat(); // generate adjacency matrix
            }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// List nodes in order of vector norm  ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            vec norm(no_nodes); // generate vector of vector norms of nodes for network ordering
            for (int i = 0; i < no_nodes; i++) {
                norm(i) = sqrt(pow(network(i, 0) + network(i, 1), 2));
            }
            network = join_horiz(network, norm);

            int end = no_nodes - 1;
            rowvec swap;
            colvec swap_col;
//        for (int j = end - 1; j > 0; j--) { // vector ordering algorithm
//            for (int i = 0; i < end; i++) {
//                if (network(i, 2) > network(i + 1, 2)) {
//                     reorder network
//                    swap = network.row(i + 1);
//                    network.row(i + 1) = network.row(i);
//                    network.row(i) = swap;

            // reorder adjacency - THIS NEEDS TO BE DEBUGGED
//                    swap = adjMat.row(i + 1);
//                    swap_col = adjMat.col(i + 1);
//                    adjMat.row(i + 1) = adjMat.row(i);
//                    adjMat.col(i + 1) = adjMat.col(i);
//                    adjMat.row(i) = swap;
//                    adjMat.col(i) = swap_col;
//                }
//            }
//            end--;
//        }
            network.shed_col(2);

            if (distMat.n_rows == 0) {
                genDistMat(); // generate distance matrix
            }
            if (adjMat.n_rows == 0) {
                genAdjMat(); // generate adjacency matrix
            }

            fVec.zeros(no_nodes); // initialize indicator vector
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

            unsigned int no_nodes_sqrt = sqrt(no_nodes);

            for (int k=1; k<=no_nodes; k++) {
                int up = k - no_nodes_sqrt;
                int down = k + no_nodes_sqrt;
                int left = k - 1;
                int right = k + 1;
                if (up > 0) {
                    adjMat(k-1, up-1) = 1.0;
                }
                if (down <= no_nodes_sqrt) {
                    adjMat(k-1, down-1) = 1.0;
                }
                if (left % no_nodes_sqrt != 0) {
                    adjMat(k-1, left-1) = 1.0;
                }
                if (k % no_nodes_sqrt != 0) {
                    adjMat(k-1, right-1) = 1.0;
                }
            }

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Sample from envVar random fields /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vec zVec, eRow;
    envMat.set_size(envVar, no_nodes); // initialize enviroment object
    for (int i=0; i<envVar; i++) {
        zVec.randn(no_nodes); // random uniform variables
        eRow = (sigEVec * sigEVal) * zVec; // sample from random field
//        eRow -= eRow.min(); // translate e_min to 0
//        eRow /= eRow.max(); // rescale to range [0,1]
        envMat.row(i) = eRow.t(); // store in environment object
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

void Topography::genDomainDecomp(mat netImprtd) {

    // summary:
        // generate a spatial network and domain decompose via recursive spectral bisection (https://en.wikipedia.org/wiki/Graph_partition)
        // minimum cost cut approach for MPI efficiency, subdomains fully connected and balanced

    // arguments:
        // netImprtd - imported network for simulation extension

    // required members:
        // bisec - number of recursive bisections
        // network - cartesian coordinates of spatial network
        // distMat - euclidean distance matrix
        // adjMat - unweighted graph adjacency matrix

    // output:
        // fVec - indicator vector, elements of which denote subdomain allocation of nodes

    // external function calls:
        // genNetwork()
        // genDistMat()
        // genAdjMat()
        // genEnvironment()

    if (bisec > 0) {
        int no_sub = pow(2, bisec); // final number of sumdomains
        printf("\nClustering %d nodes into %d subdomains of %d\n", no_nodes, no_sub, (int) no_nodes / no_sub);
        uvec f_temp, fUnq; // indicator vectors of subdomains
        bool stop; // flag for restarting do loop if any subdomain is not connected
        int unb = 0, unc = 0, attempt = 0; // counters for recording number of times do loop repeated due to unbalanced/unconnected subdomains

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Spectral bisection algorithm ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        do {
            stop = false;
            attempt++;
            fVec.reset();
            if (netImprtd.n_rows == 0) {
                network.reset();
                distMat.reset();
                adjMat.reset();
                printf("\rGenerating cartesian coordinates and distance matrix, attempt %d ... ", attempt);
                genNetwork();
                if (envVar != 0) { // generate enviroment
                    genEnvironment();
                }
                if (T_int != -1.0) {
                    genTempGrad();
                }
            } else {
                if (attempt == 1) {
                    cout << "\nDecomposing imported network... ";
                    genDistMat();
                    genAdjMat();
                } else {
                    cout << "\nDecomposition of imported network failed, check bisec parameter" << endl;
                    stop = true;
                    return;
                }
            }

            double medField; // median of Fielder vector
            mat lMat, dMat, wMat; // matrices used in computation of graph Laplacian
            mat adjMatSub, distMatSub; // storage for recursive subsetting of adjacency and distance matrices
            uvec index; // index for subsetting whole domain
            vec fieldRelx; // Fielder vector storage object
            uvec fieldDscr; // Discretized Fielder vector
            f_temp.zeros(no_nodes);

            for (int b = 0; b < bisec; b++) {
                fUnq = unique(f_temp); // current clusters
                for (int i = 0; i < fUnq.n_rows; i++) {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// Generate graph Laplician //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    index = find(f_temp == fUnq(i)); // indices of subdomain i
                    distMatSub = distMat.submat(index, index); // subset distance matrix
                    wMat = distMatSub;
                    wMat(find(wMat)) //locate non-zero elements (off-d)
                            = pow(wMat(find(wMat)), -1); // edge weights defined as inverse of distance for cut selection
                    adjMatSub = adjMat.submat(index, index);
                    wMat %= adjMatSub; // weighted adjacency matrix
                    dMat.zeros(index.n_rows, index.n_rows);
                    dMat.diag() = sum(wMat, 1); // diagonal weighted degree matrix
                    lMat = dMat - wMat; // graph Laplacian

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Compute Fielder vector and use to bisec (sub)domain /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    vec eigval; // eigenvalue storage object
                    mat eigvec; // eigenvector storage object
                    eig_sym(eigval, eigvec, lMat); // eigen decomposition
                    fieldRelx = eigvec.col(1);
                    medField = median(fieldRelx);
                    fieldDscr.zeros(fieldRelx.n_rows);
                    fieldDscr(find(fieldRelx >= medField)).ones(); // discretize Fielder vector by median
                    f_temp(find(f_temp == fUnq(i))) = 1 + max(f_temp) + fieldDscr; // define new factor level
                    if (sum(fieldDscr) != fieldDscr.n_rows / 2) { // check new subdomains are balanced
                        printf("Unbalanced subdomain");
                        fflush(stdout);
                        unb++;
                        stop = true;
                        break;
                    }

                    for (int j=0; j<2; j++) { // check new subdomains each connected

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Compute subdomain graph Laplacian and check lambda[2] > 0 ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        index = find(fieldDscr == j); // ambiguous redeclaration
                        wMat = distMatSub.submat(index, index); // subset distance matrix
                        wMat(find(wMat)) //locate non-zero elements (off-d)
                                = pow(wMat(find(wMat)), -1); // edge weights defined as inverse of distance for cut selection
                        wMat %= adjMatSub.submat(index, index); // weighted adjacency matrix
                        dMat.zeros(index.n_rows, index.n_rows);
                        dMat.diag() = sum(wMat, 1); // diagonal weighted degree matrix
                        lMat = dMat - wMat; // graph Laplacian
                        eig_sym(eigval, eigvec, lMat);
                        index = find(eigval < 1e-10); // search for zero eigenvalues
                        if (index.n_rows > 1) {
                            printf("Unconnected subdomain");
                            fflush(stdout);
                            unc++;
                            stop = true;
                            break;
                        }
                    }

                    if (stop) { break; }
                }

                if (stop) { break; }
            }

            if (stop) { continue; }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Reset fVec from zero and order network by decomposition ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            cout << "partitioning successful" << endl;
            f_temp -= min(f_temp); // reset numbering to 0:(2^sppPool.topo.bisec-1)
            network.resize(network.n_rows, network.n_cols+1);
            network = network.rows(sort_index(f_temp)); // reorder network by decomposition
            f_temp = f_temp.elem(sort_index(f_temp));
            network = network.cols(0,1);
            fVec = f_temp; // save indicator vector
            uvec unq = unique(fVec);
            uvec count;
            count.zeros(unq.n_rows);
            cFVec.set_size(fVec.n_rows);
            for (int i = 0; i < fVec.n_rows; i++) {
                cFVec(i) = count(fVec(i));
                count(fVec(i))++;
            }
            printf("Graph reseeded %d/%d times due to unbalanced/unconnected subdomains\n", unb, unc);
        } while (stop);
    } else {
        if (netImprtd.n_rows == 0) {
            network.reset();
            distMat.reset();
            adjMat.reset();
            printf("\rGenerating cartesian coordinates and distance matrix without decomposition... ");
            genNetwork();
            if (envVar != 0) { // generate enviroment
                genEnvironment();
            }
        } else {
            genDistMat();
            genAdjMat();
        }
        if (scMatFileName.length() != 0) { // import scaling matrix
            cout << "\nImporting scaling matrix " << scMatFileName << endl;
            mat Sc;
            Sc.load(scMatFileName);
            scVec_prime = Sc.col(0);
            scVec = Sc.col(1);
        }
    }
}

// Local Variables:
// c-file-style: "stroustrup"
// End:
