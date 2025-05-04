/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <regex>

#include "Utils.h"
#include <iostream>
#include <numeric>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>

/**
 * algorithm to create a volume grid at vdW radii of atoms, with a grid spacing 
 * og @c gridspacing (in Angstrom)
 * - create a bounding box at @c gridspacing
 * - for each grid point check distance to atoms, mark when whithin vdW radius of any
 * algortihm from ChatGPT
 * @param atoms
 * @param gridspacing
 * @return 
 */
std::vector<Vec3> Utils::vdw_vol_grid(const std::vector<Atom>& atoms, double gridspacing, int verbosity) {

    Vec3 bbox_min, bbox_max;
    std::vector<Vec3> grid, vdw_volume;
    double minx(9999), miny(9999), minz(9999), maxx(-9999), maxy(-9999), maxz(-9999);

    // find min and max of coordinates
    for (auto atom : atoms) {
        minx = std::min(minx, (atom.xyz()).x() - atom.vdw_radius());
        maxx = std::max(maxx, (atom.xyz()).x() + atom.vdw_radius());

        miny = std::min(miny, (atom.xyz()).y() - atom.vdw_radius());
        maxy = std::max(maxy, (atom.xyz()).y() + atom.vdw_radius());

        minz = std::min(minz, (atom.xyz()).z() - atom.vdw_radius());
        maxz = std::max(maxz, (atom.xyz()).z() + atom.vdw_radius());
    }
    bbox_min = Vec3(minx, miny, minz);
    bbox_max = Vec3(maxx, maxy, maxz);

    int stepsx = (maxx - minx) / gridspacing;
    int stepsy = (maxy - miny) / gridspacing;
    int stepsz = (maxz - minz) / gridspacing;

    //! create list of grid points
    if (verbosity > 2) {
        std::cout << Utils::prompt(2) << "Grid points for VdW surface:\n";
    }
    for (int ix = 0; ix < stepsx; ++ix) {
        for (int iy = 0; iy < stepsy; ++iy) {
            for (int iz = 0; iz < stepsz; ++iz) {
                Vec3 xyz(minx + ix*gridspacing, miny + iy*gridspacing, minz + iz * gridspacing);
                grid.push_back(xyz);
                for (auto atom : atoms) {
                    if (atom.insphere(xyz)) {
                        vdw_volume.push_back(xyz);
                        if (verbosity > 3) {
                            std::cout << xyz << '\n';
                        }
                        break;
                    }
                }
            }
        }
    }
    if (verbosity > 0) {
        std::cout << Utils::prompt(0) << "Generated bounding box with "
                << grid.size() << " points and VdW volume with "
                << vdw_volume.size() << " points\n"
                << Utils::prompt(0) << "Grid spacing: " << gridspacing << "A\n";
    }

    return vdw_volume;
}

/**
 * compute the vector that moves the centroid of @c moved to the centroid of @c fixed
 * @param fixed reference set of coordinates
 * @param moved moving set of coordinates
 * @return translation vector
 */
std::vector<Vec3> Utils::centroid(const std::vector<Vec3>& coords) {
    Vec3 ctr = 1. / coords.size() * std::accumulate(coords.begin(), coords.end(), Vec3(0, 0, 0));
    std::vector<Vec3> centred(coords);
    for (auto &x : centred) {
        x = x - ctr;
    }

    return coords;
}

/**
 * Compute the symmetric fill matrix of intramolecular distances for a set 
 * f coordinates
 * @param coords
 * @return 
 */
std::vector<double> Utils::distance_matrix(const std::vector<Vec3>& coords) {
    std::vector<double> distances;
    const int N = (coords.size() * coords.size());
    distances.reserve(N);

    for (std::vector<Vec3>::const_iterator itr = coords.begin(); itr != coords.end() - 1; ++itr) {
        for (std::vector<Vec3>::const_iterator itc = coords.begin(); itc != coords.end(); ++itc) {
            double dsqd = ((*itr) - (*itc)).lengthsq();
            distances.push_back(std::sqrt(dsqd));
        }
    }
    return distances;
}

/**
 * Compute rotation matrix that rotates moved onto fixed
 * follows  Kabsch algortithm with SVD as explained at 
 * https://en.wikipedia.org/wiki/Kabsch_algorithm
 * atoms supposed to be in aligned order and already moved to their respective centroids
 * @param fixed
 * @param moved
 * @return 
 */
Mat33 Utils::KabschR(const std::vector<Vec3>& fixed, const std::vector<Vec3> moved) {
    // prepare for H = fixed^T * moved
    gsl_matrix* H;
    int N = fixed.size();
    H = gsl_matrix_alloc(N, N);
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            double v(0.0);
            // compute matrix multiplication H = moved^T * fixed
            v += moved[r].x() * fixed[c].x() + moved[r].y() * fixed[c].y() + moved[r].z() * fixed[c].z();
            gsl_matrix_set(H, r, c, v);
            // now compute SVD of H
        }
    }
    gsl_matrix* V = gsl_matrix_alloc(N, N);
    gsl_vector* S = gsl_vector_alloc(N);
    gsl_vector* work = gsl_vector_alloc(N);

    std::cout << "---> Before SVD for H, dimensions " << H->size1 << ' ' << H->size2 << '\n';
    for (int r = 0; r < H->size1; ++r) {
        std::cout << "   ";
        for (int c = 0; c < H->size2; ++c) {
            std::cout << gsl_matrix_get(H, r, c) << ' ';

        }
        std::cout << '\n';
    }


    int result = gsl_linalg_SV_decomp(H, V, S, work);

    std::cout << "---> After SVD for H, dimensions " << H->size1 << ' ' << H->size2 << '\n';
    for (int r = 0; r < H->size1; ++r) {
        std::cout << "   ";
        for (int c = 0; c < H->size2; ++c) {
            std::cout << gsl_matrix_get(H, r, c) << ' ';

        }
        std::cout << '\n';
    }

    double det = Utils::determinant(H) * Utils::determinant(V);
    std::cout << "---> Determinant of U*VT = " << det << '\n';
    if (det < 0.0) {
        // flip sign of last column of U^T
        // = last row of U^T
        int row = H->size1 - 1;
        for (int c = 0; c < H->size2; ++c) {
            double h = -1 * gsl_matrix_get(H, row, c);
            gsl_matrix_set(H, row, c, h);
        }
    }


    /*
     SVD of H results in U Sigma V^T
     * R = V U^T where U is stored in H
     * V should be a 3xN matrix, U^T a Nx3 matrix
     * R_ij = sum_k V_ik U^T _kj = sum_k V_ik*U_jk
     */
    Mat33 R;
    gsl_matrix* rgsl = gsl_matrix_alloc(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; ++j) {
            double r(0.0);
            for (int k = 0; k < 3; ++k) {
                r += gsl_matrix_get(V, i, k) * gsl_matrix_get(H, j, k);
            }
            R(i, j) = r;
            gsl_matrix_set(rgsl, i, j, r);
        }
    }

    det = Utils::determinant(rgsl);
    std::cout << "---> Determinant of Rgsl = " << det << '\n';
    for (int r = 0; r < rgsl->size1; ++r) {
        std::cout << "   ";
        for (int c = 0; c < rgsl->size2; ++c) {
            std::cout << gsl_matrix_get(rgsl, r, c) << ' ';

        }
        std::cout << '\n';
    }

    gsl_matrix_free(H);
    gsl_matrix_free(V);
    gsl_vector_free(S);

    det = R.determinant();
    std::cout << Utils::prompt(1) << "Determinant of R = \n" << R << '\n' << det << '\n';

    return R;
}

/**
 * Computes the determinant of M , assumed to be square
 * based on gsl LU decomposition
 * @param M
 * @return 
 */
double Utils::determinant(const gsl_matrix* M) {
    int signum;
    gsl_permutation * p;
    double det = 1.0;

    // make a copy of B to call LU_decomp
    gsl_matrix* B = gsl_matrix_alloc(M->size1, M->size2);
    gsl_matrix_memcpy(B, M);

    p = gsl_permutation_alloc(B->size1);

    gsl_linalg_LU_decomp(B, p, &signum);
    gsl_linalg_LU_det(B, signum);

    gsl_matrix_free(B);
    gsl_permutation_free(p);
    return det;

}

double Utils::CC_gsl(const std::vector<double>& d1, const std::vector<double>& d2) {
    gsl_vector* g1 = gsl_vector_alloc(d1.size());
    gsl_vector* g2 = gsl_vector_alloc(d1.size());

    for (size_t i = 0; i < d1.size(); ++i) {
        gsl_vector_set(g1, i, d1[i]);
        gsl_vector_set(g2, i, d2[i]);
    }

    const double cc = gsl_stats_correlation(&(d1[0]), 1, &(d2[0]), 1, d1.size());

    gsl_vector_free(g1);
    gsl_vector_free(g2);
    return cc;
}

/**
 * return the prompt arrow ---> which 2+ verbosity number of dashes
 * @param verbosity
 * @return 
 */
std::string Utils::prompt(const unsigned short& verbosity) {
    std::string p("#--");
    for (int i = 0; i < verbosity; ++i) {
        p += '-';
    }
    p += "> ";
    return p;

}