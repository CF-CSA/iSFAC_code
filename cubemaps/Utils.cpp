/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <regex>

#include "Utils.h"
#include <iostream>
#include <numeric>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
/**
 * algorithm to create a surface grid at vdW radii of atoms, with a grid spacing 
 * og @c gridspacing (in Angstrom)
 * - creeate a bounding box at @c gridspacing
 * - for each grid point check distance to atoms, mark when whithin vdW radius of any
 * algortihm from ChatGPT
 * @param atoms
 * @param gridspacing
 * @return 
 */
std::vector<Vec3> Utils::surfacegrid (const std::vector<Atom>& atoms, double gridspacing, int verbosity) {
    
    Vec3 bbox_min, bbox_max;
    std::vector<Vec3> grid, surface;
    double minx(9999), miny(9999), minz(9999), maxx(-9999), maxy(-9999), maxz(-9999);

    // find min and max of coordinates
    for (auto atom : atoms) {
        if ((atom.xyz()).x()-atom.vdw_radius() < minx) minx = (atom.xyz()).x()-atom.vdw_radius() ;
        if ((atom.xyz()).x()+atom.vdw_radius() > maxx) maxx = (atom.xyz()).x()+atom.vdw_radius() ;

        if ((atom.xyz()).y()-atom.vdw_radius() < miny) miny = (atom.xyz()).y()-atom.vdw_radius();
        if ((atom.xyz()).y()+atom.vdw_radius() > maxy) maxy = (atom.xyz()).y()+atom.vdw_radius();
        
        if ((atom.xyz()).z()-atom.vdw_radius() < minz) minz = (atom.xyz()).z()-atom.vdw_radius();
        if ((atom.xyz()).z()+atom.vdw_radius() > maxz) maxz = (atom.xyz()).z()+atom.vdw_radius();
    }
    bbox_min = Vec3(minx, miny, minz);
    bbox_max = Vec3(maxx, maxy, maxz);
    
    int stepsx = (maxx-minx)/gridspacing;
    int stepsy = (maxy-miny)/gridspacing;
    int stepsz = (maxz-minz)/gridspacing;
    
    //! create list of grid points
    for (int ix = 0; ix < stepsx; ++ix){
        for (int iy = 0; iy < stepsy; ++iy) {
            for (int iz = 0; iz < stepsz; ++iz) {
                Vec3 xyz (minx + ix*gridspacing, miny + iy*gridspacing, minz + iz*gridspacing);
                grid.push_back(xyz);
                for (auto atom : atoms) {
                    if (atom.insphere(xyz)) {
                        surface.push_back(xyz);
                        break;
                    }
                }
            }
        }
    }
    if (verbosity > 1) {
        std::cout << "---> Generated bounding box with " << grid.size() << " and surface with " 
                << surface.size() << " grid points\n";
    }
    
    return surface;
    
}


/**
 * compute the vector that moves the centroid of @c moved to the centroid of @c fixed
 * @param fixed reference set of coordinates
 * @param moved moving set of coordinates
 * @return translation vector
 */
std::vector<Vec3> Utils::centroid (const std::vector<Vec3>& coords) {
    Vec3 ctr = 1./coords.size()*std::accumulate(coords.begin(), coords.end(), Vec3(0,0,0));
    std::vector<Vec3> centred (coords);
    for (auto &x: centred) {
        x = x-ctr;
    }
    
    return coords;
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
    H = gsl_matrix_alloc(N,N);
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            double v(0.0);
            // compute matrix multiplication H = moved^T * fixed
            v += moved[r].x() * fixed[c].x() + moved[r].y() * fixed[c].y() + moved[r].z() * fixed[c].z();
            gsl_matrix_set(H, r, c, v);
            // now compute SVD of H
        }
    }
    gsl_matrix* V = gsl_matrix_alloc(N,N);
    gsl_vector* S = gsl_vector_alloc(N);
    gsl_vector* work = gsl_vector_alloc(N);
    
    int result = gsl_linalg_SV_decomp(H, V, S, work);
 
    
    /*
     SVD of H results in U Sigma V^T
     * R = V U^T where H is stored in H
     * V should be a 3xN matrix, U^T a Nx3 matrix
     * R_ij = sum_k V_ik U^T _kj = sum_k V_ik*U_jk
     */
    Mat33 R;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; ++j) {
            double r(0.0);
            for (int k = 0; k < 3; ++k){
                r += gsl_matrix_get(V, i, k)*gsl_matrix_get(H, j, k);
            }
            R(i,j) = r;
        }
    }
    
    return R;
}
