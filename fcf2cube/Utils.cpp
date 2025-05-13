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

/**
 * 
 * @param a
 * @param b
 * @param c
 * @param alpha angle in degree
 * @param beta  angle in degree
 * @param gamma angle in degree
 * @return 
 */
std::tuple<Vec3, Vec3, Vec3> Utils::unit_cell_vector(double a, double b, double c, double alpha, double beta, double gamma) {
        const float salpha = std::sin(alpha * M_PI / 180.0);
    const float calpha = std::cos(alpha * M_PI / 180.0);
    const float sbeta  = std::sin(beta  * M_PI / 180.0);
    const float cbeta  = std::cos(beta  * M_PI / 180.0);
    const float sgamma = std::sin(gamma * M_PI / 180.0);
    const float cgamma = std::cos(gamma * M_PI / 180.0);
    Vec3 A = Vec3(a, 0.0, 0.0);
    Vec3 B = Vec3(b*cgamma, b*sgamma, 0);
        double cx = c * cbeta;
    double cy = (c * calpha - cx * cgamma) / sgamma;
    double cz = std::sqrt(c * c - cx * cx - cy * cy);
    Vec3 C = Vec3(cx, cy, cz);
    return std::tuple<Vec3, Vec3, Vec3>(A, B, C);
}