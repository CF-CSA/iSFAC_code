/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Utils.h
 * Author: tg
 *
 * Created on October 11, 2024, 2:12 PM
 */

#ifndef UTILS_H
#define UTILS_H

#include "Vec3.h"
#include "Mat33.h"
#include "Atom.h"

#include <cmath>
#include <gsl/gsl_matrix_double.h>
#include <tuple>



namespace Utils {
    // create gridpoints for van-der-Waals volume
    std::vector<Vec3> vdw_vol_grid(const std::vector<Atom>& atoms, double gridspacing, int verbosity);

    //! translate set of coordinates to their centroid
    std::vector<Vec3> centroid(const std::vector<Vec3>& coords);

    template <typename T> T CC(const std::vector<T>&, const std::vector<T>&);


    //! compute bounding box.
    template <typename T> std::pair<Vec3, Vec3> bbox3D(const std::vector<T> coords, const unsigned short& verbosity);
    
    //! compute unit cell vectors from parameters
    std::tuple<Vec3, Vec3, Vec3> unit_cell_vector (double a, double b, double c, double alpha, double beta, double gamma);

    std::string prompt(const unsigned short& verbosity);

}

/**
 * Calculates the correlation coefficient for two lists of numbers of type T.
 * It is not checked whether their sizes match
 * \param l1    first  list of numbers
 * \param l2    second list of numbers
 */
template <typename T> T Utils::CC(const std::vector<T>& l1,
        const std::vector<T>& l2) {
    T sum1, sum1_sq, sum2, sum2_sq, sum12;
    T nominator;
    T denominator;

    sum1 = sum2 = sum1_sq = sum2_sq = sum12 = 0;

    double inv_list_size = 1.0 / l1.size();

    for (unsigned int i = 0; i < l1.size(); i++) {
        sum1 += l1[i];
        sum2 += l2[i];
        sum1_sq += l1[i] * l1[i];
        sum2_sq += l2[i] * l2[i];
        sum12 += l1[i] * l2[i];
    }

    nominator = sum12 - inv_list_size * sum1 * sum2;
    denominator = (sum1_sq - inv_list_size * sum1 * sum1) *
            (sum2_sq - inv_list_size * sum2 * sum2);

    return (nominator / std::sqrt(denominator));
}


/**
 * compute 3D bounding box of @c coords. T must provide x(), y(), z()
 * @param coords
 * @param verbosity
 * @return llc, urc pair of lower left and upper right corner
 */
template <typename T> std::pair<Vec3, Vec3> Utils::bbox3D(const std::vector<T> coords, const unsigned short& verbosity) {
    const T a = coords.front();
    double minx = a.x();
    double miny = a.y();
    double minz = a.z();
    double maxx = a.x();
    double maxy = a.y();
    double maxz = a.z();
    for (auto x: coords) {
        minx = std::min(x.x(), minx);
        miny = std::min(x.y(), miny);
        minz = std::min(x.z(), minz);
        maxx = std::max(x.x(), maxx);
        maxy = std::max(x.y(), maxy);
        maxz = std::max(x.z(), maxz);
    }
    const Vec3 llc (minx, miny, minz);
    const Vec3 urc (maxx, maxy, maxz);
    
    return std::pair<Vec3, Vec3>(llc, urc);
}

#endif /* UTILS_H */

