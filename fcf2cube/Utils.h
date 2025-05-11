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



namespace Utils {
    // create gridpoints for van-der-Waals volume
std::vector<Vec3> vdw_vol_grid (const std::vector<Atom>& atoms, double gridspacing, int verbosity);

//! translate set of coordinates to their centroid
std::vector<Vec3> centroid (const std::vector<Vec3>& coords);

template <typename T> T CC(const std::vector<T>&, const std::vector<T>&);

std::string prompt(const unsigned short& verbosity);

}
/**
 * Calculates the correlation coefficient for two lists of numbers of type T.
 * It is not checked whether their sizes match
 * \param l1    first  list of numbers
 * \param l2    second list of numbers
 */
template <typename T> T Utils::CC(const std::vector<T>& l1,
                                        const std::vector<T>& l2)
{
    T sum1, sum1_sq, sum2, sum2_sq, sum12;
    T nominator;
    T denominator;

    sum1 = sum2 = sum1_sq = sum2_sq = sum12 = 0;

    double inv_list_size = 1.0/l1.size();

    for (unsigned int i = 0; i < l1.size(); i++)
    {
        sum1    += l1[i];
        sum2    += l2[i];
        sum1_sq += l1[i]*l1[i];
        sum2_sq += l2[i]*l2[i];
        sum12   += l1[i]*l2[i];
    }

    nominator   = sum12 - inv_list_size*sum1 * sum2;
    denominator = (sum1_sq - inv_list_size * sum1*sum1) *
                  (sum2_sq - inv_list_size * sum2*sum2);

    return (nominator / std::sqrt(denominator));
}

#endif /* UTILS_H */

