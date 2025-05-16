/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FCFInfo.cpp
 * Author: tg
 * 
 * Created on July 30, 2023, 9:51 PM
 */

#include "FCFInfo.h"
#include "Utils.h"
#include <cmath>
#include <complex>

FCFInfo::FCFInfo(float a, float b, float c, float alpha, float beta, float gamma, float dhighres, float F000, const std::vector<Int3x3> symops, int verbosity):
a_ (a), b_(b), c_(c), alpha_(alpha), beta_(beta), gamma_(gamma), dhighres_(dhighres), F000_(F000), symops_(symops), verbosity_(verbosity)
{
    centrosymmetric_ = check_inversion();
}

/**
 * the space group has a centre of inversion, if any of the rotation matrices has 
 * a negative determinant
 * @return 
 */
bool FCFInfo::check_inversion() const{
    bool c = false;
    for (auto it = symops_.begin(); it != symops_.end(); ++it) {
        if (verbosity_ > 2) {
            float det = it->det();
            std::cout << Utils::prompt(2) << "Checking Determinant for \n" << (*it) 
                    << "---> Det(R) = " << det << "\n\n";
        }
        if (it->det() == -1) {
            c = true;
            continue;
        }
        if (it->det() != 1) {
            std::cout << "*** Error: Determinant of rotation matrix " << *it 
                    << " != +/-1" << std::endl;
            throw std::logic_error("Error: Determinan of symmetry matrix must be +/-1");
        }
    }
    return c;
}

/**
 * Following sxfft.f, eliminate lattice and inversion operators
 * @return 
 */
std::vector<Int3x3> FCFInfo::symops_no_inv() const {
    
    std::vector<Int3x3> symops_noinv (symops_);
    for (auto it1 = symops_noinv.begin(); it1 != symops_noinv.end(); ++it1) {
        for (auto it2 = it1+1; it2 != symops_noinv.end(); ++it2){
            int u(0), v(0);
            // pairwise sum and difference of matrix entries
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    u += std::abs((*it2)(r, c) - (*it1)(r,c));
                    v += std::abs((*it2)(r, c) + (*it1)(r,c));
                }
            }
            if (u != v) continue; // next it2
            else {
                // eliminate it2 by overwriting it with the last one
                *it2 = symops_noinv.back();
                symops_noinv.pop_back();
            }
        }
            
    }
    if (verbosity_ > 1) {
        std::cout << "---> Elimination of symmetry and inversion operators \n"
                << "       results in "  << "    " << symops_noinv.size() << " rotational operators\n";
    }
    return symops_noinv;
    
}

/**
 * this factor is used during FFT in sxfft.f as C(15); not sure of its meaning
 * @return 
 */
float FCFInfo::fftscale() const {  
    const float cdalpha = std::cos(M_PI / 180.0 * alpha_);
    const float cdbeta  = std::cos(M_PI / 180.0 * beta_);
    const float cdgamma = std::cos(M_PI / 180.0 * gamma_);
    
    const float v = 1.0 - cdalpha*cdalpha - cdbeta*cdbeta - cdgamma*cdgamma +
                    2.0 * cdalpha*cdbeta*cdgamma;

    float c15 = a_*b_*c_*std::sqrt(v);
    
    return c15;
    
}