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

FCFInfo::FCFInfo(double a, double b, double c, double alpha, double beta, double gamma, double dhighres, double F000, const std::vector<Int3x3> symops, int verbosity):
a_ (a), b_(b), c_(c), alpha_(alpha), beta_(beta), gamma_(gamma), dhighres_(dhighres), F000_(F000), 
        centrosymmetric_ (false),
        symops_(symops), 
        verbosity_(verbosity)
{

}

/**
 * Following sxfft.f, eliminate lattice and inversion operators
 * @return 
 */
std::vector<Int3x3> FCFInfo::symops_no_inv() {
    std::vector<Int3x3> symops_noinv(symops_);
    for (auto it1 = symops_noinv.begin(); it1 != symops_noinv.end(); ++it1) {
        for (auto it2 = it1 + 1; it2 != symops_noinv.end(); ++it2) {
            double u(0), v(0);
            // pairwise sum and difference of matrix entries
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    // u checks for lattice parameter, U=U' for rotational part
                    u += std::abs((*it2)(r, c) - (*it1)(r, c));
                    // V checks whether it2 has an equivalent V-V'=0
                    v += std::abs((*it2)(r, c) + (*it1)(r, c));
                }
            }
            // if neither u nor v is zero, it2 has no inversion element, nor is
            //  a pure lattice part
            if (std::min(u, v) > 0.01) continue; // next it2
                // both u and v greater than 0
                // sxfft.f sets NC here in case of centrosymmetric space groups,
                // this program does the separately with the function
                // centrosymmetric().
            else {
                // eliminate it2 by overwriting it with the last one
                if (verbosity_ > 2) {
                    std::cout << "---> debug: size of symmop vector: "
                            << symops_noinv.size() << '\n'
                            << " u = " << u << '\n'
                            << " v = " << v << '\n'
                            << *it1 << '\n'
                            << *it2 << std::endl;

                }
                // do not store it2, replace it with the last entry
                *it2 = symops_noinv.back();
                if (v < 0.01) {
                    centrosymmetric_ = true;
                }
                symops_noinv.pop_back();
                // repeat search for the back_operator
                --it2;
                continue;
            }
        }

    }
    if (verbosity_ > 1) {
        std::cout << "---> Elimination of symmetry and inversion operators \n"
                << "       results in " << "    " << symops_noinv.size() << " rotational operators\n";
    }
    return symops_noinv;

}

/**
 * this factor is used during FFT in sxfft.f as C(15); v scales the cell volume 
 * for non-orthorhombic cells, hence C(15) is the cell volume.
 * @return 
 */
double FCFInfo::fftscale() const {  
    // D1-3 contains sin; D4-6 is cos
    // sxfft.f computes C(7-12), but these do not seem to be used
    const double calpha = std::cos(M_PI / 180.0 * alpha_);
    const double cbeta  = std::cos(M_PI / 180.0 * beta_);
    const double cgamma = std::cos(M_PI / 180.0 * gamma_);
    
    const double v = 1.0 - calpha*calpha - cbeta*cbeta - cgamma*cgamma +
                    2.0 * calpha*cbeta*cgamma;
    if (verbosity_ > 2) {
    	std::cout << "---> debugging: value of v = " << v << '\n';
	}

    double c15 = a_*b_*c_*std::sqrt(v);
    
    return c15;
}
