/* 
 * File:   Reflex.cpp
 * Author: tg
 * 
 * Created on July 3, 2012, 4:30 PM
 */

#include "Reflex.h"
#include <cmath>
#include <complex>
#include <iomanip>

Reflex::Reflex()
: hkl_(0,0,0), F_ (0.0), phi_(0.0){
}

Reflex::Reflex(const HKL& hkl, double F,double phi)
        : hkl_(hkl), F_(F),phi_(phi)
{}

Reflex::~Reflex() {
}

/**
 * calculates the phase shift 2\pi h^T t from a symmetry operator
 * S = R+t
 * @param symops
 */
double Reflex::phaseshift(const Int3x3& symop) const{
    double shift (1.0);
    shift = 1.0* hkl_.h() * symop(0)
            + 1.0* hkl_.k() * symop(1)
            + 1.0*hkl_.l() * symop(2);
    shift *= 2.0*M_PI;
    return shift;
}

/**
 * Transform @c reflex such that h,k,l have minimal absolute value, 
 * and transform its phase accordingly \$f phi' = phi + 2\pi i h^t t\f$
 * @param symops
 */
void Reflex::stdsetting(const std::vector<Int3x3>& symops) {
    std::vector<Int3x3>::const_iterator it;
    
    HKL minHKL (hkl_);
    
    
    double phi (M_PI/180. * phi_);
    std::complex<double> phase = std::polar<double>(1.0, phi);
    std::complex<double> minphase(phase);
    
    for (it=symops.begin(); it != symops.end(); ++it) {
        HKL hkl = *it * hkl_;
        if (hkl < minHKL) {
            minHKL = hkl;
            double dphi = phaseshift(*it);
            std::complex<double> dphase (std::cos(dphi), std::sin(dphi));
            minphase = phase * dphase;
        }
    }
    hkl_ = minHKL;
    phi_ = 180.0 / M_PI * std::arg(minphase);
}

double Reflex::dstarsq(const Vec3& astar, const Vec3& bstar, const Vec3& cstar) const {
    return hkl_.dstarsq( astar, bstar, cstar);
}

/* comparison based on indices */
bool operator <  (const Reflex& R1, const Reflex& R2) {
    return R1.hkl_ <  R2.hkl_;
}

/* equality based on indices */
bool operator == (const Reflex& R1, const Reflex& R2) {
    return R1.hkl_ == R2.hkl_;
}

std::ostream& operator<< (std::ostream& outp, const Reflex& R) {
    outp    << std::setw(3) << R.hkl_.h()
            << std::setw(3) << R.hkl_.k()
            << std::setw(3) << R.hkl_.l()
            << std::fixed
            << std::setw(12)
            << std::setprecision(3)
            << R.F_
            << std::fixed
            << std::setw(12)
            << std::setprecision(3)
            << R.phi_;
    return outp;            
}