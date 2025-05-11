#include "HKL.h"
#include "HKLops.h"
#include "Int3x3.h"
#include "FCFitem.h"
// #include <cstdlib>
#include <iostream>
#include <iomanip>

/**
 * Comparison of reflections runs in the order h, k, l
 * i.e. (0 1 2) < (1 0 0) etc.
 * @param h1
 * @param h2
 * @return 
 */
bool operator< (const HKL& h1, const HKL& h2) {
    if (std::abs(h1.h_) < std::abs(h2.h_)) return true;
    else if (std::abs(h1.h_) > std::abs(h2.h_) ) return false;
    else if (std::abs(h1.k_) < std::abs(h2.k_) ) return true;
    else if (std::abs(h1.k_) > std::abs(h2.k_) ) return false;
    else return (std::abs(h1.l_) < std::abs(h2.l_) );
}

/**
 * Equality of (h k l) if all indices identical
 */
bool operator==(const HKL& h1, const HKL& h2) {
    return (h1.h_ == h2.h_ && h1.k_ == h2.k_ && h1.l_ == h2.l_);
}

/**
 * inverse comparison, see operator <
 */
bool operator >(const HKL& h1, const HKL& h2) {
    return (h2 < h1);
}

/**
 * add two reflections 
 */
HKL operator+ (const HKL& h1, const HKL& h2) {
    HKL h (h1.h_ + h2.h_, h1.k_+h2.k_, h1.l_+h2.l_);
    return h;
}

/**
 * difference of two reflections 
 */
HKL operator- (const HKL& h1, const HKL& h2) {
    HKL h (h1.h_-h2.h_, h1.k_-h2.k_, h1.l_-h2.l_);
    return h;
}

/**
 * NB: this calculates R^t * h
 * @param R
 * @param h
 * @return 
 */
HKL operator *(const Int3x3& R, const HKL& h) {
    HKL hkl (0,0,0);
    
    hkl.h_ = R.R_[0] * h.h_ + R.R_[3] * h.k_ + R.R_[6] * h.l_;
    hkl.k_ = R.R_[1] * h.h_ + R.R_[4] * h.k_ + R.R_[7] * h.l_;
    hkl.l_ = R.R_[2] * h.h_ + R.R_[5] * h.k_ + R.R_[8] * h.l_;
    
    return hkl;
    
}

/**
 * output operator
 */
std::ostream& operator<<(std::ostream& outp, const FCFitem& fcfitem) {
    outp << fcfitem.hkl_
            << std::setw(10) << std::setprecision(2) << fcfitem.Imeas_
            << std::setw(8) << std::setprecision(2) << fcfitem.sigImeas_
            << std::setw(10) << std::setprecision(2) << fcfitem.Fcalc_
            << std::setw(7) << std::setprecision(1) << fcfitem.phicalc_;
    return outp;
}

std::ostream& operator<<(std::ostream& outp, const HKL& hkl) {
    outp << std::setw(5) << hkl.h_
            << std::setw(5) << hkl.k_
            << std::setw(5) << hkl.l_;
    return outp;
}