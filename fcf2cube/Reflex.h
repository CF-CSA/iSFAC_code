/* 
 * File:   Reflex.h
 * Author: tg
 *
 * Created on July 3, 2012, 4:30 PM
 */

#ifndef REFLEX_H
#define	REFLEX_H

#include "HKL.h"

/**
 * \class Reflex
 * \brief Storage and handling of reflection
 * 
 * The 'Reflex' in this project only contains the material used
 * for the FFT, i.e., all modifications like scaling or (weighted) 
 * differences must be done prior to creating a Reflex
 */
class Reflex {
private:
    HKL hkl_;
    double F_;
    double phi_;
public:
    Reflex();
    Reflex(const HKL& hkl, double F, double phi);
    ~Reflex();
    
    double F() const { return F_; }
    double phi() const { return phi_; }
    HKL hkl() const { return hkl_; }
    
    void stdsetting (const std::vector<Int3x3>& symops);
    double phaseshift (const Int3x3& symop) const;
    double dstarsq (const Vec3& astar, const Vec3& bstar, const Vec3& cstar) const;
    
    friend bool operator <  (const Reflex& R1, const Reflex& R2);
    friend bool operator == (const Reflex& R1, const Reflex& R2);
    friend std::ostream& operator<< (std::ostream& outp, const Reflex& R);
};

#endif	/* REFLEX_H */

