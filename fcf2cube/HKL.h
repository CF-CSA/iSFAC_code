/* 
 * File:   HKL.h
 * Author: tg
 *
 * Created on June 28, 2012, 10:32 PM
 */

#ifndef HKL_H
#define	HKL_H

#include <vector>
#include <ostream>
#include "Int3x3.h"
#include "Vec3.h"

class HKL {
private:
    int h_;
    int k_;
    int l_;
public:
    
    HKL(){}
    HKL(int h, int k, int l): h_(h), k_(k), l_(l) {}
    
    int h() const { return h_; }
    int k() const { return k_; }
    int l() const { return l_; }
    
    int operator()(unsigned char idx) const;
    /* inverse square resolution */
    double dstarsq(const Vec3& astar, 
        const Vec3& bstar, const Vec3& cstar) const;
    
    /* defined in 'HKLops.cpp' */
    friend bool operator< (const HKL& h1, const HKL& h2);
    friend bool operator== (const HKL& h1, const HKL& h2);
    friend bool operator> (const HKL& h1, const HKL& h2);
    
    friend HKL operator+(const HKL& h1, const HKL& h2);
    friend HKL operator-(const HKL& h1, const HKL& h2);
    friend HKL operator*(const Int3x3& R, const HKL& h);
    
    // output operator
    friend std::ostream& operator<<(std::ostream& outp, const HKL& hkl);
};


#endif	/* HKL_H */

