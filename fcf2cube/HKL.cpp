/* 
 * File:   HKL.cpp
 * Author: tg
 * 
 * Created on June 28, 2012, 10:32 PM
 */

#include "HKL.h"
#include "myExceptions.h"
double HKL::dstarsq(const Vec3& astar, const Vec3& bstar, const Vec3& cstar) const {
    Vec3 p = 1.0*h_*astar + 1.0*k_*bstar + 1.0*l_*cstar;
    double invdsq (p.lengthsq() );
    return invdsq;
}

int HKL::operator ()(unsigned char idx) const {
    switch (idx) {
        case 0: return h_;
            break;
        case 1: return k_;
        break;
        case 2: return l_;
        break;
        default:
            throw myExcepts::Programming("Index operator for HKL out of range.");
    }

}