/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Mat33.h
 * Author: tg
 *
 * Created on July 25, 2023, 9:59 PM
 */

#ifndef MAT33_H
#define MAT33_H

#include <iostream>
class XYZ;

class Mat33 {
private:
    double matrix_[9];
public:
    Mat33();
    Mat33(const double& m00, const double& m01, const double& m02,
            const double& m10, const double& m11, const double& m12,
            const double& m20, const double& m21, const double& m22);
    ~Mat33();

    // write access
    double& operator() (unsigned short idx1, unsigned short idx2) &{
        return matrix_ [3*idx1+idx2];
    }
    
    // read access
    const double& operator() (unsigned short idx1, unsigned short idx2) const&{
        return matrix_[3*idx1+idx2];
    }
    friend Mat33 operator*(const Mat33& m1, const Mat33& m2);
    friend XYZ operator* (const Mat33& m, const XYZ& xyz);
    friend std::ostream& operator<<(std::ostream& outp, const Mat33& m);
};

#endif /* MAT33_H */

