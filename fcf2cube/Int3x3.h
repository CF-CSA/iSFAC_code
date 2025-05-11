/* 
 * File:   Int3x3.h
 * Author: tg
 *
 * Created on June 28, 2012, 10:37 PM
 */

#ifndef INT3X3_H
#define	INT3X3_H

class HKL;
#include <iostream>

/**
 * 3x3 matrix with integer coefficients: symmetry operation for hkl's
 * also stores translational part in order to correct for phase transformation
 */
class Int3x3 {
private:
    int R_[9];
    double T_[3];
public:
    Int3x3();
    Int3x3(int a11, int a12, int a13, 
    int a21, int a22, int a23, 
    int a31, int a32, int a33,
    double t1, double t2, double t3);
    ~Int3x3(){}
    friend HKL operator*(const Int3x3& R, const HKL& h);
    /* read and write matrix entries */
    int operator() (int i, int j) const {
        return R_[3*i+j];
    }
    int& operator() (int i, int j) {
        return R_[3*i+j];
    }
    /* read and write vector entries */
    double operator() (int i) const {
        return T_[i];
    }
    double& operator() (int i) {
        return T_[i];
    }
    
    //! return determinant of R_
    int det() const;

    friend std::ostream& operator<< (std::ostream& outp, const Int3x3& R);
    friend bool operator== (const Int3x3 R1, const Int3x3 R2);
    friend bool operator<  (const Int3x3 R1, const Int3x3 R2);

};

#endif	/* INT3X3_H */

