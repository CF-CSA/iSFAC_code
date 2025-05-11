/* 
 * File:   Int3x3.cpp
 * Author: tg
 * 
 * Created on June 28, 2012, 10:37 PM
 */

#include "Int3x3.h"
//#include "HKLops.h"

#include <iomanip>

Int3x3::Int3x3() {
    R_[0] = 0;
    R_[1] = 0;
    R_[2] = 0; 
    R_[3] = 0;
    R_[4] = 0;
    R_[5] = 0;
    R_[6] = 0;
    R_[7] = 0;
    R_[8] = 0;
    T_[0] = 0.0;
    T_[1] = 0.0;
    T_[2] = 0.0;
}

Int3x3::Int3x3(int a11, int a12, int a13, 
    int a21, int a22, int a23, 
    int a31, int a32, int a33,
    double t1, double t2, double t3) {
    R_[0] = a11;
    R_[1] = a12;
    R_[2] = a13; 
    R_[3] = a21;
    R_[4] = a22;
    R_[5] = a23;
    R_[6] = a31;
    R_[7] = a32;
    R_[8] = a33;
    T_[0] = t1;
    T_[1] = t2;
    T_[2] = t3;
}

int Int3x3::det() const {
    int D = 0;
    D = R_[0]*(R_[4] * R_[8] - R_[7]*R_[5]) 
        - R_[3]*(R_[1]*R_[8] - R_[7]*R_[2])
        + R_[6]*(R_[1]*R_[5] - R_[4]*R_[2]);
    return D;
}


std::ostream& operator<< (std::ostream& outp, const Int3x3& R) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            outp << std::setw(3) << R.R_[3*i+j];
        }
        outp << std::fixed
                << std::setw(7) 
                << std::setprecision(2)
                << R.T_[i];
        outp << "\n";
    }
    return outp;
}

bool operator== (const Int3x3 R1, const Int3x3 R2) {
    bool res (true);
    
    for (int i=0; i <9; ++i) {
        res &= (R1.R_[i] == R2.R_[i]);
    }
    return res;
}

bool operator< (const Int3x3 R1, const Int3x3 R2) {
    for (int i = 0; i < 9; ++i) {
        if (R1.R_[i] < R2.R_[i]) return true;
    }
    return false;
}