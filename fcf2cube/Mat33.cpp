/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Mat33.cpp
 * Author: tg
 * 
 * Created on July 25, 2023, 9:59 PM
 */

#include <iomanip>

#include "Mat33.h"
#include "Vec3.h"

/**
 * create a zero filled matrix
 */
Mat33::Mat33() {
    for (int i = 0; i < 9; ++i) {
        matrix_[i] = 0.0;
    }
}

Mat33::Mat33(bool unit) : Mat33() {
    if (unit) {
        matrix_[0] = matrix_[4] = matrix_[8] = 1.0;
    }
}

Mat33::Mat33(const double& m00, const double& m01, const double& m02,
        const double& m10, const double& m11, const double& m12,
        const double& m20, const double& m21, const double& m22) {
    matrix_[0] = m00;
    matrix_[1] = m01;
    matrix_[2] = m02;
    matrix_[3] = m10;
    matrix_[4] = m11;
    matrix_[5] = m12;
    matrix_[6] = m20;
    matrix_[7] = m21;
    matrix_[8] = m22;
}

Mat33::~Mat33() {
}

Mat33 operator*(const Mat33& m1, const Mat33& m2) {
    Mat33 M;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            for (int i = 0; i < 3; ++i) {
                M(r, c) += m1(r, i) * m2(i, c);
            }
        }
    }
    return M;
}

Vec3 operator*(const Mat33& m, const Vec3& xyz) {
    double x[3];
    for (int r = 0; r < 3; r++) {
        x[r] = m(r, 0) * xyz.x() + m(r, 1) * xyz.y() + m(r, 2) * xyz.z();
    }
    return Vec3(x[0], x[1], x[2]);
}

std::ostream& operator<<(std::ostream& outp, const Mat33& m) {
    outp << std::setw(8) << std::setprecision(2)
            << m(0, 0) << ' ' << m(0, 1) << ' ' << m(0, 2) << '\n'
            << "    " << m(1, 0) << ' ' << m(1, 1) << ' ' << m(1, 2) << '\n'
            << "    " << m(2, 0) << ' ' << m(2, 1) << ' ' << m(2, 2);
    return outp;
}

double Mat33::determinant() const {
    double d = 0;
    d  = matrix_[0] * (matrix_[4]*matrix_[8] - matrix_[5]*matrix_[7]);
    d -= matrix_[1] * (matrix_[3]*matrix_[8] - matrix_[5]*matrix_[6]);
    d += matrix_[2] * (matrix_[3]*matrix_[7] - matrix_[4]*matrix_[6]);
    return d;
}