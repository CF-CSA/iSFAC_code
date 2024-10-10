/* 
 * File:   Vec3.cpp
 * Author: tg
 * 
 * Created on July 24, 2012, 11:12 PM
 */

#include "Vec3.h"

double Vec3::lengthsq() const {
    return *this * *this;
}

double operator*(const Vec3& v1, const Vec3& v2) {
    double s (0.0);
    s += v1.x_ * v2.x_;
    s += v1.y_ * v2.y_;
    s += v1.z_ * v2.z_;
    
    return s;
}

Vec3 operator*(const double s, const Vec3& v) {
    Vec3 r(v);
    r.x_ *= s;
    r.y_ *= s;
    r.z_ *= s;
    
    return r;
}

Vec3 cross (const Vec3& v1, const Vec3& v2) {
    Vec3 v;
    v.x_ = v1.y_ * v2.z_ - v1.z_ * v2.y_;
    v.y_ = v1.z_ * v2.x_ - v1.x_ * v2.z_;
    v.z_ = v1.x_ * v2.y_ - v1.y_ * v2.x_;
    
    return v;
}

Vec3 operator+ (const Vec3& v1, const Vec3& v2) {
    Vec3 v(v1);
    v.x_ += v2.x_;
    v.y_ += v2.y_;
    v.z_ += v2.z_;
    
    return v;
}

Vec3 operator- (const Vec3& v1, const Vec3& v2) {
    Vec3 v(v1);
    v.x_ -= v2.x_;
    v.y_ -= v2.y_;
    v.z_ -= v2.z_;
    
    return v;
}
