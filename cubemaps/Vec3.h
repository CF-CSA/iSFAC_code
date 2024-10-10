/* 
 * File:   Vec3.h
 * Author: tg
 *
 * Created on July 24, 2012, 11:12 PM
 */

#ifndef VEC3_H
#define	VEC3_H

class Vec3 {
private:
    double x_, y_, z_;
public:
    Vec3() {}
    Vec3(double x, double y, double z): x_(x), y_(y), z_(z) {}
    ~Vec3() {}
    
    double lengthsq() const;
    
    double x() const throw() { return x_; }
    double y() const throw() { return y_; }
    double z() const throw() { return z_; }
    
    friend double operator*(const Vec3& v1, const Vec3& v2);
    friend Vec3 operator* (const double s, const Vec3& v);
    friend Vec3 operator+ (const Vec3& v1, const Vec3& v2);
    friend Vec3 operator- (const Vec3& v1, const Vec3& v2);
    friend Vec3 cross( const Vec3& v1, const Vec3& v2);

};

double operator*(const Vec3& v1, const Vec3& v2);
Vec3 operator* (const double s, const Vec3& v);
Vec3 operator+ (const Vec3& v1, const Vec3& v2);
Vec3 cross( const Vec3& v1, const Vec3& v2);

#endif	/* VEC3_H */

