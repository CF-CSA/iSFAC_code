/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MapValues.cpp
 * Author: tg
 * 
 * Created on September 5, 2023, 4:25 PM
 */

#include <cmath>

#include "MapValues.h"

MapValues::MapValues(const std::vector<double>& map, int gridx, int gridy, int gridz, int verbosity) :
map_(map), 
        gridx_(gridx), gridy_(gridy), gridz_(gridz), 
        dx_(1./gridx_), dy_(1./gridy_), dz_(1./gridz_), 
        verbosity_(verbosity){

}

int MapValues::idx(int j, int k, int l) const {
    int m = j + gridx_* (k+ gridy_*l);
    return m;
}

/**
 * ensue 0<= x < 1, i.e. x inside unit cell
 * @param x
 * @return 
 */
double MapValues::in_unitcell(const double& x) const {
    float X = std::fmod(x, 1.0);
    if (X<0) X+= 1.0;
    return X;
}

Vec3 MapValues::in_unitcell(const Vec3& x) const 
{
    const double xu = in_unitcell(x.x());
    const double yu = in_unitcell(x.y());
    const double zu = in_unitcell(x.z());
    Vec3 X = Vec3(xu, yu, zu);
    return X;
}

/**
 * compute the lower left corner grid point for xyz
 * @param xyz
 * @param gx resulting x- index for grid point
 * @param gy resulting y- index for grid point
 * @param gz resulting z- index for grid point
 */
void MapValues::gridpoint(const Vec3& xyz, int& gx, int& gy, int& gz) const{
    double x = in_unitcell(xyz.x());
    double y = in_unitcell(xyz.y());
    double z = in_unitcell(xyz.z());
    gx = int (gridx_ * x);
    gy = int (gridy_ * y);
    gz = int (gridz_ * z);
}
/**
 * computes the value of the map at @c coord as distance-weighted 
 * average of the surrounding corner points around coord
 * Trilinear interpolation as described in "A Survey of Voxel Interpolation Methods and an Evaluation of Their
 * Impact on Volumetric Map-Based Visual Odometry", D. R. Canelhas, T. Stoyanov, A. J. Lilienthal, Sweden
 * http://iliad-project.eu/wp-content/uploads/2018/03/dlcs_icra_2018.pdf
 * @param coords
 * @return 
 */
double MapValues::mapvalue(const Vec3& XYZ) const {
    // offsets for 8 corner points surrounding X
    const double a[] = {0, 0, 0};
    const double b[] = {1, 0, 0};
    const double c[] = {0, 1, 0};
    const double d[] = {0, 0, 1};
    const double e[] = {0, 1, 1};
    const double f[] = {1, 0, 1};
    const double g[] = {1, 1, 0};
    const double h[] = {1, 1, 1};
    // get grid coordinate of a, lower bound for @c XYZ
    int gx, gy, gz;
    Vec3 ucXYZ = Vec3(in_unitcell(XYZ));
    gridpoint (ucXYZ, gx, gy, gz);
    // fractional distance from A (gx, gy, gz)
    const double dx = ucXYZ.x() - gx*dx_;
    const double dy = ucXYZ.y() - gy*dy_;
    const double dz = ucXYZ.z() - gz*dz_;
    
    // there is most likely much room for optimisation here.
    
    // indices and map values for corner points
    const int idx_a = idx(gx+a[0], gy+a[1], gz+a[2]);
    const double map_a = map_[idx_a];
    const int idx_b = idx(gx+b[0], gy+b[1], gz+b[2]);
    const double map_b = map_[idx_b];
    const int idx_c = idx(gx+c[0], gy+c[1], gz+c[2]);
    const double map_c = map_[idx_c];
    const int idx_d = idx(gx+d[0], gy+d[1], gz+d[2]);
    const double map_d = map_[idx_d];
    const int idx_e = idx(gx+e[0], gy+e[1], gz+e[2]);
    const double map_e = map_[idx_e];
    const int idx_f = idx(gx+f[0], gy+f[1], gz+f[2]);
    const double map_f = map_[idx_f];
    const int idx_g = idx(gx+g[0], gy+g[1], gz+g[2]);
    const double map_g = map_[idx_g];
    const int idx_h = idx(gx+h[0], gy+h[1], gz+h[2]);
    const double map_h = map_[idx_h];
    
    const double map_p = map_a * (1 - dx) + map_b * dx;
    const double map_q = map_c * (1 - dx) + map_g * dx;
    const double map_r = map_d * (1 - dx) + map_f * dx;
    const double map_s = map_e * (1 - dx) + map_h * dx;
    
    const double map_t = map_p * (1 - dy) + map_q * dy;
    const double map_u = map_r * (1 - dy) + map_s * dy;
    
    const double valxyz = map_t*(1 - dz) + map_u*dz;
    
    return valxyz;
}