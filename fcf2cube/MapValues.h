/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MapValues.h
 * Author: tg
 *
 * Created on September 5, 2023, 4:25 PM
 */

#ifndef MAPVALUES_H
#define MAPVALUES_H

#include "Vec3.h"
#include <vector>

class MapValues {
private:
    std::vector<double> map_;
    const int gridx_, gridy_, gridz_;
    const double dx_, dy_, dz_;
    
    int verbosity_;
    
    //! ensure @c x is inside unit cell
    double in_unitcell(const double& x) const;
    Vec3 in_unitcell(const Vec3& x) const;
    
    //! compute the index into map from grid coordinates
    int idx(int j, int k, int l) const;
    
    //! calculate coordinates of lower bound grid point
    void gridpoint(const Vec3& xyz, int& j, int& k, int& l) const;
public:
    MapValues(const std::vector<double>& map, int gridx, int gridy, int gridz, int verbosity);
    MapValues(const MapValues& orig) = default;
    ~MapValues() = default;
    
    //! return interpolated value at position @c coords
    double mapvalue (const Vec3& XYZ) const;

};

#endif /* MAPVALUES_H */

