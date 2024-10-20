/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Cube.h
 * Author: tg
 *
 * Created on October 8, 2024, 3:35 PM
 */

#ifndef CUBE_H
#define CUBE_H

#include "cbAtom.h"
#include "Vec3.h"
#include <vector>
#include <string>

/**
 *Cube map format description from 
 * https://paulbourke.net/dataformats/cube/
 * -> Z is the fast changing direction, X is slowest, i.e.
 * idx = iz + Vz_*(iy + Vy_*ix)
 * positioning:
 * xyz = origin + ix*ex_ + iy*ey_ + iz* ez_;
 */
class Cube {
private:
    //! two header lines
    std::string h1_, h2_;
    //! number of atoms
    int numAtoms_;
    //! origin
    Vec3 origin_;
    //! number of voxels in each direction
    int Vx_, Vy_, Vz_;
    //! directions of grid
    Vec3 ex_, ey_, ez_;
    std::vector<cbAtom> cbatoms_;
    std::vector<double> gridvalues_;
    short verbosity_;

    void readMap(const std::string& fname);

public:
    Cube() = default;
    Cube(const std::string& filename, short verbosity_);
    ~Cube() = default;
    
    //! getter functions
    int numAtoms() const { return numAtoms_; }  
    double mapValue(const Vec3&) const;
    size_t gridIndex(int ix, int iy, int iz) const;
    double mapValue(int ix, int iy, int iz) const;
    // return position of atom idx
    Vec3 pos(unsigned short idx) const;

    std::vector<Atom> atoms() const;
    
    // return list of distances to position
    std::vector<double> distances_sq(const Vec3& pos) const;
    // compute the trace of moduli of distances^2
    double deltaTrace(const Cube& cubemap) const;
    
    //! compute pearson coefficient with a second grid
    double CC(const Cube& cube) const;

};

#endif /* CUBE_H */

