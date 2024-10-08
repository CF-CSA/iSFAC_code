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
class Cube {
private:
    //! two header lines
    std::string h1_, h2_;
    //! origin
    int numAtoms_;
    //! number of atoms
    float orgx_, orgy_, orgz_;
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
    

};

#endif /* CUBE_H */

