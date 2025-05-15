/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Cubefile.h
 * Author: tg
 *
 * Created on May 9, 2025, 5:57 PM
 */

#ifndef CUBEFILE_H
#define CUBEFILE_H

#include "Atom.h"
#include "Vec3.h"
#include "myExceptions.h"
#include "ResFile.h"

#include <vector>
#include <string>
#include <array>
#include <fstream>
/**
 * Create a cube-file from SHELX project with
 * FCF file, RES file..
 */
class Cubefile {
private:
    std::vector<Atom> atoms_;
    std::vector<double> map_;
    int Nx_, Ny_, Nz_;
    Vec3 origin_;
    Vec3 ex_, ey_, ez_;
    std::ofstream outp_;
    short verbosity_;
    
    int cubeHeader(const ResFile& resfile);
    
public:
    Cubefile() = delete;
    Cubefile(const std::string& filename, const std::array<std::string, 2>& header, short verbosity);
    ~Cubefile();
    
    int writeCube(const ResFile& resfile, const std::array<int, 3>& grid, 
        const std::array<Vec3, 3>& unitvecs, const std::vector<double> data);

};

#endif /* CUBEFILE_H */

