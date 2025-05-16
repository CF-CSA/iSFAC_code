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
#include "MapValues.h"

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
    ResFile resfile_;
    MapValues mapvals_;
    int Nx_, Ny_, Nz_;
    //! lower left corner of box is origin of map
    Vec3 llc_, urc_;
    Vec3 ex_, ey_, ez_;
    float margin_;
    short verbosity_;
    
    void makemap();
    
public:
    Cubefile() = delete;
    Cubefile(const ResFile& resfile, const MapValues& mapvals, float margin, short verbosity);
    ~Cubefile() = default;
    
    //! write info to cube file, including two header lines
    int writeCube(const std::string filename, 
        const std::array<std::string, 2>& header) const;

};

#endif /* CUBEFILE_H */

