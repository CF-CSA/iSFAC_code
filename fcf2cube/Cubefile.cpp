/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Cubefile.cpp
 * Author: tg
 * 
 * Created on May 9, 2025, 5:57 PM
 */

#include "Cubefile.h"
#include "Utils.h"
#include "defines.h"

#include <fstream>
#include <iomanip>

Cubefile::Cubefile(const std::string& filename, const std::array<std::string, 2>& header, short verbosity):
verbosity_(verbosity)
{
    outp_.open(filename);
    if (! outp_.is_open()) {
        if (verbosity_ > 0) {
            std::cout << Utils::error(1) 
                    << "Cannot open file " << filename << " for writing\n"
                    << "    Please check permissions!\n";
            throw myExcepts::FileIO("Cubefile " + filename);
        }
    }
    outp_ << header.front() << '\n'
            << header[1] << '\n';
}

Cubefile::~Cubefile() {
    if (outp_.is_open()) {
        outp_.close();
    }
}

/**
 * write information in Cube format. resfile includes information about coorainte system and
 * the origin of its VdW volume bounding box. 
 * @param resfile
 * @param grid
 * @param unitvecs
 * @return 
 */
int Cubefile::writeCube(const ResFile& resfile, const std::array<int, 3>& grid, 
        const std::array<Vec3, 3>& unitvecs, const std::vector<double> data) {
    const double s = 1.0/Physics::a0;
    outp_ << resfile.atom_list_for_cube();
    /*
    outp_ << std::setw(5) << resfile.num_atoms()
            << std::setw(12) << std::setprecision(6) << s*resfile.vdw_llc().x()
            << std::setw(12) << std::setprecision(6) << s*resfile.vdw_llc().y()
            << std::setw(12) << std::setprecision(6) << s*resfile.vdw_llc().z()
            << '\n';
    outp_ << std::setw(5) << grid[0] 
            << std::setw(12) << std::setprecision(6) << s*unitvecs[0].x()
            << std::setw(12) << std::setprecision(6) << s*unitvecs[0].y()
            << std::setw(12) << std::setprecision(6) << s*unitvecs[0].z()
            << '\n';
    outp_ << std::setw(5) << grid[1] 
            << std::setw(12) << std::setprecision(6) << s*unitvecs[1].x()
            << std::setw(12) << std::setprecision(6) << s*unitvecs[1].y()
            << std::setw(12) << std::setprecision(6) << s*unitvecs[1].z()
            << '\n';
    outp_ << std::setw(5) << grid[2] 
            << std::setw(12) << std::setprecision(6) << s*unitvecs[2].x()
            << std::setw(12) << std::setprecision(6) << s*unitvecs[2].y()
            << std::setw(12) << std::setprecision(6) << s*unitvecs[2].z()
            << '\n';
    
    for (int i = 0; i < resfile.num_atoms(); ++i) {
        outp_ << std::setw(5) << resfile.atomZ(i)
                << std::setprecision(6) << 1.0*resfile.atomZ(i)
                << std::setw(12) << std::setprecision(6) << resfile.atomxyz.x()
                << std::setw(12) << std::setprecision(6) << resfile.atomxyz.y()
                << std::setw(12) << std::setprecision(6) << resfile.atomxyz.z()
                << '\n';
    }
    
    // write out the map, assumed to be in the correct order
    // Nx slow, Ny mid, Nz fast
    for (size_t idxy = 0; idxy < Nx_*Ny_; ++idxy) {
        for (size_t idz = 0; idz < Nz_; ++idz) {
            size_t index = idz + Nz_*idxy;
            outp_ << std::setw(13) << std::setprecision(5) << data[index];
            
        }
        outp_ << '\n';
    }
    */
    outp_.close();
    return  0;
}
