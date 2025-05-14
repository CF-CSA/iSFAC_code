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

int Cubefile::writeCube(const ResFile& resfile, const std::array<int, 3>& grid, const std::array<Vec3, 3>& unitvecs) const {
    const double s = 1.0/Physics::a0;
    outp_ << std::setw(5) << resfile.num_atoms()
            << std::setw(12) << std::setprecision(6) << s*resfile.orig().x()
            << std::setw(12) << std::setprecision(6) << s*resfile.orig().y()
            << std::setw(12) << std::setprecision(6) << s*resfile.orig().z()
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
            
            
    outp << std::setw(5) << atoms_.size()
            << std::setw(12) << std::setprecision(6) << origin_.x()
            << std::setw(12) << std::setprecision(6) << origin_.y()
            << std::setw(12) << std::setprecision(6) << origin_.z()
            << '\n';
    outp << std::setw(5) << Nx_ 
            << std::setw(12) << std::setprecision(6) << ex_.x()
            << std::setw(12) << std::setprecision(6) << ex_.y()
            << std::setw(12) << std::setprecision(6) << ex_.z()
            << '\n';
    outp << std::setw(5) << Ny_ 
            << std::setw(12) << std::setprecision(6) << ey_.x()
            << std::setw(12) << std::setprecision(6) << ey_.y()
            << std::setw(12) << std::setprecision(6) << ey_.z()
            << '\n';
    outp << std::setw(5) << Nz_ 
            << std::setw(12) << std::setprecision(6) << ez_.x()
            << std::setw(12) << std::setprecision(6) << ez_.y()
            << std::setw(12) << std::setprecision(6) << ez_.z()
            << '\n';
    for (auto a: atoms_) {
        outp << std::setw(5) << a.Z()
            << std::setw(12) << std::setprecision(6) << a.x()
            << std::setw(12) << std::setprecision(6) << a.y()
            << std::setw(12) << std::setprecision(6) << a.z()
            << '\n';
    }
    
    // write out the map, assumed to be in the correct order
    // Nx slow, Ny mid, Nz fast
    for (size_t idxy = 0; idxy < Nx_*Ny_; ++idxy) {
        for (size_t idz = 0; idz < Nz_; ++idz) {
            size_t index = idz + Nz_*idxy;
            outp << std::setw(13) << std::setprecision(5) << map_[index];
            
        }
        outp << '\n';
    }
 
    outp.close();
    return  0;
}
