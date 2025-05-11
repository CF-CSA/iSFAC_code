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

#include <fstream>
#include <iomanip>

Cubefile::Cubefile() {
}

int Cubefile::writeCube(const std::string& filename, const std::array<std::string,2>& header) const {
    std::ofstream outp(filename);
    if (!outp.is_open()) {
        std::cout << Utils::prompt(1) << "Error opening file " << filename 
                << std::endl;
        throw myExcepts::FileIO("Writing "+filename);
    }
    outp << header[1] << '\n'
            << header[2] << '\n';
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
