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
#include <complex>
#include <algorithm>
#include <numeric>
#include <limits>
#include <gsl/gsl_statistics_double.h>

/**
 * Initiates ex, ey, ez with unit cell vectors. Need to be multiplied with 
 * delta of grid fragments
 * @param resfile
 * @param mapvals
 * @param header
 * @param verbosity
 */
Cubefile::Cubefile(const ResFile& resfile, const MapValues& mapvals, double margin, short verbosity) :
resfile_(resfile),
mapvals_(mapvals),
margin_(margin),
ex_(resfile.A()),
ey_(resfile.B()),
ez_(resfile.C()),
verbosity_(verbosity) {
    Nx_ = mapvals_.gridx();
    Ny_ = mapvals_.gridy();
    Nz_ = mapvals_.gridz();
    std::tie(llc_, urc_) = resfile_.bbox_frac();

    if (verbosity > 1) {
        const double s = 1. / Physics::a0;
        std::cout << Utils::prompt(1) << "Setting up Cubefile. Initial unit cell vectors:\n"
                << Utils::prompt(1) << " A[A] = " << ex_ << ", A[a0] = " << s * ex_ << '\n'
                << Utils::prompt(1) << " B[A] = " << ey_ << ", B[a0] = " << s * ey_ << '\n'
                << Utils::prompt(1) << " C[A] = " << ez_ << ", C[a0] = " << s * ez_ << '\n'
                << Utils::prompt(1) << "BBox: [" << llc_ << ", " << urc_ << "]\n";
    }

    // expand bbox by margin
    Vec3 diagonal = urc_ - llc_;
    llc_ = llc_ - margin*diagonal;
    urc_ = urc_ + margin*diagonal;
    if (verbosity_ > 1) {
        std::cout << Utils::prompt(1) << "expanding bbox diagonal by " <<
                margin * diagonal << " to [" << llc_ << ", " << urc_ << "]\n";
    }
    makemap();

}

void Cubefile::prepAtomlist() {
    resfile_.skip_close_atoms();
}

/**
 * write information in Cube format. resfile includes information about coorainte system and
 * the origin of its VdW volume bounding box. 
 * All conversion from A to Bohr should take place here!
 * @param resfile
 * @param grid
 * @param unitvecs
 * @return 
 */
int Cubefile::writeCube(const std::string filename,
        const std::array<std::string, 2>& header) const {
    std::ofstream outp_(filename);
    if (!outp_.is_open()) {
        if (verbosity_ > 0) {
            std::cout << Utils::error(1)
                    << "Cannot open file " << filename << " for writing\n"
                    << "    Please check permissions!\n";
            throw myExcepts::FileIO("Cubefile " + filename);
        }
    }
    outp_ << header.front() << '\n'
            << header[1] << '\n';
    int Natoms(0);
    std::string atomlist = resfile_.atom_list_for_cube(Natoms);

    const double s = 1. / Physics::a0;
    Vec3 origin(llc_.x() * resfile_.A() + llc_.y() * resfile_.B() + llc_.z() * resfile_.C());

    if (verbosity_ > 2) {
        std::cout << Utils::prompt(2) << "Natoms, s, and origin: "
                << Natoms << ' ' << s << ' ' << origin << '\n'
                << " converted to Bohr with scale factor " << s * origin << '\n'
                << Utils::prompt(2) << "Grid covers in Bohr:\n"
                << Utils::prompt(2) << s*origin << " to " 
                << s*origin + Nx_*s*ex_ + Ny_*s*ey_ + Nz_*s*ez_ << '\n';
    }
    outp_ << std::setw(5) << Natoms
            << std::setw(12) << std::setprecision(6) << std::fixed << s * origin.x()
            << std::setw(12) << std::setprecision(6) << std::fixed << s * origin.y()
            << std::setw(12) << std::setprecision(6) << std::fixed << s * origin.z()
            << " 1\n";

    outp_ << std::setw(5) << Nx_
            << std::setw(12) << std::setprecision(6) << s * ex_.x()
            << std::setw(12) << std::setprecision(6) << s * ex_.y()
            << std::setw(12) << std::setprecision(6) << s * ex_.z()
            << '\n';
    outp_ << std::setw(5) << Ny_
            << std::setw(12) << std::setprecision(6) << s * ey_.x()
            << std::setw(12) << std::setprecision(6) << s * ey_.y()
            << std::setw(12) << std::setprecision(6) << s * ey_.z()
            << '\n';
    outp_ << std::setw(5) << Nz_
            << std::setw(12) << std::setprecision(6) << s * ez_.x()
            << std::setw(12) << std::setprecision(6) << s * ez_.y()
            << std::setw(12) << std::setprecision(6) << s * ez_.z()
            << '\n';

    outp_ << atomlist;

    // write out the map, assumed to be in the correct order
    // Nx slow, Ny mid, Nz fast
    for (size_t idxy = 0; idxy < Nx_ * Ny_; ++idxy) {
        for (size_t idz = 0; idz < Nz_; ++idz) {
            size_t index = idz + Nz_*idxy;
            outp_ << std::scientific << std::setw(13) << std::setprecision(5) << map_[index];

        }
        outp_ << '\n';
    }
    outp_.close();
    return 0;
}

/**
 * fills @c map_ with data from FFT map. Updates grid vectors
 * ex_, ey_, ez_ left in Angstrom
 * @param llc
 * @param urc
 */
void Cubefile::makemap() {
    // this function works in fractional coordinates
    const Vec3 diag = urc_ - llc_;
    const double deltaX = diag.x() / Nx_;
    const double deltaY = diag.y() / Ny_;
    const double deltaZ = diag.z() / Nz_;

    if (deltaX <= 0.0 || deltaY <= 0.0 || deltaZ <= 0.0) {
        if (verbosity_ > 0) {
            std::cout << Utils::error(1) << "Error making map: upper corner is not a positive diagional\n";
            throw myExcepts::Programming("make sure llc < urc in x,y,z");
        }
    }

    map_.resize(Nx_ * Ny_ * Nz_, 0.0);

    // cube maps: X is slow, Z is fast
    for (size_t ix = 0; ix < Nx_; ++ix) {
        const double x = llc_.x() + ix*deltaX;
        for (size_t iy = 0; iy < Ny_; ++iy) {
            const double y = llc_.y() + iy*deltaY;
            for (size_t iz = 0; iz < Nz_; ++iz) {
                const double z = llc_.z() + iz*deltaZ;
                Vec3 point (x, y, z);
                double val = mapvals_.mapvalue(point);
                size_t idx = iz + Nz_ * (iy + Ny_ * ix);
                map_[idx] = val;
            }
        }
    }
    // update setps in X, Y, Z for grid
    ex_ = deltaX*ex_;
    ey_ = deltaY*ey_;
    ez_ = deltaZ*ez_;
}

void Cubefile::printState() const {
    std::vector<double> map(map_);
    std::sort(map.begin(),map.end());
    double sum = std::accumulate(map.begin(), map.end(), 0.0);
    const double mu = sum/map.size();
    double var(0.0);
    for (auto x: map) {
        var += (x-mu)*(x-mu);
    }
    var = var/map.size();
    
    std::cout << Utils::prompt(verbosity_) << " mu  = " << mu << '\n'
            << Utils::prompt(verbosity_)   << " sd  = " << std::sqrt(var) << '\n'
            << Utils::prompt(verbosity_)   << " min = " << map.front() << '\n'
            << Utils::prompt(verbosity_)   << " max = " << map.back() << '\n';
}