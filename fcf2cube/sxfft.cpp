/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   sxfft.cpp
 * Author: tg
 * 
 * Created on July 29, 2023, 9:40 PM
 */

#include <cmath>
#include <algorithm>
#include <kissfft/kiss_fftnd.h>

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <iomanip>

#include "sxfft.h"
#include "myExceptions.h"

/**
 * magicFC are multiples of 2 3 and 5 
 * FRom sxfft.f: "The first line of the map file gives Nc, Nx, Ny, Nz, p, q and s,
 * followed on subsequent lines by density values scaled so that they can
 * be represented'/' by integers N in the range 0 to 9999. The
 * density in e/A^3 can be calculated'/' by rho=p*N+q, s is the
 * square root of the map variance in e/A^3"
 * @param fcfdata
 * @param fcfinfo
 * @param maptype
 * @param verbosity
 */
sxfft::sxfft(const std::vector<FCFitem>& fcfdata, const FCFInfo& fcfinfo, int maptype, int verbosity) :
fcfdata_(fcfdata),
fcfinfo_(fcfinfo),
magicFC_({2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45,
    48, 50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128, 135, 144, 150,
    160, 162, 180, 192, 200, 216, 225, 240, 243, 250, 256, 270, 288, 300, 320, 324,
    360, 375, 384, 400, 405, 432, 450, 480, 486, 500, 512, 540, 576, 600, 625, 640,
    648, 675, 720, 729, 750, 768, 800, 810, 864, 900, 960, 972, 1000, 1024, 1080,
    1125, 1152, 1200, 1215, 1250, 1280, 1296, 1350, 1440, 1458, 1500, 1536, 1600,
    1620, 1728, 1800, 1875, 1920, 1944, 2000, 2025, 2048, 2160, 2187, 2250, 2304,
    2400, 2430, 2500, 2560, 2592, 2700, 2880, 2916, 3000, 3072, 3125, 3200, 3240,
    3375, 3456, 3600, 3645, 3750, 3840, 3888, 4000, 4050, 4096, 4320, 4374, 4500,
    4608, 4800, 4860, 5000}),
maptype_(maptype),
verbosity_(verbosity),
centrosymmetric_(false) {
    symops_ = fcfinfo_.symops_no_inv();
    if (verbosity_ > 0) {
        std::cout << "---> Setting to standard indices (h>=0, etc)\n"
                << "---> Sorting reflections, and merging identical indices\n";
    }
    standard_hkl();
    std::sort(fcfdata_.begin(), fcfdata_.end());
    merge_data();

}

/**
 * from sxfft.f, with aid from fourxle.cpp (SHELXLE); prepares and carries out
 * FFT for reflection data
 */
void sxfft::fft(const double& weakWeight, const double& gridresol) {

    // find the largest miller indices, including symmetry equivalents
    int mh = 0, mk = 0, ml = 0;
    for (size_t n = 0; n < fcfdata_.size(); n++) {
        double u = fcfdata_[n].hkl().h();
        double v = fcfdata_[n].hkl().k();
        double w = fcfdata_[n].hkl().l();

        for (auto sym = symops_.begin(); sym != symops_.end(); sym++) {
            int a, b, c;
            a = abs((int) (u * (*sym)(0, 0) + v * (*sym)(1, 0) + w * (*sym)(2, 0)));
            b = abs((int) (u * (*sym)(0, 1) + v * (*sym)(1, 1) + w * (*sym)(2, 1)));
            c = abs((int) (u * (*sym)(0, 2) + v * (*sym)(2, 1) + w * (*sym)(2, 2)));
            mh = (mh < a) ? a : mh;
            mk = (mk < b) ? b : mk;
            ml = (ml < c) ? c : ml;
        }
    }
    if (verbosity_ > 2) {
        std::cout << "---> maximal miller indices for fourier grid:\n"
                << "   " << mh << ", " << mk << ", " << ml << '\n';
    }
    // define the FFT grid in 3D as n1, n2, n3
    // resgrid_ increases the resolution of the output grid (default: 5.0)
    grid_n1_ = magicTop(int(gridresol * mh + .5));
    grid_n2_ = magicTop(int(gridresol * mk + .5));
    grid_n3_ = magicTop(int(gridresol * ml + .5));
    if (verbosity_ > 1) {
        std::cout << "---> Preparation of map grid. \n"
                << "     Original number of gridpoints:         "
                << mh << " x " << mk << " x " << ml << "\n"
                << "     multiplied with " << gridresol << "\n"
                << "     final ``magic'' number of grid points: "
                << grid_n1_ << " x " << grid_n2_ << " x " << grid_n3_ << "\n";
    }

    const int n4 = grid_n2_*grid_n1_;
    const int n5 = grid_n3_*n4;

    // define the grid of the map
    const double DX = 1.0 / grid_n1_;
    const double DY = 1.0 / grid_n2_;
    const double DZ = 1.0 / grid_n3_;

    if (verbosity_ > 1) {
        std::cout << "---> grid steps in x, y, z:"
                << "   " << std::fixed << std::setprecision(5) << DX << ", " << DY << ", " << DZ << '\n';
    }

    /**
     * sxfft.f here sets up a Hermitian 3D coefficients for FFT. fourxle.cpp
     * does not do this - maybe not required for kissfft
     */

    // prepare for FFT (kissfft)

    int nbytes;
    kiss_fft_cpx* B = (kiss_fft_cpx*) KISS_FFT_MALLOC(nbytes = (sizeof (kiss_fft_cpx) * n5));
    if (B == NULL) {
        throw std::runtime_error("Cannot allocate memory for KISS FFT");
    }
    // zero array for kiss fft
    for (int i = 0; i < n5; i++) {
        B[i].r = 0;
        B[i].i = 0;
    }
    // Set up data in B (complex Hermitian) for FFT with KISSFFT
    const double c15 = fcfinfo_.fftscale();
    for (int i = 0; i < fcfdata_.size(); i++) {
        double u, v, w;
        u = fcfdata_[i].hkl().h();
        v = fcfdata_[i].hkl().k();
        w = fcfdata_[i].hkl().l();
        double s = 0, t = 0, q, p;
        // count number of symmetry equivalents
        for (auto sym = symops_.begin(); sym != symops_.end(); ++sym) {

            int j = (int) (u * (*sym)(0, 0) + v * (*sym)(1, 0) + w * (*sym)(2, 0));
            int k = (int) (u * (*sym)(0, 1) + v * (*sym)(1, 1) + w * (*sym)(2, 1));
            int l = (int) (u * (*sym)(0, 2) + v * (*sym)(1, 2) + w * (*sym)(2, 2));
            // (j,k,l) == fcfdata_[i].hkl() ??
            // weighting for identical reflections
            if (HKL(j, k, l) == fcfdata_[i].hkl()) {
                s += 1.0;
            }
            // (j,k,l) == -fcfdata_[i].hkl() ??
            // t: weighting for Bijvoet pair
            if (HKL(-j, -k, -l) == fcfdata_[i].hkl()) {
                t += 1.0;
            }
        }

        if (i == 0) {//printf("v%f s%f t%f\n",C[14],s,t);
            s = 1;
            t = 0; //f000 
        }
        // typ = 0 is difference map, otherwise MFo - (M-1)Fc map
        if (maptype_ == 0) {
            s = (fcfdata_[i].Imeas() - fcfdata_[i].Fcalc()) / (c15 * (s + t));
        } else {
            s = (maptype_ * fcfdata_[i].Imeas() - (maptype_ - 1) * fcfdata_[0].Fcalc())
                    / (c15 * (s + t));
        }
        if (fcfdata_[i].Fcalc() > 1.E-6) {
            s = s / (1. + weakWeight * std::pow(fcfdata_[i].sigImeas() / fcfdata_[i].Fcalc(), 4));
        }
        for (auto sym = symops_.begin(); sym != symops_.end(); ++sym) {
            int j, k, l, m;
            j = (int) (u * (*sym)(0, 0) + v * (*sym)(1, 0) + w * (*sym)(2, 0));
            k = (int) (u * (*sym)(0, 1) + v * (*sym)(1, 1) + w * (*sym)(2, 1));
            l = (int) (u * (*sym)(0, 2) + v * (*sym)(1, 2) + w * (*sym)(2, 2));
            //          q=(-2*M_PI*(u*sy[9][n]+v*sy[10][n]+w*sy[11][n]))-M_PI*(j*DX+k*DY+l*DZ);
            // q seems to fit alright, debugged
            q = (fcfdata_[i].phicalc() - 2 * M_PI * (u * (*sym)(0) + v * (*sym)(1) + w * (*sym)(2))) - M_PI * (j * DX + k * DY + l * DZ);
            j = (999 * grid_n1_ + j) % grid_n1_;
            k = (999 * grid_n2_ + k) % grid_n2_;
            l = (999 * grid_n3_ + l) % grid_n3_;
            // compared with sxfft.f: start counting at 0, and half indx, because B is complex
            m = j + grid_n1_ * (k + grid_n2_ * l);
            p = s * std::cos(q);
            B[m].r = p;
            q = s * std::sin(q);
            B[m].i = q;
            j *= -1;
            if (j < 0)j += grid_n1_;
            k *= -1;
            if (k < 0)k += grid_n2_;
            l *= -1;
            if (l < 0)l += grid_n3_;
            m = j + grid_n1_ * (k + grid_n2_ * l);
            B[m].r = p;
            B[m].i = -q;
            // debug: up to here seems ok
        }
    }
    
    int dims[3];
    
    // note inverse definition - grid_n1_ is for h-index
    dims[0] = grid_n3_;
    dims[1] = grid_n2_;
    dims[2] = grid_n1_;

    kiss_fftnd_cfg fwd_plan = kiss_fftnd_alloc(dims, 3, 0, 0, 0);
    // this overwrites B with FFT data
    kiss_fftnd(fwd_plan, B, B);
    free(fwd_plan);
    // make room for the map and difference map
    map_.resize(n5);
    if (verbosity_ > 0) {
        std::cout << "---> FFT done. Filling " << n5 << " map points into map\n";
    }

    // map min and max are only used in sxfft.f for scaling the output data when written to a textfile
    for (size_t i = 0; i < n5; i++) {
        map_[i] = B[i].r;
    }
    // ensure map statistics are being calculated
    mapstats();
    free(B);
}

/**
 * Merges sorted list of reflections on identical indices,
 * and converts I to F by weighted sqrt
 */
void sxfft::merge_data() {
    int n = -1; // count number of reflections after merging
    size_t idx = 0;
    while (idx < fcfdata_.size()) {
        int t = 0; // number of equivalent reflections
        double u = 0.; // sums up Imeas
        double v = 0.; // sums up sigmaImeas statistically correct
        double z = 0.; // sums up Fcalc
        double p = 0.; // store an arbitrary phase within each equivalence class
        int m;
        size_t k = idx;
        // data have been sorted; loop through identical indices and merge
        while ((idx < fcfdata_.size()) &&
                (fcfdata_[idx].hkl() == fcfdata_[k].hkl())) {
            ++t;
            u += fcfdata_[idx].Imeas();
            v += 1. / (fcfdata_[idx].sigImeas() * fcfdata_[idx].sigImeas());
            z += fcfdata_[idx].Fcalc();
            p = fcfdata_[idx].phicalc();
            idx++;
        }
        m = n + 1;
        // store merged reflections at next available slot;
        // convert intensities (Imeas) to amplitudes
        fcfdata_[m].Imeas(std::sqrt(fmax(0., u / t)));
        fcfdata_[m].sigImeas(sqrt(fcfdata_[m].Imeas() * fcfdata_[m].Imeas() + sqrt(1. / v)) - fcfdata_[m].Imeas());
        fcfdata_[m].Fcalc(z / t);
        fcfdata_[m].phicalc(p);
        fcfdata_[m].hkl(fcfdata_[k].hkl());
        n = m;
        // store the smallest index
    }
    n++; // this leaves room for F000
    // data have been merged - only keep the first n
    fcfdata_ = std::vector<FCFitem>(fcfdata_.begin(), fcfdata_.begin() + n);
}

/**
 * find equivalent reflections, set to standard indices, and transform 
 * phases by translational part of symop
 */
void sxfft::standard_hkl() {
    for (auto it = fcfdata_.begin(); it != fcfdata_.end(); ++it) {
        // original index: i
        double u = it->hkl().h();
        double v = it->hkl().k();
        double w = it->hkl().l();
        int mh = it->hkl().h();
        int mk = it->hkl().k();
        int ml = it->hkl().l();
        double p, q = M_PI / 180. * (it->phicalc());
        it->phicalc(std::fmod(4 * M_PI + q, 2 * M_PI));

        for (auto itsym = symops_.begin(); itsym != symops_.end(); ++itsym) {
            // original running index: k
            int nh = (int) (u * (*itsym)(0, 0) + v * (*itsym)(1, 0) + w * (*itsym)(2, 0));
            int nk = (int) (u * (*itsym)(0, 1) + v * (*itsym)(1, 1) + w * (*itsym)(2, 1));
            int nl = (int) (u * (*itsym)(0, 2) + v * (*itsym)(1, 2) + w * (*itsym)(2, 2));

            double t = 1.0;

            // ensure the smallest non-zero index i positive
            // my sorting is h > k > l - should the following be inverted?
            if ((nh < 0) || ((nh == 0)&&(nk < 0)) || ((nh == 0)&&(nk == 0)&&(nl < 0))) {
                nh *= -1; // nh == 0 or now positive
                nk *= -1; // nk == 0 or now positive
                nl *= -1; // nh 
                t = -1.0;
            }
            // sort by indices to standard order
            if ((nh < mh) || ((nh == mh)&&(nk < mk)) || ((nh == mh)&&(nk == mk)&&(nl <= ml))) continue;
            mh = nh;
            mk = nk;
            ml = nl;
            // phase transform with translational part
            p = 0.0;
            if (it->Fcalc() > -998.0) {
                p = u * (*itsym)(0) + v * (*itsym)(1) + w * (*itsym)(2);
            }
            // we are still in degree, not rad
            it->phicalc(std::fmod(719.99 + t * std::fmod(q - 360 * p, 360), 360) + 0.01);
        }
        it->hkl(HKL(mh, mk, ml));
    }

}

int sxfft::magicTop(int j) const {
    std::array<int, 143>::const_iterator it;
    for (it = magicFC_.begin(); it < magicFC_.end() && *it < j; ++it) {
    }

    if (it == magicFC_.end()) {
        std::cout << "*** Error: Resolution for grid is excessively high > "
                << magicFC_.back() << "\n";
        throw std::runtime_error("Excessively large grid resolution for n1");
    }
    return (*it);
}

/**
 * write map data to text file - for comparison with original sxfft.f
 * @param outfile
 */
void sxfft::asciimap(const std::string outfile) const {
    std::ofstream outp(outfile);
    if (!outp.is_open()) {
        std::cout << "*** -> Error: cannot open file " << outfile << "for writing"
                << std::endl;
        throw myExcepts::FileIO("Output for ASCII mao");
    }
    double scale = 9999.0 / (maxpix_ - minpix_);

    outp << std::setw(5) << (centrosymmetric_ ? 1 : 0)
            << std::setw(5) << grid_n1_
            << std::setw(5) << grid_n2_
            << std::setw(5) << grid_n3_
            << std::setw(12) << std::setprecision(8)
            << std::fixed << 1. / scale
            << std::setw(12) << std::setprecision(6) << minpix_
            << std::setw(12) << std::setprecision(8) << map_variance_
            << std::endl;
    // todo: for centrosymmetric map, only write half map to be consistent with 
    // SXFFT.f
    for (auto it = map_.begin(); it != map_.end();) {
        // write 16 integers per line
        for (int i = 0; i < 16; ++i) {
            int val = std::round(scale * (*it - minpix_));
            outp << std::setw(5) << val;
            ++it;
        }
        outp << std::endl;
    }
    outp.close();
}

/**
 * find maximum and minimum value in map array, as well as mean and variance
 */
void sxfft::mapstats() {
    maxpix_ = -std::numeric_limits<double>::infinity();
    minpix_ = std::numeric_limits<double>::infinity();
    double dm, ds;
    for (auto it = map_.begin(); it != map_.end(); ++it) {
        if (maxpix_ < *it) maxpix_ = *it;
        if (minpix_ > *it) minpix_ = *it;
        dm += *it;
        ds += (*it)*(*it);
    }
    ds = 2.0 * ds / map_.size() - (2.0 * dm / map_.size());
    map_mean_ = dm / map_.size();
    map_variance_ = std::sqrt(ds);
}
