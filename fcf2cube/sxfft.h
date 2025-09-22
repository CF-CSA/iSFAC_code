/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   sxfft.h
 * Author: tg
 *
 * Created on July 29, 2023, 9:40 PM
 */

#ifndef SXFFT_H
#define SXFFT_H

#include "FCFitem.h"
#include "FCFInfo.h"

#include <array>

/**
 * Fast Fourier Transform for FCF data; derived from original sxfft.h as well as  ShelXle 
 * code (which is most likely based on sxfft
 */
class sxfft {
    FCFInfo fcfinfo_;
    std::vector<FCFitem> fcfdata_;
    // 143 fixed entries for fast FFT grids
    const std::array<int, 143> magicFC_; 

    // storage for map
    std::vector<double> map_;
    double maxpix_, minpix_;
    double map_mean_;
    double map_variance_;
    int maptype_;
    int grid_n1_, grid_n2_, grid_n3_;
    
    // reciprocal space symmetry operators
    std::vector<Int3x3> symops_;
    
    int verbosity_;
    
    //! convert to some standard hkls
    void standard_hkl();
    
    //! merge sorted fcfdata_ and convert Imeas to Fmeas;
    void merge_data();

    //! find smallest entry of magic coefficients greater than j
    int magicTop(int j) const;
    
    
    void mapstats();
    
public:
    sxfft(const std::vector<FCFitem>& fcfdata, const FCFInfo& fcfinfo, int maptype, int verbosity);
    sxfft() = delete;
    ~sxfft() = default;
    // carry out fft - map data are stored in map_Fo and map_Delta_Fo_Fc
    void fft(const double& weakWeight, const double& gridresol);
    void asciimap(const std::string outfile) const;
    void datamap (const std::string outfile, const std::vector<double>& Br) const;
    
    std::vector<double> map() const { return map_; }
    int gridn1() const { return grid_n1_; }
    int gridn2() const { return grid_n2_; }
    int gridn3() const { return grid_n3_; }

};

#endif /* SXFFT_H */

