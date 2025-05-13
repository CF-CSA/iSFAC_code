/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ResFile.h
 * Author: tg
 *
 * Created on September 3, 2023, 8:56 PM
 */

#ifndef RESFILE_H
#define RESFILE_H

#include "Vec3.h"

#include <vector>
#include <string>


/**
 * This is a simplified parser for SHELX res file. It is only
 * interested in the atom coordinates (for now), and therefore e.g. does not 
 * take proper care of continuation lines, or RESI, ...
 * ResFile produces a list of coordinates present in the res-file
 * It also reads symmetry information (CELL, LATT, and SYMM) to confirm with 
 * FCF file
 * See project cellopt for a more sophisticated parser
 */
class ResFile {
private:
    std::string filename_;
    std::vector<std::string> resfile_;
    std::vector<std::string> sfacs_;
    std::vector<Vec3> atomcoordinates_;
    std::vector<std::string> atomnames_;
    std::vector<short> atomsfacs_;
    
    float lambda_;
    float a_, b_, c_;
    // angles are kept as angles here
    float alpha_, beta_, gamma_;
    
    short lattice_;
    std::vector<std::string> symmcards_;
    bool centrosymmetric_;
    bool read_qs_;
    
    short verbosity_;

    int addsfac(const std::string& sfacvalue);
    bool is_xcmd(std::string& word) const;
    std::string upcase(const std::string& word) const;

public:
    ResFile(const std::string& filename, bool read_qs, short verbosity);
    ResFile(const ResFile& orig) = default;
    ~ResFile() = default;

    //! read in lines from res file
    int readres();
    //! extract what we are interested in
    int resinfo();
    //! extract all fractional coordinats
    int getxyzs();
    
    std::string filename() const { return filename_; }
    std::vector<Vec3> coordinates() const { return atomcoordinates_; }
    std::vector<std::string> atomnames() const { return atomnames_; }
    std::vector<short> atomsfacs() const { return atomsfacs_; }
    
    unsigned int sfac2Z(const std::string& sfac) const;
};

#endif /* RESFILE_H */

