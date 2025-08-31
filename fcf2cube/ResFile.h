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
    struct XAtom {
        std::string name_ = "";
        unsigned short sfac_idx_ = 0;
        double x_, y_, z_;
        double occupancy_;
        bool skip_ = false;
    } xatom;
    struct SfacZ {
        std::string element_;
        unsigned short Z_;
    } sfacz;
    
    //! original input file name
    std::string filename_;
    
    //! all lines of the resfile 
    std::vector<std::string> resfile_;
    
    //! atomic numbers of sfacs_
    std::vector<SfacZ> sfacsZ_;
    std::vector<XAtom> xatoms_;
    
    double lambda_;
    double a_, b_, c_;
    // angles are kept as angles (degree) here
    double alpha_, beta_, gamma_;
    
    //! unit cell vectors
    Vec3 A_, B_, C_;
    
    short lattice_;
    std::vector<std::string> symmcards_;
    bool centrosymmetric_;
    bool read_qs_;
    
    short verbosity_;

    //! read in lines from res file
    int readres();
    //! process header / commands from res file
    int resheader();
    //! extract all fractional coordinats
    int getxyzs();
    
    int addsfac(const std::string& sfacvalue);
    bool is_xcmd(std::string& word) const;
    std::string upcase(const std::string& word) const;
    // true if SFAC for this XAtom is a proper PSE element name
    bool proper_PSE_element(const XAtom& xatom) const;

public:
    ResFile(const std::string& filename, bool read_qs, short verbosity);
    ResFile(const ResFile& orig) = default;
    ~ResFile() = default;

    //! get input filename
    std::string filename() const { return filename_; }
    
    //! convert the element part of SFAC to PSE Z number
    unsigned int sfac2Z(const std::string& sfac) const;

    //! getters
    Vec3 A() const { return A_; }
    Vec3 B() const { return B_; }
    Vec3 C() const { return C_; }
    
    //! compute a tight bounding box, converted to Angstroem, no margin added
    std::pair<Vec3, Vec3> bbox_frac() const;
    
    //! mark atom pairs too close to one another to skip
    void skip_close_atoms();
    
    //! cube formatted string of atom coordinates
    std::string atom_list_for_cube(int& num_atoms) const;

};

#endif /* RESFILE_H */

