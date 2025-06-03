/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ResFile.cpp
 * Author: tg
 * 
 * Created on September 3, 2023, 8:56 PM
 */

#include <fstream>
#include <iomanip>
#include <cstring>

#include "ResFile.h"
#include "myExceptions.h"
#include "FCFInfo.h"
#include "Utils.h"
#include "defines.h"

ResFile::ResFile(const std::string& filename, bool read_qs = true, short verbosity = 1) :
filename_(filename),
centrosymmetric_(false),
read_qs_(read_qs),
verbosity_(verbosity) {
    // add empty sfac so that shelxl index numbers match
    SfacZ s;
    s.Z_ = 0;
    s.element_ = "";
    sfacsZ_.push_back(s);
    // read res file
    readres();
    // process commands from header
    resheader();
    // extract atoms with sfac and coordinates
    getxyzs();
    // compute unit cell vectors
    std::tie (A_, B_, C_) = Utils::unit_cell_vectors(a_, b_, c_, alpha_, beta_, gamma_);
}

/**
 * read in res file line by line, storing content in vector
 * resfile_
 * @return number of lines read
 */
int ResFile::readres() {
    std::ifstream inp(filename_);
    if (!inp.is_open()) {
        if (verbosity_ > 0) {
            std::cout << "---> Error opening file " << filename_ << " for reading.\n"
                    << "    check file permisions! Does file exist?\n";
        }
        throw myExcepts::FileIO("Cannot open " + filename_);
    }
    while (!inp.eof()) {
        std::string myline;
        std::getline(inp, myline);
        if (inp.fail() || inp.eof()) break;
        else {
            resfile_.push_back(myline);
        }
    }
    inp.close();
    return resfile_.size();
}

/**
 * extract info for
 * cell
 * lattice
 * symmetry cards
 * @return 
 */
int ResFile::resheader() {
    for (auto it = resfile_.begin(); it != resfile_.end(); ++it) {
        if (it->empty()) continue;
        std::string key = it->substr(0, 4);
        if (key == "CELL") {
            std::istringstream conv(it->substr(4));
            conv >> lambda_ >> a_ >> b_ >> c_ >> alpha_ >> beta_ >> gamma_;
            if (conv.fail()) {
                if (verbosity_ > 0) {
                    std::cout << "---> Error: cannot extract cell from "
                            << conv.str() << '\n';
                }
                throw myExcepts::Format("resfile CELL");
            }
            continue;
        }
        if (key == "LATT") {
            std::istringstream conv(it->substr(4));
            conv >> lattice_;
            if (lattice_ > 0) {
                centrosymmetric_ = false;
            } else centrosymmetric_ = true;
            continue;
        }
        if (key == "SYMM") {
            const std::string mysymm = it->substr(4);
            symmcards_.push_back(mysymm);
            continue;
        }
        if (key == "SFAC") {
            addsfac(it->substr(4));
            continue;
        }
    }

    if (verbosity_ > 0) {
        std::cout << Utils::prompt(1) << "Scattering factors:";
        size_t idx(0);
        for (auto s: sfacsZ_) {
            std::cout << std::setw(4) << s.element_ << " (Z= " << s.Z_ << ")";
            ++idx;
        }
        std::cout << '\n';
    }

    return lattice_;
}

/**
 * checks whether word is n command list
 * @param word
 * @return 
 */
bool ResFile::is_xcmd(std::string& word) const {
    const std::vector<std::string> cmdlist = {"ABIN", "ACTA", "AFIX", "ANIS", "ANSC", "ANSR",
        "BASF", "BIND", "BLOC", "BOND", "BUMP", "CELL",
        "CGLS", "CHIV", "CONF", "CONN", "DAMP", "DANG",
        "DEFS", "DELU", "DFIX", "DISP", "EADP", "END ",
        "EQIV", "EXTI", "EXYZ", "FEND", "FLAT", "FMAP",
        "FRAG", "FREE", "FVAR", "GRID", "HFIX", "HKLF",
        "HTAB", "ISOR", "LATT", "LAUE", "LIST", "L.S.",
        "MERG", "MORE", "MOVE", "MPLA", "NCSY", "NEUT",
        "OMIT", "PART", "PLAN", "PRIG", "REM ", "RESI",
        "RIGU", "RTAB", "SADI", "SAME", "SFAC", "SHEL",
        "SIMU", "SIZE", "SPEC", "STIR", "SUMP", "SWAT",
        "SYMM", "TEMP", "TITL", "TWIN", "TWST", "UNIT",
        "WGHT", "WIGL", "WPDB", "XNPD", "ZERR"};

    //! checking END and REM separately, just in case
    if (upcase(word) == "END" || upcase(word) == "REM") {
        return true;
    }
    for (auto it = cmdlist.begin(); it != cmdlist.end(); ++it) {
        if (verbosity_ > 4) {
            std::cout << "---> Checking " << word << " against list of SHELX commands\n";
        }
        if (*it == upcase(word)) {
            return true;
        } else continue;
    }
    return false;
}

/**
 * extract sfac values from string sfacvalue.
 * First one is a string. If 2nd one is a number, rest is ignores, otherwise 
 * remaining strings are extracted and added
 * @param sfacvalue
 * @return 
 */
int ResFile::addsfac(const std::string& sfacvalue) {
    std::istringstream conv(sfacvalue);
    int newsfacs(0);

    std::string sfac;
    conv >> sfac;
    if (conv.fail()) {
        if (verbosity_ > 0) {
            std::cout << "---> Error: Malformatted SFAC command " << sfacvalue
                    << '\n';
        }
        throw myExcepts::Format("SFAC card " + sfacvalue);
    }
    if (verbosity_ > 2) {
        std::cout << "---> Extracted SFAC " << sfac << " from " << sfacvalue
                << '\n';
    }

    sfacz.element_ = sfac;
    sfacz.Z_ = sfac2Z(sfac);
    sfacsZ_.push_back(sfacz);
    ++newsfacs;

    float dummy;
    conv >> dummy;
    if (conv.good()) return newsfacs;
    conv.clear();
    while (conv >> sfac) {
        if (verbosity_ > 2) {
            std::cout << "---> Extracted SFAC " << sfac << " from " << sfacvalue
                    << '\n';
        }

    sfacz.element_ = sfac;
    sfacz.Z_ = sfac2Z(sfac);
    sfacsZ_.push_back(sfacz);
    }
    return newsfacs;
}

/**
 * convert SFAC to Z; assumes that E starts with a sensible chemical name, 
 * like Ca+ or O2-, followed by a non-letter
 * @param sfac
 * @return Z for this element
 */
unsigned int ResFile::sfac2Z(const std::string& sfac) const {
    char E[3] = "0\0";
    if (!isupper(sfac.at(0))) {
        if (verbosity_ > 0) {
            std::cout << Utils::error(1) 
                    << "Illegal Sfac does not start with upper case character: " << sfac << "\n";
            throw myExcepts::Format("SFAC card not upper case\n");
        }
    }
    E[0] = sfac.at(0);
    if (sfac.length() > 1 && islower(sfac.at(1))) {
        E[1] = sfac.at(1);
    }
    // search for element in string
    unsigned int Z = 0;
    for (int z = 0; z < PSE::NUM_ELEMENTS; ++z) {
        if (strcmp(E, PSE::Elements[z]) == 0) {
            Z = z;
            break;
        }
    }
    if (Z == 0) {
        if (verbosity_ > 0) {
            std::cout << Utils::error(1) 
                    << "Unknown element " << E << " in SFAC " << sfac << std::endl;
            throw myExcepts::Format("SFAC card"+sfac);
        }
    }
    if (verbosity_ > 1) {
        std::cout << Utils::prompt(1) 
                << "Z = " << Z << " extracted for sfac " << sfac << '\n';
    }
    return Z;
}

/**
 * extract all coordinates from RES file, ignores continuation lines
 * @return 
 */
int ResFile::getxyzs() {
    for (auto it = resfile_.begin(); it != resfile_.end(); ++it) {
        if (it->empty()) continue;
        std::string key = it->substr(0, 4);
        if ( (key == "END" || key == "END ") && !read_qs_) break; 
        if (is_xcmd(key)) continue;
        if (key.at(0) == ' ') continue;

        // these should be all options - this should be a atom line
        unsigned short sfac;
        float x, y, z;
        float occupancy;
        std::istringstream conv(it->substr(4));
        conv >> sfac >> x >> y >> z >> occupancy;
        if (conv.fail()) {
            std::cout << "---> Error: expected atom line but cannot extract coordinates from \n"
                    << "    \'" << conv.str() << "\'\n";
            throw myExcepts::Format("No SHELX atom line: " + conv.str());
        }
        
        xatom.name_ = key;
        xatom.sfac_idx_ = sfac;
        xatom.x_ = x;
        xatom.y_ = y;
        xatom.z_ = z;
        xatom.occupancy_ = occupancy;
        xatoms_.push_back(xatom);
    }
    return xatoms_.size();
}

std::string ResFile::upcase(const std::string& word) const {
    std::string up(word);
    for (std::string::iterator it = up.begin(); it != up.end(); ++it) {
        *it = std::toupper(*it);
    }
    return up;
}

/**
 * Compute the bounding box of the coordinates (opposite corners)
 * converted to Angstroem
 * @return 
 */
std::pair<Vec3, Vec3> ResFile::bbox_frac() const {
    std::vector<Vec3> coords;
    for (auto a: xatoms_) {
        coords.push_back(Vec3(a.x_, a.y_, a.z_));
    }
    std::pair<Vec3, Vec3> bb = Utils::bbox3D(coords, verbosity_);
    
    if (verbosity_ > 2) {
        std::cout << Utils::prompt(2) << "Bounding box determined in A:\n"
                << bb.first.x()*A_ + bb.first.y()*B_ + bb.first.z()*C_
                << " " 
                << bb.second.x()*A_ + bb.second.y()*B_ + bb.second.z()*C_
                << std::endl;
        const double s = 1./Physics::a0;
        std::cout << Utils::prompt(2) << "Bounding box determined in Bohr:\n"
                << s*bb.first.x()*A_ + s*bb.first.y()*B_ + s*bb.first.z()*C_
                << " " 
                << s*bb.second.x()*A_ + s*bb.second.y()*B_ + s*bb.second.z()*C_
                << std::endl;
    }
    return bb;
}

bool ResFile::proper_PSE_element(const XAtom& xatom) const {
   unsigned short Z = sfacsZ_[xatom.sfac_idx_].Z_;
    std::string sfacname = sfacsZ_[xatom.sfac_idx_].element_;

    if (PSE::Elements[Z] == sfacname) {
        return true;
    }
    else {
        return false;
    }
}

/**
 * make atoms closer than 0.5A to a previous one as skipped atom
 * @return 
 */
void ResFile::skip_close_atoms() {
    // mark all atoms closer than 0.5A is too close and do not include in list
    for (std::vector<XAtom>::iterator it = xatoms_.begin();
            it != xatoms_.end()-1; ++it) {
        if (it->skip_) {
            continue;
        }
        for (std::vector<XAtom>::iterator it2 = it+1; it2 != xatoms_.end(); ++it2) {
            const Vec3 delta = (it2->x_ - it->x_)*A_ 
                                + (it2->y_ - it->y_)*B_ 
                                + (it2->z_ - it->z_)*C_;
            const double d2 = delta*delta;
            if (d2 < 0.5*0.5) {
                it->skip_ = true;
            }
        }
    }
    
}
/**
 * create a string of the proper elements in atom list
 * in format suitable for CUBE file, with coordinates converted
 * to Bohr from Angstroem
 * @return 
 */
std::string ResFile::atom_list_for_cube(int& num_atoms) const{
    std::ostringstream outp;
    const double s = 1./Physics::a0;
    num_atoms = 0;
    for (auto a: xatoms_) {
        if (a.skip_ == true) continue;
        if (proper_PSE_element(a)) {
            unsigned short Z = sfacsZ_[a.sfac_idx_].Z_;
            Vec3 x = a.x_*A_ + a.y_*B_ + a.z_*C_;
            outp << std::setw(5) << Z
                    << std::setw(12) << std::setprecision(6) << std::fixed << double(Z)
                    << std::setw(12) << std::setprecision(6) << s*x.x()
                    << std::setw(12) << std::setprecision(6) << s*x.y()
                    << std::setw(12) << std::setprecision(6) << s*x.z()
                    << '\n';
		++num_atoms;
        }
    }
    return outp.str();
}
