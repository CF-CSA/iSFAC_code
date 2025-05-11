/* 
 * File:   fcffile.cpp
 * Author: tg
 * 
 * Created on July 9, 2012, 8:05 PM
 */

#include "FCFfile.h"
#include "myExceptions.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>


FCFfile::FCFfile(const std::string& filename, const double& f000,
        int verbosity)
: filename_(filename), verbosity_(verbosity) {
    try {
        read_fcf(filename, f000);
    }
    catch (...) {
        throw;
    }
}

void FCFfile::read_fcf(const std::string& filename, const double& f000) {
    std::ifstream inp (filename);
    
    // check whether .fcf suffix is present
    if (!inp.is_open()) {
        std::string fn (filename);
        fn += ".fcf";
        inp.open(fn.c_str());
        if (!inp.is_open()) {
            std::cerr << "*** Error: unable open "
                    << filename 
                    << " or " << fn << "\n";
            throw myExcepts::FileMissing(filename + " as FCF file");
        return;
        }
        /* update filename to what could have been opened */
        this->filename(fn);
    }
    else if (verbosity_ > 0) {
        std::cout << "---> Reading FCF file " << filename << '\n';
    }

    read_fcfheader(inp, f000);
    read_fcfdata(inp);
    if (verbosity_ > 1) {
        std::cout << "---> Number of reflections from FCF file: " 
                << fcfdata_.size() << '\n';
    }
    
    inp.close();
}

void FCFfile::read_fcfheader(std::ifstream& inp, const double& f000) {
    const unsigned int linelength (81);
    char line[linelength];

    int listcode (-1);
    double Fcalcmax;
    float F000;

    /*
     *  some initialisations
     */
    float a, b, c, alpha, beta, gamma;
    float dhighres = -1.0;
    std::vector<Int3x3> symops;
    
    while (! inp.fail() ) {
        inp.getline(line, linelength);
        std::istringstream conv(line);
        std::string dummy;
        conv >> dummy;
        // last keyword: break and read data
        /* some short-cuts */
        if (dummy.empty() || dummy.at(0) == '#') continue;
        if (dummy == "_refln_phase_calc") break;
        if (dummy == "_shelx_refln_list_code") {
            conv >> listcode;
            if (listcode != 6) {
                std::ostringstream ostr;
                ostr << listcode;
                throw myExcepts::FCF::ListCode(ostr.str().c_str());
            }
            else continue;
        }
        if (dummy == "_shelx_F_calc_maximum") {
            conv >> Fcalcmax; 
            continue;
        }
        if (dummy == "_exptl_crystal_F_000") {
            if (std::isnan(f000)) {
                conv >> F000;
                if (verbosity_ > 0) {
                    std::cout << "---> F(000) set from FCF header to " << F000 << '\n';
                }
            }
            else {
                F000 = f000;
                if (verbosity_ > 0) {
                    std::cout << "---> F(000) set from command line to " << F000 << '\n';
                }
            }
            continue;
        }
        if (dummy == "_reflns_d_resolution_high") {
            conv >> dhighres; 
            continue;
        }
        if (dummy == "_cell_length_a") {
            conv >> a; 
            continue;
        }
        if (dummy == "_cell_length_b") {
            conv >> b; 
            continue;
        }
        if (dummy == "_cell_length_c") {
            conv >> c; 
            continue;
        }
        if (dummy == "_cell_angle_alpha") {
            conv >> alpha; 
            continue;
        }
        if (dummy == "_cell_angle_beta") {
            conv >> beta; 
            continue;
        }
        if (dummy == "_cell_angle_gamma") {
            conv >> gamma; 
            continue;
        }
        if (dummy == "_symmetry_equiv_pos_as_xyz" ||
                dummy == "_space_group_symop_operation_xyz") {
            try {
                readsymop(inp, symops);
            }
            catch (myExcepts::FCF::Symop& e){
                continue;
            }
            catch(...) {
                throw;
            }
            continue;
        }
    }
    if (inp.fail() || inp.eof()) {
        throw myExcepts::FileIO("Error reading FCF file header from file " + filename_ );
    }
    if (verbosity_ > 1) {
        std::cout << "---> Info from FCF header:\n"
                << "---> Unit cell parameters: " << a << ' ' << b << ' ' 
                << c << ' ' << alpha << ' ' << beta << ' ' << gamma << '\n'
                << "---> dmin = " << dhighres << " F000 = " << F000 
                << ", #symops: " << symops.size() << "\n"; 
    }
    fcfinfo_ = FCFInfo(a, b, c, alpha, beta, gamma, dhighres, F000, symops, verbosity_);
    // add 000 - experimental
    FCFitem datum = { HKL(0,0,0), F000*F000, std::abs(F000), F000, 0.0  };
    fcfdata_.push_back(datum);
}

void FCFfile::read_fcfdata(std::ifstream& inp) {
    /* data read from LIST 6 - fcf file */
    int h,k,l;
    double Imeas; //! actually F^2
    double sigImeas;
    double Fcalc;
    float  phicalc;

    // finally read the data
    if (verbosity_ > 1) {
        std::cout << "---> Reading FCF data\n";
    }
    while (inp.good()) {
        inp >> h >> k >> l
                >> Imeas
                >> sigImeas
                >> Fcalc
                >> phicalc;
        if (inp.eof()) {
            break;
        }
        if (inp.fail()) {
            throw myExcepts::FileIO("Error reading FCF file header from file " + filename_ );
        };
        HKL hkl (h,k,l);
        FCFitem datum = {hkl, Imeas, sigImeas, Fcalc, phicalc};
        fcfdata_.push_back(datum);
    }
    
}

/**
 * This expects inp to be placed just after _space_group_symop_operation_xyz
 * and the next lines are symmetry operators
 * @param inp
 */
void FCFfile::readsymop(std::ifstream& inp, std::vector<Int3x3>& symops) {
    Int3x3 R;
    const unsigned int linewidth (81);
    char line[linewidth];
    unsigned int firsttick, lasttick;
    unsigned int firstcomma, lastcomma;
    
    inp.getline(line, linewidth);
    std::string symop (line);
    std::string entry;
    
    /* sanity checks */
    if (std::count(symop.begin(), symop.end(), '\'') != 2 || 
            std::count(symop.begin(), symop.end(), ',') != 2) {
        throw myExcepts::FCF::Symop(symop.c_str());
    }
    
    // full expression between first and last tick
    firsttick = symop.find('\'');
    lasttick = symop.rfind('\'');
    // firsttick - firstcomma: x', firstcomma - lastcomma: y', lastcomma - lasttick: z'
    firstcomma = symop.find(',');
    lastcomma = symop.rfind(',');

    std::string symopString = symop.substr(firsttick + 1, lasttick - firsttick - 1);
    symopsStrings_.push_back(symopString);

    /* extract first row data */
    entry = symop.substr(firsttick + 1, firstcomma - firsttick - 1);
    tosymop(entry, R, 0);
    /* extract second row data */
    entry = symop.substr(firstcomma + 1, lastcomma - firstcomma - 1);
    tosymop(entry, R, 1);
    /* extract third row data */
    entry = symop.substr(lastcomma + 1, lasttick - 1);
    tosymop(entry, R, 2);

        symops.push_back(R);
    try {
        readsymop(inp, symops);
    }
    catch (myExcepts::FCF::Symop& s) {
        throw;
    }
}

/**
 * converts a string to symop including translational part
 * @param symstring
 * @param symop
 * @param row
 */
void FCFfile::tosymop (const std::string& symstring, Int3x3& symop, int row) const{
    if (symstring.find("-x") != std::string::npos) {
        symop(row, 0) = -1;
    } else if (symstring.find("x") != std::string::npos) {
        symop(row, 0) = 1;
    }
    if (symstring.find("-y") != std::string::npos) {
        symop(row, 1) = -1;
    } else if (symstring.find("y") != std::string::npos) {
        symop(row, 1) = 1;
    }
    if (symstring.find("-z") != std::string::npos) {
        symop(row, 2) = -1;
    } else if (symstring.find("z") != std::string::npos) {
        symop(row, 2) = 1;
    }
    size_t dividx = symstring.find('/');
    
    if ( dividx != std::string::npos) {
        std::istringstream inp (symstring.substr(dividx-2,4));
        int nominator, denominator;
        char slash;
        inp >> nominator >> slash >> denominator;
        symop(row) = 1.0*nominator/denominator;
    }
}

