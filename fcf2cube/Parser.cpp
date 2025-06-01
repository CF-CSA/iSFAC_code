/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Parser.cpp
 * Author: tg
 * 
 * Created on May 9, 2025, 9:19 PM
 */

#include "Parser.h"
#include "myExceptions.h"
#include "Utils.h"

#include <iostream>
#include <cmath>
#include <vector>

Parser::Parser(int argc, const char* const argv[]) :
name_(""),
fcffile_(""),
resfile_(""),
cubefile_(""),
f000_(std::nan("")),
maptype_(1),
read_qs_(false),
margin_(0.25),
weight_weaks_(1.0),
hklgridres_(5.0),
verbosity_(1) {

    std::string outfile;
    for (int i = 1; i < argc; ++i) {
        std::string option = argv[i];
        if (option.substr(0, 2) == "-h" || option.substr(0, 2) == "-?") {
            throw myExcepts::Usage("Help message");
            return;
        }
        if (option.substr(0, 2) == "-q") {
            read_qs_ = true;
            continue;
        }

        if (getoption(option, "-f", fcffile_, i, argc, argv)) continue;
        if (getoption(option, "-r", resfile_, i, argc, argv)) continue;
        if (getoption(option, "-o", cubefile_, i, argc, argv)) continue;
        if (getoption(option, "-0", f000_, i, argc, argv)) continue;
        if (getoption(option, "-m", maptype_, i, argc, argv)) continue;
        if (getoption(option, "-b", margin_, i, argc, argv)) continue;
        if (getoption(option, "-v", verbosity_, i, argc, argv)) continue;

        // when reaching this points, unknown option character
        if (option.at(0) == '-') {
            std::cout << "Unknown option string " << option << std::endl;
            throw myExcepts::Usage(option.c_str());
        }
        // non-option term, e.g. project root
        non_options_.push_back(argv[i]);
    }
    if (non_options_.size() == 1) {
        // name of project, ignore fcf-file root
        name_ = non_options_.front();
        if (verbosity_ > 0) {
            std::cout << Utils::prompt(1) << "Project name " << name_
                    << ", will read " << name_ + ".fcf"
                    << " and " << name_ + ".hkl" << '\n'
                    << "     output will be written to " << name_ + ".cub\n";
        }
        if (fcffile_.empty()) {
            fcffile_ = name_ + ".fcf";
        }
        if (resfile_.empty()) {
            resfile_ = name_ + ".res";
        }
        if (cubefile_.empty()) {
            cubefile_ = name_ + ".cub";
        }

    }
    if (fcffile_.empty() || resfile_.empty() || cubefile_.empty()) {
        if (verbosity_ > 0) {
            std::cout << Utils::prompt(1) << "Missing option for fcf-, res-, or cubefile\n"
                    << "      Either specifiy project root, or explicitly all three files\n";
        }
        // error
        throw myExcepts::Usage("not all filenames set");
    }
}

/**
 * Interpret option string
 * @param option
 * @param opt
 * @param optval
 * @param idx
 * @param argc
 * @param argv
 * @return 
 */
template <typename T> bool Parser::getoption(const std::string& option,
        const std::string& opt, T& optval, int& idx, int argc,
        const char* const argv[]) {

    if (option.substr(0, 2) != opt) return false;

    std::string convertee;
    if (option.length() > 2) {
        convertee = option.substr(2, option.length());
    } else {
        if (idx + 1 >= argc) {
            std::cout << "*** Error: option " << argv[idx] << " requires "
                    << "an argument.\n";
            throw myExcepts::Usage("Option requires argument");
        }
        convertee = argv[idx + 1];
        ++idx;
    }
    std::istringstream inp(convertee);
    inp >> optval;
    if (inp.fail()) {
        std::cout << "*** Error: Cannot convert " << convertee
                << std::endl;
        throw myExcepts::Usage("Conversion of command line parameter");
    }
    return true;
}

