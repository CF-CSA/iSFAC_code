/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: tg
 *
 * Created on May 9, 2025, 9:11 PM
 */

#include <iostream>
#include "Parser.h"
#include "FCFfile.h"
#include "Usage.h"
#include "ResFile.h"
#include "sxfft.h"
#include "myExceptions.h"
#include "MapValues.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    try {
        Parser parser(argc, argv);
        if (parser.verbosity() > 0) {
            hello(std::cout);
        }

        // resfile defines the grid points for the map surrouding molecule
        ResFile resfile(parser.resfile(), parser.read_qs(), parser.verbosity());
        if (parser.verbosity() > 1) {
            std::cout << resfile.atom_list_for_cube();
        }

        // read fcf file
        FCFfile fcffile(parser.fcffile(), parser.f000(), parser.verbosity());
        FCFInfo fcfinfo = fcffile.fcfinfo();
        std::vector<FCFitem> fcfdata = fcffile.fcfdata();

        if (parser.verbosity() > 0) {
            std::cout << "---> Converting FCF into Map by FFT\n";
        }
        // compute fft
        sxfft myfft(fcfdata, fcfinfo, parser.maptype(), parser.verbosity());
        myfft.fft(parser.weight(), parser.hklgridres());
        
        // create the map for the Cubefile surrounding molecule
        MapValues mapvals(myfft.map(), myfft.gridn1(), myfft.gridn2(), myfft.gridn3(), parser.verbosity());

    }    catch (myExcepts::Usage e) {
        usage();
        return -1;
    }

    return 0;
}

