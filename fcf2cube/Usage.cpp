#include "Usage.h"

#include <iostream>
#include <cstdio>

/**
 *print a short welcome notice
 */
std::ostream& hello(std::ostream& outp) {
    outp 
            << "!---------------------------------------------------------------!\n"
            << "! fcf2cube: Convert SHELXL FCF file to map in cube format       !\n"
            << "!           Covering all atoms in RES file                      !\n"
            << "! Copyright:    Tim Gruene, 2025                                !\n"
            << "!---------------------------------------------------------------!\n"
            ;
    return outp;
}

void usage() {
    std::cout << "!---------------------------------------------------------------!\n"
              << "! Usage: fcf2cube [Options] [name]                              !\n"
              << "!         OPTIONS:                                              !\n"
              << "!        name      : read name.fcf, name.res and write name.cub !\n"
              << "!                    overwritten by -f | -r | -o                !\n"
	      << "!        -f file   : read FCF file                              !\n"
	      << "!        -r file   : read RES file                              !\n"
              << "!        -o file   : write Cube map to file                     !\n"
              << "!        -0 num    : set F(000) to num                          !\n" 
              << "!        -g grid   : adjust grid sampling (higher for finer) [5.0]!\n"
              << "!        -b num    : expand margin around molecule [0.35]       !\n"
	      << "!                    (add fraction to space diagonal)           !\n"
              << "!        -v num : set verbosity to num                          !\n"
              << "!        -h / -?: print this help message and exit              !\n"
              << "!---------------------------------------------------------------!\n"
	      << "\n";
   
}
