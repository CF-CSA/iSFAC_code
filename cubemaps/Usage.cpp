#include "Usage.h"

#include <iostream>
#include <cstdio>

/**
 *print a short welcome notice
 */
std::ostream& hello(std::ostream& outp) {
    outp << "!---------------------------------------------------------------!\n"
            << "! cubemaps: Compute Pearson CC between two cube files        !\n"
            << "! Copyright:    Tim Gruene, 2024                             !\n"
            << "!------------------------------------------------------------!\n"
            ;
    return outp;
}

void usage() {
    std::cout << "!---------------------------------------------------------------!\n"
              << "! Usage: cubemaps [Options]                                     !\n"
              << "!        maps one cube file onto another and computes Pearson CC!\n"
	      << "!        between those two                                      !\n"
              << "!         OPTIONS:                                              !\n"
	      << "!        -R     :  reference (non-moving) map                   !\n"
	      << "!        -M     :  moving map                                   !\n"
              << "!        -g d   : set grid spacing for VdW surface [0.4A]       !\n"
              << "!        -v num : set verbosity to num                          !\n"
              << "!        -h / -?: print this help message and exit              !\n"
              << "!---------------------------------------------------------------!\n"
	      << "\n";
   
}
