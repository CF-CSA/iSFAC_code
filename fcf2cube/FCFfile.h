/* 
 * File:   fcffile.h
 * Author: tg
 *
 * Created on July 9, 2012, 8:05 PM
 */

#ifndef FCFFILE_H
#define	FCFFILE_H

#include "Reflex.h"
#include "Int3x3.h"
#include "FCFitem.h"
#include "FCFInfo.h"

#include <fstream>
#include <list>
#include <string>

/**
 * Possibilities for fcf-files:
 * - 
 */
class FCFfile {
private:
    
    std::string filename_;

    FCFInfo fcfinfo_;
    std::vector<FCFitem> fcfdata_;
    
    int verbosity_;
    
    std::vector<std::string> symopsStrings_;
    
    void readsymop (std::ifstream& inp, std::vector<Int3x3>& symops);
    void tosymop (const std::string& symstring, Int3x3& symop, int row) const;

    void read_fcf (const std::string& filename, const double& f000);
    void read_fcfheader (std::ifstream& inp, const double& f000);
    void read_fcfdata (std::ifstream& inp);

    void filename (const std::string& fn) { filename_ = fn; }
public:
    FCFfile(const std::string& filename, const double& f000, int verbosity);

    ~FCFfile(){}
    
    std::vector<FCFitem> fcfdata () const { return fcfdata_; }
    FCFInfo fcfinfo() const { return fcfinfo_; }
    
    
};

#endif	/* FCFFILE_H */

