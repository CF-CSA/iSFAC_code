/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Parser.h
 * Author: tg
 *
 * Created on May 9, 2025, 9:19 PM
 */

#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <vector>

class Parser {
private:
    std::string fcffile_;
    std::string resfile_;
    std::string cubefile_;
    std::string name_;
    std::vector<std::string> non_options_;
    double f000_;
    bool read_qs_;
    unsigned short verbosity_;
  
    template <typename T> bool getoption(const std::string& option,
            const std::string& opt, T& optval, int& idx, int argc,
            const char* const argv[]);
    

public:
    Parser(int argc, const char* const argv[]);
    ~Parser(){};
    
    //! return name of fcffile
    std::string fcffile() const { return fcffile_; }
    std::string resfile() const { return resfile_; }
    //! return name of cubefule
    std::string cubefile() const { return cubefile_; };
    
    double f000() const { return f000_; }
    bool   read_qs() const { return read_qs_; }
    
    unsigned char verbosity() const { return verbosity_; }
};

#endif /* PARSER_H */

