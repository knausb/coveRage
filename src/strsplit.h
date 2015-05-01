#ifndef __STRSPLIT_H_INCLUDED__   // if x.h hasn't been included yet...
#define __STRSPLIT_H_INCLUDED__   //   #define this so the compiler knows it has been included

// Forward declared dependencies
//class mystring;
//class vec_o_strings;
//class split;

#include <string>
#include <vector>

class strsplit{
  public:
    //static 
//    strsplit();
//    void 
    strsplit(std::string&, std::vector<std::string>&, char&);
//    static void strsplit(std::string&, std::vector<std::string>&, char&);
//  static void str_split(std::string&, std::vector<std::string>&, char&);
//  static void strsplit(unsigned char[]&, std::vector<std::string>&, char&);
};

#endif 