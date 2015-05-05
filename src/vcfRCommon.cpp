#include "vcfRCommon.h"
#include <string>

void vcfRCommon::strsplit(std::string& mystring, std::vector<std::string>& vec_o_strings, char& split){

  int start = 0;
  int i=0;

  for(i = 1; i < mystring.size(); i++){
    if( mystring[i] == split){
      if( mystring[i] == split & mystring[i-1] == split){
        vec_o_strings.push_back("NA");
      } else {
        std::string temp = mystring.substr(start, i - start);
        vec_o_strings.push_back(temp);
      }
      start = i+1;
      i = i+1;
    }
  }

  // Handle the last element.
  std::string temp = mystring.substr(start, i - start);
  vec_o_strings.push_back(temp);
  
}


