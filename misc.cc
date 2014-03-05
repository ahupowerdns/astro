#include "misc.hh"
#include <vector>
#include <string.h>

//! read a line of text from a FILE* to a std::string, returns false on 'no data'
bool stringfgets(FILE* fp, std::string* line)
{
  char buffer[1024];   
  line->clear();
  
  do {
    if(!fgets(buffer, sizeof(buffer), fp))
      return !line->empty();
    
    line->append(buffer);
  } while(!strchr(buffer, '\n'));
  return true;
}

//! make sure the average of a vector is 1
void normalize(std::vector<double>& values)
{
  double tot = 0;
  for(const auto& val : values)
    tot+=val;
  const double average = tot/values.size();

  for(auto& val : values)
    val /= average;  
}

