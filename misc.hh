#pragma once
#include <stdio.h>
#include <string>
#include <stdint.h>
#include <vector>

/** Rapid estimator of variance and mean of a series of doubles. 
    API compatible with a, sadly, far slower boost::accumulator_set
    doing the same thing*/
class VarMeanEstimator
{
public:
  VarMeanEstimator() : N(0), xTot(0), x2Tot(0) {}
  void operator()(double val) 
  {
    ++N;
    xTot += val;
    x2Tot += val*val;
  }
  friend double mean(const VarMeanEstimator& vme);
  friend double variance(const VarMeanEstimator& vme);
private:
  uint64_t N;
  double xTot;
  double x2Tot;
};

//! extract 'mean' from a VarMeanEstimator
inline double mean(const VarMeanEstimator& vme)
{
  return vme.xTot/vme.N;
}

//! extract 'variance' from a VarMeanEstimator
inline double variance(const VarMeanEstimator& vme) 
{
  return (vme.x2Tot - vme.xTot*vme.xTot/vme.N)/vme.N;
}
bool stringfgets(FILE* fp, std::string* line);

template<typename T>
std::vector<std::pair<typename T::iterator, typename T::iterator> >
splitRange(T& container, int parts)
{
  typedef std::pair<typename T::iterator, typename T::iterator> range_t;
  std::vector<range_t > ret;
  range_t range;
  typename T::size_type stride = container.size()/parts;

  for(int n=0; n < parts; ++n) {
    if(!n)
      range.first = container.begin();
    else
      range.first = range.second;

    if(n+1 == parts)
      range.second = container.end();
    else
      range.second = range.first + stride;
      
    ret.push_back(range);
  }
  
  return ret;
}

void normalize(std::vector<double>& values);
uint64_t filesize(const char* name);
