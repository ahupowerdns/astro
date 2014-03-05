#pragma once
#include <string>
#include <vector>

struct Observation
{
  double t;
  double rflux;
  double flux;
  double fixedFlux;
  bool operator<(const Observation& rhs) const 
  {
    return t < rhs.t;
  }

  bool operator<(double theirT) const 
  {
    return t < theirT;
  }

};

class KeplerLightCurve
{
public:
  void addFits(const std::string& fname);
  void addTxt(const std::string& fname);
  void sort();
  void removeJumps();
  void removeDC();
  void plot(const std::string& fname);
  void dump(const std::string& fname);
  std::vector<Observation> d_obs;
};
