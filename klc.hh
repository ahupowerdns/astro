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


struct CSplineSignalInterpolator
{
  explicit CSplineSignalInterpolator(const std::vector<Observation>& obs)
    : d_obs(&obs)
  {
    double p, qn, sig, un;

    int n=obs.size();
    std::vector<double> u(n-1);
    d_y2.resize(n);

    d_y2[0]=u[0]=0.0;
    for (int i=1;i<n-1;i++) {
      sig=(obs[i].t - obs[i-1].t)/(obs[i+1].t - obs[i-1].t);
      p=sig*d_y2[i-1]+2.0;
      d_y2[i]=(sig-1.0)/p;
      u[i]=(obs[i+1].fixedFlux - obs[i].fixedFlux)/(obs[i+1].t-obs[i].t) - 
	(obs[i].fixedFlux - obs[i-1].fixedFlux)/(obs[i].t - obs[i-1].t);
      u[i]=(6.0*u[i]/(obs[i+1].t- obs[i-1].t)-sig*u[i-1])/p;
    }
    
    qn=un=0.0;
    d_y2[n-1]=(un-qn*u[n-2])/(qn*d_y2[n-2]+1.0);
    for (int k=n-2;k>=0;k--)
      d_y2[k]=d_y2[k] * d_y2[k+1]+u[k];
  }

  double operator()(double t) const {

    auto hi=lower_bound(d_obs->begin(), d_obs->end(), t);
    if(hi == d_obs->begin() || hi == d_obs->end())
      return 0;
    
    auto lo = hi - 1;
    if(lo == d_obs->begin()) {
      //      cout << t << '\t' << 0 << '\n';
      return 0;
    }

    double h=hi->t - lo->t;
    if(h > 10) {// makes little sense to interpolate over days
      //      cout << t << '\t' << 0 << '\n';
      return 0;
    }

    //    if (h == 0.0) nrerror("Bad xa input to routine splint");
    double a=(hi->t - t)/h;
    double b=(t - lo->t)/h;
    double ret= a*lo->fixedFlux + b*hi->fixedFlux +((a*a*a-a)* d_y2[lo - d_obs->begin()]
						    +(b*b*b-b)*d_y2[hi - d_obs->begin()])*(h*h)/6.0;
    // cout << t << '\t\ << ret << '\n';
    return ret;
  }

  const std::vector<Observation>* d_obs;
  std::vector<double> d_y2;
};
