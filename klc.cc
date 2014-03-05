#include "klc.hh"
#include <iostream>
#include <fstream>
#include <CCfits/CCfits>
using namespace CCfits;


using namespace std;
void KeplerLightCurve::addFits(const std::string& fname)
{
  cerr<<fname<<endl;
  CCfits::FITS pInfile(fname, Read); 

  ExtHDU& table = pInfile.extension(1);
    
  // read all the keywords, excluding those associated with columns.
    
  table.readAllKeys();
  
  std::vector <double> flux;
  table.column("PDCSAP_FLUX").read(flux, 0, table.rows());
  
  std::vector <double> rflux;
  table.column("SAP_FLUX").read(rflux, 0, table.rows());
    
  std::vector <double> jd;
  table.column("TIME").read(jd, 0, table.rows());
    
  for(unsigned int i=0; i < jd.size(); ++i) {
    if(std::isnan(jd[i]) || std::isnan(rflux[i]) || std::isnan(flux[i])) 
      continue;
    d_obs.push_back({86.0*jd[i], rflux[i], flux[i]});
  }
}

void KeplerLightCurve::addTxt(const std::string& fname)
{
  cerr<<fname<<endl;
  ifstream ifs(fname);
  if(!ifs)
    throw runtime_error("Error reading file '"+fname+"'");
  
  double t,ppm;
  string comment;
  while(!ifs.eof()) {
    if(!(ifs >> t)) {
      ifs.clear();
      getline(ifs, comment);
      cerr<<"Read: "<<comment<<endl;
      continue;
    }
    ifs >> ppm;
    cout << 1000*t << " -> "<< ppm <<endl;
    d_obs.push_back({1000*t, ppm, ppm});
  }
}

void KeplerLightCurve::sort()
{
  ::sort(d_obs.begin(), d_obs.end());
  cerr<< "Kilosecond timespan: "<<d_obs.begin()->t << " - "<< d_obs.rbegin()->t<<endl;
}

void KeplerLightCurve::plot(const std::string& fname)
{
  ofstream of(fname);
  of << std::fixed;
  for(const auto& o : d_obs) 
    of << o.t*1000.0 << '\t' << o.fixedFlux<<'\n';
}

void KeplerLightCurve::dump(const std::string& fname)
{
  ofstream r(fname);
  for(auto& o : d_obs) {
    r << o.t << '\t' << o.fixedFlux<<'\n';
  }
  r.flush();

  //  ofstream inter(fname+".inter");
  //  for(auto t = d_obs.begin()->t; t < d_obs.rbegin()->t; t+=0.5) {
  //    inter << t << '\t' << si(t) << '\n';
  //  }
}

void KeplerLightCurve::removeJumps()
{
  Observation* prev=0;
  double reflux=0, diff;
  unsigned int skips=0;
  for(auto& o : d_obs) {
    if(!prev) {
      prev = &o;
      continue;
    }
    diff = o.flux - prev->flux;
    if(fabs(diff) < 2000) {
      reflux+=diff;
    }
    else {
      skips++;
    }
    prev->fixedFlux=reflux;
    prev = &o;
  }
  prev->fixedFlux = reflux;
  cerr<<"Skipped "<<skips<<" implausible jumps"<<endl;

}

void KeplerLightCurve::removeDC()
{
  double tot=0;
  for(auto& o : d_obs) {
    tot += o.fixedFlux;
  }
  double average=tot/d_obs.size();
  cerr<<"Average: "<<average<<endl;

  for(auto& o : d_obs) {
    o.fixedFlux -= average;
  }
}
