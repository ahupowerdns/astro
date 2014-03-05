#include <boost/lexical_cast.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include "oscil.hh"
#include <fenv.h> 
#include <future>
#include <utility>
#include <fftw3.h>
#include <boost/algorithm/string.hpp>
#include <boost/circular_buffer.hpp>
#include "klc.hh"
#include "misc.hh"

using namespace std;

typedef HarmonicOscillatorFunctor<CSplineSignalInterpolator> oscil_t;

struct Status 
{
  double t;
  double frequency;
  double power;
};

vector<Status> doFreq(oscil_t& o, const vector<double>& otimes)
{
  cerr<<"freq: "<<o.d_freq<<endl;
  typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
  
  state_type x(2);
  x[0] = 0.0; // start at x=0.0, p=0.0
  x[1] = 0.0;
  
  vector<Status> column;

  integrate_times(make_controlled<error_stepper_type>(1e-1, 1.0e-8), 
		  o, x, otimes.begin(), otimes.end(), 0.1,
		  [&](const state_type& s, double t) {
		    double power = o.d_freq*o.d_freq*(0.5*x[1]*x[1] + 0.5*x[0]*x[0]*o.d_b)/(o.d_q*o.d_q);
		    column.push_back({t, o.d_freq, power});
		  });
  
  return column;
}


void emitPowerGraph(double start, double stop, int index,  const vector<vector<Status> >& results)
{
  vector<pair<double, double> > ret;
  ofstream pplot("power."+boost::lexical_cast<string>(index));
  vector<pair<double, double> > avg(results.size());

  unsigned int count=0;
  for(unsigned int i =0 ; i < results.size(); ++i) {
    const auto& column = results[i];
      
    for(unsigned int j = start*column.size() ; j < stop*column.size() && j < column.size() ; ++j) {
      avg[i].first=column[j].frequency; 
      avg[i].second+=column[j].power;
      count++;
    }    
  }
    
  for(unsigned int i = 0 ; i < avg.size(); ++i) {
    pplot << avg[i].first << "\t";
    pplot<<avg[i].second/count <<endl;
    //    ret.push_back({avg[i].first,  avg[i].second/count});
  }
  //  return ret;
}


int main(int argc, char**argv)
{   
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW); 
  KeplerLightCurve klc;

  if(argc==1) {
    cerr<<"Syntax: hw 1.fits|1.txt [2.fits|2.txt] ..."<<endl;
    return EXIT_FAILURE;
  }

  for(int fnum = 1; fnum < argc; ++fnum) {
    if(boost::ends_with(argv[fnum], ".txt"))
      klc.addTxt(argv[fnum]);
    else
      klc.addFits(argv[fnum]);
  }

  klc.sort();
  //  klc.reverse();

  klc.removeJumps();
  klc.removeDC();
  klc.plot("lightcurve.plot");
  CSplineSignalInterpolator si(klc.d_obs); // provides continuous view on our sampled data

  // these are the times we want results for
  vector<double> otimes;
  for(double t = klc.d_obs.begin()->t ; t < klc.d_obs.rbegin()->t; t += 7.2) {
    otimes.push_back(t);
  }
  otimes.push_back(klc.d_obs.rbegin()->t);

  // these are the frequencies we want results for
  vector<oscil_t> oscillators;
  for(double f = 0.09; f< 0.18; f+=0.00001) {
    oscillators.push_back({f, 10*f/0.0001, si});
  }

  // this is a single thread that works for us
  auto doPart=[&](vector<oscil_t>::iterator begin, vector<oscil_t>::iterator end) {
    vector<vector<Status > > ret;
    for(auto iter= begin; iter != end; ++iter) {
      auto column = doFreq(*iter, otimes);
      ret.emplace_back(column);
    }
    return ret;
  };

  // launch 4 threads - astoundingly ugly
  vector<decltype(std::async(std::launch::async, doPart, oscillators.begin(), oscillators.end()))> futures;
  for(auto& s : splitRange(oscillators, 4)) {
    futures.emplace_back(std::async(std::launch::async, doPart, s.first, s.second));
  }
  vector<vector<Status> > results;
  
  // reap the results
  for(auto& future : futures) {
    auto part=future.get();
    for(auto& p: part)
      results.push_back(p);
  }

  sort(results.begin(), results.end(), [](const vector<Status>& a, const vector<Status>&b) 
       {
	 return a.cbegin()->frequency < b.cbegin()->frequency;
       });


  for(auto& column : results) {
    string filename = "perfreq/"+ (boost::format("%f") % column.begin()->frequency).str();
    FILE* fp=fopen(filename.c_str(), "w");
    if(!fp)
      throw runtime_error("Unable to open file "+filename+": "+string(strerror(errno)));
    for(auto& val: column) { 
      fwrite(&val.t, sizeof(val.t), 1, fp);
      fwrite(&val.power, sizeof(val.power), 1, fp);
    }
    fclose(fp);
  }
      

  double maxVal=0;
  vector<double> pixelhisto;
  for(auto& column : results) 
    for(auto& val: column) { 
      maxVal=std::max(maxVal, val.power*val.power);
      pixelhisto.push_back(val.power);
    }

  sort(pixelhisto.begin(), pixelhisto.end());
  ofstream pixelhistofile("pixelhisto");
  double limitstep = *pixelhisto.rbegin()/1000;
  double nowlimit=limitstep;
  unsigned int limitcount=0;
  for(auto d : pixelhisto) {
    limitcount++;
    if(d > nowlimit) {
      pixelhistofile<<nowlimit<<'\t'<<limitcount<<'\n';
      limitcount=0;
      nowlimit += limitstep;
    }
  }
  pixelhistofile<<nowlimit<<'\t'<<limitcount<<'\n';
  pixelhistofile.flush();
  ofstream plot("waterfall.ppm");
  
  plot << "P6"<<"\n"<<results.size()<<" "<<results.begin()->size()<<"\n255\n";

  for(unsigned y=0; y < results.begin()->size(); ++y) {  
    for(unsigned x=0; x < results.size(); ++x) {
      auto val = results[x][y].power;
      double r = 2*val*val/maxVal;
      double g = 2*val*val/maxVal - 0.5;
      double b=  2*val*val/maxVal  -1;
      if(r < 0) r=0;
      if(r > 1) r=1;
      if(g < 0) g=0;
      if(g > 1) g=1;
      if(b < 0) b=0;
      if(b > 1) b=1;

      plot << (char)(255*r);
      plot << (char)(255*g);
      plot << (char)(255*b);
    }
  }

  
  for(int n=1; n <11; ++n) 
    emitPowerGraph((n-1)*0.1, n*0.1, n, results);
  emitPowerGraph(0.0, 1.0, 0, results);
  emitPowerGraph(0.1, 1.0, 100, results);
  emitPowerGraph(0.1, 1.0, 10100, results);
  emitPowerGraph(0.2, 1.0, 20100, results);
  emitPowerGraph(0.4, 1.0, 40100, results);
  return 0;       
}
