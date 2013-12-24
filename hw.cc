#include <CCfits/CCfits>
#include <boost/lexical_cast.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include "oscil.hh"
#include <fenv.h> 
#include <utility>

using namespace std;
using namespace CCfits;

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

struct SignalInterpollator
{
  explicit SignalInterpollator(const vector<Observation>& obs)
    : d_obs(obs)
  {
  }

  double operator()(double t) const {
    auto b=lower_bound(d_obs.begin(), d_obs.end(), t);
    if(b == d_obs.begin() || b == d_obs.end())
      return 0;
    
    auto a = b - 1;
    if(a == d_obs.begin())
      return 0;

    double dt = (b->t - a-> t);
    if(dt > 10) // makes little sense to interpollate over days
      return 0;
    double frac = (t - a->t) / (b->t - a-> t);
    double val =  frac * b->fixedFlux + (1-frac)*a->fixedFlux;
    //    cerr<<"t = "<<t<<", Frac: "<<frac<<", aVal: "<<a->fixedFlux<<", bVal: "<<b->fixedFlux<<", returning: "<<val<<endl;
    return val;
  }

  const vector<Observation>& d_obs;
};

int main(int argc, char**argv)
{   
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW); 

  
  vector<Observation> obs;

  for(int fnum = 1; fnum < argc; ++fnum) {
    const char* fname=argv[fnum];
    cerr<<fname<<endl;
    CCfits::FITS pInfile(fname, Read); 

    // define a reference for clarity. (std::auto_ptr<T>::get returns a pointer

    ExtHDU& table = pInfile.extension(1);
    
    // read all the keywords, excluding those associated with columns.
    
    table.readAllKeys();
    
    // print the result.
    
    //  std::cout << table << std::endl;
    std::vector <double> flux;
    table.column("PDCSAP_FLUX").read(flux, 0, table.rows());

    std::vector <double> rflux;
    table.column("SAP_FLUX").read(rflux, 0, table.rows());


    std::vector <double> jd;
    table.column("TIME").read(jd, 0, table.rows());
    
    for(unsigned int i=0; i < jd.size(); ++i) {
      if(std::isnan(jd[i]) || std::isnan(rflux[i]) || std::isnan(flux[i])) 
	continue;
      obs.push_back({86.0*jd[i], rflux[i], flux[i]});
    }

  }
  sort(obs.begin(), obs.end());
  cerr<< "Kilosecond timespan: "<<obs.begin()->t << " - "<< obs.rbegin()->t<<endl;

  Observation* prev=0;
  double reflux=0, diff;
  unsigned int skips=0;
  for(auto& o : obs) {
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

  double tot=0;
  for(auto& o : obs) {
    tot += o.fixedFlux;
  }
  double average=tot/obs.size();
  cerr<<"Average: "<<average<<endl;

  for(auto& o : obs) {
    o.fixedFlux -= average;
  }

  vector<vector<pair<double, double> > > results;

  SignalInterpollator si(obs);
  vector<HarmonicOscillatorFunctor<SignalInterpollator>> oscillators;
  for(double f = 0.09; f< 0.160; f+=0.00005) {
    oscillators.push_back({f, 10*f/0.00005, si});
  }
  
  vector<double> otimes;
  for(double t = obs.begin()->t ; t < obs.rbegin()->t; t += 86) {
    otimes.push_back(t);
  }
  otimes.push_back(obs.rbegin()->t);

  for(auto& o : oscillators) {
    cerr<<"freq: "<<o.d_freq<<endl;

    ofstream perfreq("perfreq/"+boost::lexical_cast<string>(o.d_freq));
    typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
    
    state_type x(2);
    x[0] = 0.0; // start at x=0.0, p=0.0
    x[1] = 0.0;

    vector<pair<double, double> > column;
    
    integrate_times(make_controlled<error_stepper_type>(1.0e-7, 1.0e-5), 
		    o, x, otimes.begin(), otimes.end(), 0.1,
		    [&](const state_type& s, double t) {
		      column.emplace_back(make_pair(
						 o.d_freq, 
						 o.d_freq*o.d_freq*(0.5*x[1]*x[1] + 0.5*x[0]*x[0]*o.d_b)/(o.d_q*o.d_q)
						 ));
		      
		      perfreq << t << '\t' << o.d_freq*o.d_freq*(0.5*x[1]*x[1] + 0.5*x[0]*x[0]*o.d_b)/(o.d_q*o.d_q) << endl;
		      perfreq.flush();
		      
		    });

    
    results.emplace_back(column);
  }

  double maxVal=0;
  vector<double> pixelhisto;
  for(auto& column : results) 
    for(auto& val: column) { 
      maxVal=std::max(maxVal, val.second*val.second);
      pixelhisto.push_back(val.second);
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
  ofstream plot("plot");
  
  plot << "P6"<<"\n"<<results.size()<<" "<<results.begin()->size()<<"\n255\n";

  for(unsigned y=0; y < results.begin()->size(); ++y) {  
    for(unsigned x=0; x < results.size(); ++x) {
      auto val = results[x][y].second;
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
 
  vector<pair<double, double> > avg(results.size());
  for(unsigned int i =0 ; i < results.size(); ++i) {
    auto& column = results[i];
    
    for(unsigned int j = 0.1*column.size() ; j < column.size() ; ++j) {
      avg[i].second+=column[j].second;
      avg[i].first=column[j].first; 
    }
  }
  for(unsigned int i = 0 ; i < avg.size(); ++i) {
    cout << avg[i].first << "\t";
    cout<<avg[i].second/results.size()<<endl;
  }
  

  ofstream td("topdistances");
  sort(avg.begin(), avg.end(), [](const pair<double, double>& first, 
				  const pair<double, double>& second) {
	 return first.second < second.second;
       });

  vector<double> topFreqs;
  unsigned int i=0;
  for(auto iter = avg.rbegin(); iter != avg.rend() && i< 20; ++i, ++iter) 
    topFreqs.push_back(iter->first);

  vector<double> distances;
  for(auto a : topFreqs) {
    for(auto b : topFreqs) {
      if(a<b)
	distances.push_back(b-a);
    }
  }
  sort(distances.begin(), distances.end());
  for(auto d : distances)
    td << d << endl;      

  return 0;       
}



#if 0

  OscillatorBank ob(0.09, 0.160, 0.00002);

  vector<vector<pair<double, double> > > results;
  
  double y;
  double beginTime = obs.begin()->t/1000;
  double endTime = obs.rbegin()->t/1000; // kiloseconds

  double step = 0.1; // 500 seconds
  double t = obs.begin()->t;
  auto o = obs.begin();
  for(unsigned int n = 0; t < endTime; t=beginTime + n++*step) {
    y=0;
    while(o != obs.end() && t <= o->t  && o->t < t+step) {
      //      cout << o->t/1000.0 <<'\t'<<o->flux<<'\t'<<o->fixedFlux<<'\n';
      y+= o->fixedFlux;
      o++;
    } 

    //    y += 1000*(sin(2*M_PI*t*0.10)+sin(2*M_PI*t*0.11) + sin(2*M_PI*t*0.12));
    ob.feed(5*y, step);
  
    if(!(n%150)) {
      results.emplace_back(ob.get());
    }
  }

#endif
