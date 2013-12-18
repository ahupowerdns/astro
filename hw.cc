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

class OscillatorBank
{
public:
  OscillatorBank(double startFreq, double endFreq, double gap);

  void feed(double val, double interval);
  vector<pair<double, double> > get();
private:
  double d_startFreq;
  double d_endFreq;
  double d_gap;
  vector<Oscillator<double>> d_bank;
};

OscillatorBank::OscillatorBank(double startFreq, double endFreq, double gap)
  : d_startFreq(startFreq), d_endFreq(endFreq), d_gap(gap)
{
  for(auto freq = startFreq; freq < endFreq; freq += gap) {
    d_bank.push_back({freq, 10*freq/gap});
  }
}

void OscillatorBank::feed(double val, double interval)
{
  for(auto& o : d_bank)
    o.feed(val, interval);
}

vector<pair<double, double>> OscillatorBank::get()
{
  vector<pair<double,double>> ret;
  for(auto& o : d_bank)
    ret.push_back({o.d_freq, o.get()});
  return ret;
}

int main(int argc, char**argv)
{   
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); 

  struct Observation
  {
    double jd;
    double rflux;
    double flux;
    double fixedFlux;
    bool operator<(const Observation& rhs) const 
    {
      return jd < rhs.jd;
    }
  };
  
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
      obs.push_back({jd[i], rflux[i], flux[i]});
    }

  }
  sort(obs.begin(), obs.end());
  cerr<< 86.4*obs.begin()->jd << " - "<< 86.4*obs.rbegin()->jd<<endl;

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

  OscillatorBank ob(0.09, 0.160, 0.00005);

  vector<vector<pair<double, double> > > results;
  
  double y;
  double beginTime = obs.begin()->jd*86.400;
  double endTime = obs.rbegin()->jd*86.400; // kiloseconds

  double step = 0.5; // 500 seconds
  double t = obs.begin()->jd * 86.400;
  auto o = obs.begin();
  for(unsigned int n = 0; t < endTime; t=beginTime + n++*step) {
    y=0;
    while(o != obs.end() && t <= o->jd*86.400  && o->jd*86.400 < t+step) {
      //      cout << 86.400*o->jd <<'\t'<<o->flux<<'\t'<<o->fixedFlux<<'\n';
      y+= o->fixedFlux;
      o++;
    } 

    //    y += 1000*(sin(2*M_PI*t*0.10)+sin(2*M_PI*t*0.11) + sin(2*M_PI*t*0.12));
    ob.feed(y, step);
  
    if(!(n%150)) {
      results.emplace_back(ob.get());
    }
  }

  
  
  double maxVal=0;
  vector<double> pixelhisto;
  for(auto& row : results) 
    for(auto& val: row) { 
      if(val.second < 7e8)
	val.second=0;
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
  
  plot << "P6"<<"\n"<<results.begin()->size()<<" "<<results.size()<<"\n255\n";
  
  for(const auto& row : results) {
    for(auto valpair: row) {
      auto val = valpair.second;
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
 
  reverse(results.begin(), results.end());
  
  vector<pair<double, double> > avg(results.begin()->size());
  for(unsigned int i =0 ; i < results.size()*0.9; ++i) {
    auto& row = results[i];
   
    for(unsigned int j = 0 ; j < row.size() ; ++j) {
      avg[j].second+=row[j].second;
      if(!i)
	avg[j].first=row[j].first; 
    }
  }
  for(unsigned int i = 0 ; i < avg.size(); ++i) {
    cout << avg[i].first << "\t";
    cout<<avg[i].second/results.size()<<endl;
  }
  

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
    cout << d << endl;
      

  return 0;       
}
