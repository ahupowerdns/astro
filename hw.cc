#include <CCfits/CCfits>
#include <boost/lexical_cast.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include "oscil.hh"

using namespace std;
using namespace CCfits;

int main(int argc, char**argv)
{        
  struct Observation
  {
    double jd;
    double rflux;
    double flux;

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
  for(unsigned int tries=0; tries < 10; ++tries) {
    vector<Oscillator<double>> bank;
    for(double f=0.0900 + 0.0001*tries; f < 0.160; f+=0.00005)
      bank.push_back({f, 20000, 1000.0/(1765.408765)});
    
    vector<vector<double>> results;
    
    const Observation* prev=0;
    double reflux=0, diff;
    unsigned int skips=0;
    ofstream plot("plot"+boost::lexical_cast<string>(tries));
    unsigned int count=0;
    for(const auto& o : obs) {
      cout << std::fixed << 86.400*o.jd <<"\t"<<o.flux<<"\t";
      if(prev) {
	diff = prev->flux - o.flux;
	if(fabs(diff) < 2000) {
	  reflux+=diff;
	}
	else {
	  skips++;
	}
	
	cout << diff << '\t' << (prev->jd - o.jd)*86400<<'\t';
	for(auto& osc : bank)
	  osc.feed(reflux, (prev->jd - o.jd)*86.400);
	
      }
      else
	cout << 0 << "\t0\t";
      cout<<reflux<<endl;
      prev = &o;
      
      if(!(count%10)) {
	vector<double> row;
	for(auto& osc : bank) {
	  row.push_back(osc.get());
	}
	results.emplace_back(row);
      }
      count++;
    }
    cerr<<"Skipped "<<skips<<" implausible jumps"<<endl;
    
    double maxVal=0;
    for(const auto& row : results) 
      for(const auto& val: row) 
	maxVal=std::max(maxVal, val);
    
    plot << "P3"<<"\n"<<" "<<results.begin()->size()<<" "<<results.size()<<"\n255\n";
    
    for(const auto& row : results) {
      for(auto val: row) {
	unsigned int color=255*val/maxVal;
	if(color > 255)
	  color=255;
	plot << (color) <<' ';
	plot << (color) <<' ';
	plot << (color) <<' ';
	plot<<' ';
      }
      plot<<endl;
    }
  }
  /*
  reverse(results.begin(), results.end());

  vector<double> avg(results.begin()->size());
  for(unsigned int i =0 ; i < results.size()*0.9; ++i) {
    auto& row = results[i];
    for(unsigned int j = 0 ; j < row.size() ; ++j)
      avg[j]+=row[j];
  }
  for(unsigned int i = 0 ; i < avg.size(); ++i) {
    cout << i << "\t";
    cout<<avg[i]/results.size()<<endl;
  }
  */

  return 0;       
}
