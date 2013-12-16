#include <CCfits/CCfits>
#include <string>
#include <iostream>

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
  const Observation* prev=0;
  double reflux=0, diff;
  unsigned int skips=0;
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
    }
    else
      cout << 0 << "\t0\t";
    cout<<reflux<<endl;
    prev = &o;
  }
  cerr<<"Skipped "<<skips<<" implausible jumps"<<endl;
  return 0;       
}
