#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include "misc.hh"
#include <errno.h>
#include <string.h>
#include <boost/circular_buffer.hpp>
using namespace std;

struct Datum
{
  double t;
  double power;
};

int main(int argc, char**argv)
{
  if(argc==1) {
    cerr<<"Syntax: peaker perfreq/*"<<endl;
    return EXIT_FAILURE;
  }
  map<double, vector<Datum> > data;

  for(int n = 1; n < argc ; ++n) {
    string fname(argv[n]);
    auto pos = fname.find_last_of("/");
    if(pos == string::npos)
      continue;
    double freq = atof(fname.c_str()+pos+1);
    vector<Datum>& dest = data[freq];
    if(!data.empty())
      dest.reserve(data.begin()->second.size());

    FILE* fp=fopen(fname.c_str(), "r");
    if(!fp)
      throw runtime_error("Unable to open file "+fname+" for reading: "+string(strerror(errno)));

    ifstream ifs(fname);
    
    auto size=filesize(fname.c_str());
    dest.resize(size/sizeof(Datum));
    if(fread(&dest[0], 1, size, fp) != size) {
      throw runtime_error("Error reading data from '"+fname+"'");
    }
    fclose(fp);
  }
  cerr<<"Read data for "<<data.size() <<" frequencies, "<< data.begin()->second.size()*data.size()<<" elements"<<endl;

  ofstream ofs("power.average");
  int numrows = data.begin()->second.size();
  for(int n = 0 ; n < numrows; ++n) {
    VarMeanEstimator vme;
    for(auto& val : data) {
      vme(val.second[n].power);
    }
    ofs << data.begin()->second[n].t << '\t' << mean(vme)<<endl;
  }
  cerr<<"Done with averages, now unlikelies... "<<endl;

  ofstream unlikelies("unlikelies");
  for(auto& val : data) {
    double freq = val.first;
    vector<Datum>& d = val.second;

    boost::circular_buffer<double> ringbuf(12*30);   // 30 days
    double unlikely=0;
    VarMeanEstimator sig;
    for(const auto& pw : d) {
      sig(pw.power);
      ringbuf.push_back(pw.power);
      VarMeanEstimator vme;
      for(auto& val : ringbuf) {
	vme(val);
      }
      if(ringbuf.size() > 10) {
	double sigma = mean(vme)/sqrt(2*variance(vme));

	if(sigma > 4)
	  unlikely += 10;
	else if(sigma > 3)
	  unlikely+= 1;
	else if(sigma < 2)
	  unlikely -= 1;	
      }
    }
    unlikelies << freq << '\t' << unlikely<<'\t'<<mean(sig)<<'\t'<<sqrt(variance(sig))<<'\n';

  }
  cerr<<"Done"<<endl;
  
}

#if 0
  vector<double> forfft;
  forfft.reserve(column.size());
  for(const auto& ent : column)
   forfft.push_back(ent.second);

  int nc = (forfft.size()/2+1) ;
  fftw_complex *out = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * nc );
  
  pthread_mutex_lock(&g_lock); // seriously, they need this
  auto plan_forward = fftw_plan_dft_r2c_1d ( forfft.size(), &forfft[0], out, FFTW_ESTIMATE );
  pthread_mutex_unlock(&g_lock); // so sad
  fftw_execute ( plan_forward );
  vector<double> power;
  for(int n = 0 ; n < nc ; ++n) {
   power.push_back(sqrt(out[n][0]*out[n][0] + out[n][1]*out[n][1]));
  }
  fftw_free(out);
  normalize(power);
  ofstream powfile("perfreq/"+freqname+".power");
  double q[5]={0,0,0,0,0};
  for(unsigned int n =0 ; n < power.size() ; ++n) {
    powfile << n <<'\t' << power[n]<<'\n';
    if(n < 5)
      q[0]+=power[n];
    else
      if(n < 10)
	q[1]+=power[n];
      else if(n<50)
	q[2]+=power[n];
      else if(n<100)
	q[3]+=power[n];
      else
	q[4]+=power[n];
  }
  {   
    pthread_mutex_lock(&g_lock);

    op.freq = o.d_freq;
    op.unlikely = unlikely;
    op.mean = mean(totvme);
    op.stddev = sqrt(variance(totvme));
    g_postos.push_back(op);
    pthread_mutex_unlock(&g_lock); 
  }
  perfreq << " # unlikely score: "<<unlikely<<endl;
  perfreq << " # quartile powers ";
  for(auto n : {0,1,2,3,4})
    perfreq << q[n] << '\t';
  perfreq <<endl;
  
  sort(g_postos.begin(), g_postos.end(), 
       [](const OscillatorPerformance& a, const OscillatorPerformance& b) {
	 return a.freq < b.freq;
       }
       );
//////////////

  int winlen = g_postos.size()/20;
  ofstream unlikelies("unlikelies");
  for(auto iter = g_postos.begin(); iter != g_postos.end(); ++iter) {
    VarMeanEstimator meanVme, stddevVme;
    auto inner = iter - winlen/2;
    if(inner < g_postos.begin())
      inner = g_postos.begin();
    for(  ; inner < g_postos.end() && inner < iter + winlen/2 ; ++inner) {
      meanVme(inner->mean);
      stddevVme(inner->stddev);
    }
    iter->smoothedMean = mean(meanVme);
    unlikelies << iter->freq <<'\t' << iter->unlikely << '\t' << iter->mean << '\t' << iter->smoothedMean <<
      '\t' << iter->stddev  << '\t' << mean(stddevVme);
    for(auto& p : iter->powers) 
      unlikelies << '\t' << p;
    unlikelies << '\n';
  }
#endif
