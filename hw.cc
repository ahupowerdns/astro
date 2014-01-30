#include <CCfits/CCfits>
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
#include "misc.hh"

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

struct SignalInterpolator
{
  explicit SignalInterpolator(const vector<Observation>& obs)
    : d_obs(&obs)
  {
  }

  double operator()(double t) const {
    auto b=lower_bound(d_obs->begin(), d_obs->end(), t);
    if(b == d_obs->begin() || b == d_obs->end())
      return 0;
    
    auto a = b - 1;
    if(a == d_obs->begin())
      return 0;

    double dt = (b->t - a-> t);
    if(dt > 10) // makes little sense to interpolate over days
      return 0;
    double frac = (t - a->t) / (b->t - a-> t);
    double val =  frac * b->fixedFlux + (1-frac)*a->fixedFlux;
    //    cerr<<"t = "<<t<<", Frac: "<<frac<<", aVal: "<<a->fixedFlux<<", bVal: "<<b->fixedFlux<<", returning: "<<val<<endl;
    return val;
  }

  const vector<Observation>* d_obs;
};


struct CSplineSignalInterpolator
{
  explicit CSplineSignalInterpolator(const vector<Observation>& obs)
    : d_obs(&obs)
  {
    double p, qn, sig, un;

    int n=obs.size();
    vector<double> u(n-1);
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

  const vector<Observation>* d_obs;
  vector<double> d_y2;
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
  vector<Observation> d_obs;
};

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


void normalize(vector<double>& values)
{
  double tot = 0;
  for(const auto& val : values)
    tot+=val;
  const double average = tot/values.size();

  for(auto& val : values)
    val /= average;  
}

/* we have time values and for each time value numbers per frequency.
   We can promise we'll only add new time values in ascending order.
*/

class TwoDLabelledArray
{
private:
  struct Row
  {
    double t;
    vector<double> values;
  };
  vector<Row> d_rows;
};

typedef HarmonicOscillatorFunctor<CSplineSignalInterpolator> oscil_t;
pthread_mutex_t g_lock = PTHREAD_MUTEX_INITIALIZER;
vector<pair<double, double> >  doFreq(oscil_t& o, const vector<double>& otimes)
{
  cerr<<"freq: "<<o.d_freq<<endl;
  
  ofstream perfreq("perfreq/"+boost::lexical_cast<string>(o.d_freq));
  typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
  
  state_type x(2);
  x[0] = 0.0; // start at x=0.0, p=0.0
  x[1] = 0.0;
  
  vector<pair<double, double> > column;
  boost::circular_buffer<double> ringbuf(150);  
  double unlikely=1;
  integrate_times(make_controlled<error_stepper_type>(1e-1, 1.0e-8), 
		  o, x, otimes.begin(), otimes.end(), 0.1,
		  [&](const state_type& s, double t) {
		    double power = o.d_freq*o.d_freq*(0.5*x[1]*x[1] + 0.5*x[0]*x[0]*o.d_b)/(o.d_q*o.d_q);
		    ringbuf.push_back(power);

		    column.emplace_back(make_pair(
						  o.d_freq, 
						  power
						  ));
		    

		    VarMeanEstimator vme;
		    for(auto& val : ringbuf) {
		      vme(val);
		    }

		    perfreq << t << '\t' << power << '\t';
		    if(ringbuf.size() > 5) {
		      double sigma = mean(vme)/sqrt(2*variance(vme));
		      perfreq << mean(vme) << '\t' << variance(vme) << '\t' << sigma << endl;
		      if(sigma > 8)
			unlikely += 4;
		      else if(sigma > 4)
			unlikely+= 2;
		      else if(sigma < 2)
			unlikely -= 1;
		      else if(sigma < 1)
			unlikely -= 2;

		    }
		    else
		      perfreq << "0\t0\t0"<<endl;
		    perfreq.flush();		    
		  });
  
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
  ofstream powfile("perfreq/"+boost::lexical_cast<string>(o.d_freq)+".power");
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
    ofstream ofs("unlikelies", std::fstream::ate | std::fstream::app);
    ofs<<o.d_freq << '\t' << unlikely<<endl;
  }
  perfreq << " # unlikely score: "<<unlikely<<endl;
  perfreq << " # quartile powers ";
  for(auto n : {0,1,2,3,4})
    perfreq << q[n] << '\t';
  perfreq <<endl;
  
  return column;
}


vector<pair<double, double> > emitPowerGraph(double start, double stop, int index,  const vector<vector<pair<double, double> > >& results)
{
  vector<pair<double, double> > ret;
  ofstream pplot("power."+boost::lexical_cast<string>(index));
  vector<pair<double, double> > avg(results.size());

  unsigned int count=0;
  for(unsigned int i =0 ; i < results.size(); ++i) {
    const auto& column = results[i];
      
    for(unsigned int j = start*column.size() ; j < stop*column.size() && j < column.size() ; ++j) {
      avg[i].second+=column[j].second;
      avg[i].first=column[j].first; 
      count++;
    }
    
  }
    
  for(unsigned int i = 0 ; i < avg.size(); ++i) {
    pplot << avg[i].first << "\t";
    pplot<<avg[i].second/count <<endl;
    ret.push_back({avg[i].first,  avg[i].second/count});
  }
  return ret;
}

template<typename T>
vector<pair<typename T::iterator, typename T::iterator> >
splitRange(T& container, int parts)
{
  typedef pair<typename T::iterator, typename T::iterator> range_t;
  vector<range_t > ret;
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

int main(int argc, char**argv)
{   
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW); 
  KeplerLightCurve klc;
  unlink("unlikelies");
  for(int fnum = 1; fnum < argc; ++fnum) {
    if(boost::ends_with(argv[fnum], ".txt"))
      klc.addTxt(argv[fnum]);
    else
      klc.addFits(argv[fnum]);
  }
  klc.sort();

  klc.removeJumps();
  klc.removeDC();
  klc.plot("fits.plot");
  vector<vector<pair<double, double> > > results;

  CSplineSignalInterpolator si(klc.d_obs);
#if 0
  ofstream r("real");
  for(auto& o : klc.d_obs) {
    r << o.t << '\t' << o.fixedFlux<<'\n';
  }
  r.flush();

  ofstream inter("inter");
  for(auto t = klc.d_obs.begin()->t; t < klc.d_obs.rbegin()->t; t+=0.5) {
    inter << t << '\t' << si(t) << '\n';
  }
 #endif


  vector<double> otimes;
  for(double t = klc.d_obs.begin()->t ; t < klc.d_obs.rbegin()->t; t += 86) {
    otimes.push_back(t);
  }
  otimes.push_back(klc.d_obs.rbegin()->t);
#if 0
  double f=0.15791000000000191;
  oscil_t o(f, 20.0*f/0.00005, si);
  doFreq(o, otimes);
  return 0;
#endif
  bool insDone=false;
  vector<oscil_t> oscillators;
  for(double f = 0.09; f< 0.18; f+=0.00001) {
    oscillators.push_back({f, 10*f/0.0001, si});
  }

  
  auto doPart=[&](vector<oscil_t>::iterator begin, vector<oscil_t>::iterator end) {
    vector<vector<pair<double, double> > > ret;
    for(auto iter= begin; iter != end; ++iter) {
      auto column = doFreq(*iter, otimes);
      ret.emplace_back(column);
    }
    return ret;
  };

  vector<decltype(std::async(std::launch::async, doPart, oscillators.begin(), oscillators.end()))> futures;

  //  doPart(oscillators.begin(), oscillators.end());

  for(auto& s : splitRange(oscillators, 4)) {
    futures.emplace_back(std::async(std::launch::async, doPart, s.first, s.second));
  }

  
  for(auto& future : futures) {
    auto part=future.get();
    for(auto& p: part)
      results.push_back(p);
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

  
  for(int n=1; n <11; ++n) 
    emitPowerGraph((n-1)*0.1, n*0.1, n, results);
  emitPowerGraph(0.0, 1.0, 0, results);
  emitPowerGraph(0.1, 1.0, 100, results);
  emitPowerGraph(0.1, 1.0, 10100, results);
  emitPowerGraph(0.2, 1.0, 20100, results);
  auto graph = emitPowerGraph(0.4, 1.0, 40100, results);
  sort(graph.begin(), graph.end(), 
       [](const pair<double, double>& a, const pair<double, double>& b) 
       {
	 return a.second < b.second;
       }
       );

  ofstream peaks("peaks");
  for(auto peak : graph) {
    peaks << peak.first<<'\t' << peak.second<< '\n';
  }
  return 0;       
}



#if 0

  OscillatorBank ob(0.09, 0.160, 0.00002);

  vector<vector<pair<double, double> > > results;
  
  double y;
  double beginTime = klc.d_obs.begin()->t/1000;
  double endTime = klc.d_obs.rbegin()->t/1000; // kiloseconds

  double step = 0.1; // 500 seconds
  double t = klc.d_obs.begin()->t;
  auto o = klc.d_obs.begin();
  for(unsigned int n = 0; t < endTime; t=beginTime + n++*step) {
    y=0;
    while(o != klc.d_obs.end() && t <= o->t  && o->t < t+step) {
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
