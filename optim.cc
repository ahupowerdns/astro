#include <stdexcept>
#include <algorithm>
#include <string>
#include "misc.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <errno.h>
#include <boost/circular_buffer.hpp>
#include <boost/lexical_cast.hpp>
#include "nr3.h"
#include "amoeba.h"

using namespace std;

class MiniLineReader
{
public:
  MiniLineReader(const std::string& fname) 
  {
    d_fp=fopen(fname.c_str(), "r");
    if(!d_fp)
      throw runtime_error("Unable to open file '"+fname+"': "+string(strerror(errno)));
  }

  bool getLine(std::string* line)
  {
    return stringfgets(d_fp, line);
  }

  ~MiniLineReader()
  {
    fclose(d_fp);
  }
private:
  FILE* d_fp;
};

class TrueSource
{
public:
  TrueSource(const string& fname)
  {
    MiniLineReader mlr(fname);
    string line;
    double val;
    while(mlr.getLine(&line)) {
      if(line.empty() || line[0]=='#')
	continue;
      
      if((val=atof(line.c_str()))) {
	d_freqs.push_back({val/1000.0, false});
      }
    }
    sort(d_freqs.begin(), d_freqs.end());
    cerr<<"Have "<<d_freqs.size()<< " entries\n";
  }

  bool isPresent(double f, double delta, double* diff)
  {
    pair<double, bool> query(f, false);
    auto pos = lower_bound(d_freqs.begin(), d_freqs.end(), query);
    if(pos == d_freqs.end())
      return false;
    if(pos != d_freqs.begin())
      --pos;
    for(; pos < d_freqs.end(); ++pos) {
      //      cout<<"Comparing "<<*pos<<" and " << f << ": "<<fabs(f-*pos)<<endl;
      if(fabs(f-pos->first) < delta) {
	*diff = f - pos->first;
	pos->second=true;
	return true;
      }
      if(pos->first - f > delta)
	return false;
    }
    return false;
  }

  vector<double> getUnmatched()
  {
    vector<double> ret;
    for(const auto& f : d_freqs) 
      if(!f.second)
	ret.push_back(f.first);

    return ret;
  }
  
  unsigned int size()
  {
    return d_freqs.size();
  }

  void reset()
  {
    for(auto& f : d_freqs) 
      f.second=false; // unmatched
  }
private:
  vector<pair<double, bool> > d_freqs;
};

struct Possibility
{
  double freq;
  double unlikely;
  double mean;
  double smoothedMean;
  double stddev;
  double smoothedStddev;
  vector<double> powers;
};

vector<Possibility> filter(const vector<Possibility>& poss, double highsig, double sig, double nosig, double unl, double meanfact, double stdfact, double fact, double sigmonths)
{
  vector<Possibility> ret;
  for(const auto& p : poss) {
    boost::circular_buffer<double> ringbuf(sigmonths*30);  
    double unlikely=0;
    for(const auto& pw : p.powers) {
      ringbuf.push_back(pw);
      VarMeanEstimator vme;
      for(auto& val : ringbuf) {
	vme(val);
      }
      if(ringbuf.size() > 5) {
	double sigma = mean(vme)/sqrt(2*variance(vme));
       
	if(sigma > highsig)
	  unlikely += 3*fact;
	else if(sigma > sig)
	  unlikely+= fact;
	else if(sigma < nosig)
	  unlikely -= fact;	
      }
    }
    if(unlikely > 500*unl && p.mean/p.smoothedMean > meanfact && p.stddev/p.smoothedStddev > stdfact)
      ret.push_back(p);
  }
  return ret;
}

double score(TrueSource& ts, vector<Possibility>& cands, double* matches, double* mismatches, double* missed)
{
  *matches = *mismatches = *missed=0;
  double ret=0;
  for(const auto& c: cands) {
    double delta;
    if(ts.isPresent(c.freq, 0.00009, &delta)) {
      //      cout<<"Correct: "<<c.freq<< " (delta = "<<delta<<", actual="<<c.freq-delta<<")\n";
      // (*matches)++;
      // ret+=3;
    }
    else {
      //      cout<<"Wrong: "<<c.freq<<endl;
      if(c.freq > 0.1) {
	(*mismatches)++;
	ret-=35;
      }
    }
  }

  for(const auto& u : ts.getUnmatched()) {
       //    cout<<"Unmatched: "<<u<<endl;
    if(&u) { // silence warning
      (*missed)++;
      ret-=1;
    }
  }
  (*matches) = ts.size()-(*missed);
  ret+=4*(*matches);
  ts.reset();
  return -ret;
}

Doub func(VecDoub_I& v)
{
  return 0;
}

int main(int argc, char**argv)
{
  TrueSource ts(argv[1]);

  ifstream unl("unlikelies");
  vector<Possibility> poss;
  string line;
  while(getline(unl, line)) {
    istringstream istr(line);
    Possibility p;
    istr >> p.freq;
    istr >> p.unlikely;
    istr >> p.mean;
    istr >> p.smoothedMean;
    istr >> p.stddev;
    istr >> p.smoothedStddev;
    double pwr;
    while(istr >> pwr) 
      p.powers.push_back(pwr);
    poss.push_back(p);
  }
  cerr<<"Have "<<poss.size()<<" possibilities"<<endl;

  struct func
  {
    func(TrueSource& ts, vector<Possibility>& poss) 
      : d_ts(ts), d_poss(poss), d_output(false), d_count(0)
    {}
    
    double operator()(VecDoub_I& v)
    {
      auto cands = filter(d_poss, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
      double matches, mismatches, missed;
      auto s = score(d_ts, cands, &matches, &mismatches, &missed);
      cerr<<"Score "<<d_count<<": "<<s<< ", candidates: "<<cands.size()<<" matched: "<<matches<<", MISMATCHES: "<<mismatches<<", missed: "<<missed<<" @ highsig=" << 
	v[0]<<", sig="<<v[1]<<", nosig="<<v[2]<<", unl="<<500*v[3]<<", avgrat="<<v[4]<<", stdrat="<<v[5]<<", fact="<<v[6]<<", sigdays="<<30*v[7]<<endl;

      ofstream fit("fit."+boost::lexical_cast<string>(d_count++));
     
      for(auto& c: cands) {
	fit<<c.freq<<"\t"<<c.mean<<'\n';
      }

      return s;
    }
    TrueSource& d_ts;
    vector<Possibility>& d_poss;
    bool d_output;
    int d_count;
  };

  func f(ts, poss);

  Amoeba a(0.001);
  VecDoub point(8);
  point[0]=5; // highsig
  point[1]=4; // sig
  point[2]=2; // nosig
  point[3]=1; // unlikely threshold unl (/1000)
  point[4]=1; // avgrat (meanrat)
  point[5]=1; // stdrat
  point[6]=4; // fact points
  point[7]=2; // number of months

  auto v = a.minimize(point, 1, f);
  cerr<<"Evaluations: "<<a.nfunc<<endl;
  f.d_output=true;
  f(v);
  

}
