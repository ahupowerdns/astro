#include "oscil.hh"
#include <fenv.h>

int main(int argc, char** argv)
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW); 

  double step=0.00011;
  OscillatorBank ob(100, 200, 1);

  vector<double> freqs;
  for(int n=1; n < argc; ++n)
    freqs.push_back(atoi(argv[n]));

  double force;
  for(double t = 0 ; t < 1; t+=step) {
    force=0;
    for(auto f : freqs)
      force += exp(-t/3.0)*sin(2*M_PI*t*f);
    ob.feed(force, step);
  }
  auto res = ob.get();
  for(auto& v: res) {
    cout<<v.first<<'\t'<<v.second<<'\n';
  }
}
