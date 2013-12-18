#include "oscil.hh"

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
    ret.push_back({o.d_freq, sqrt(o.get())});
  return ret;
}
