#ifndef OSCIL_HH
#define OSCIL_HH
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include "boost/numeric/odeint.hpp"

using namespace std;

template<class T> 
class Oscillator
{
public:
  typedef T value_type;

  T d_r;     // friction, newton/meter/second
  T d_b;     // restoring force, newton/meter
  //  T d_m;
  //  T d_rate, d_irate;  // configured feeding rate

  T d_x;
  T d_v;

  T d_freq;
  T d_brake;
  //  Averager d_pavg, d_vavg, d_xavg;
  unsigned int d_ticks;

public:
  Oscillator() 
  {}

  Oscillator(T freq, T q)
  {
    init(freq, q);
  }

  void init(T freq, T q)
  {
    /* Math:
       Linear oscillator with linear friction.
    */
    d_freq = freq;
    freq *= 2 * M_PI;
    //    d_m = 1.0;
    d_r = freq * 1.0 / q;                        // 1.0 was d_m
    d_b = freq * freq * 1.0 / (1-(0.5/(q*q)));   // 1.0 was d_m
    d_v = d_x =0;
    /*
    d_pavg.init(0.05, (unsigned int)freq);
    d_vavg.init(0.05, (unsigned int)freq);
    d_xavg.init(0.05, (unsigned int)freq);
    */
    d_brake=1;
  }
  
  T getFreq() const
  {
    return d_freq;
  }

  /*
  float power()
  {
    return d_pavg.get();
  }

  float msgFreq()
  {
    return d_vavg.max() / (6.28*d_xavg.max());
  }
  */

  T force(T x, T v, T input)
  {
    return -d_r * v - d_b * x + input;
  }

  void feed(T value, T interval)
  {
    
    /*
      d2x                    dx
      ----   =  value  - r * --  - b *x
      dt^2                   dt


    T force = 0;
    force -= d_r * d_v;  // friction
    force -= d_b * d_x;  // restoring force
    force += value;      // our input
    */

    /* Euler - for shame 
    const T F = force(d_x, d_v, value);
    d_v += F * interval;       // account speed    
    d_x += d_v * interval; */

    T k1, k2, k3, k4, l1, l2, l3, l4;
    //    const T F = force(d_x, d_v, value);

    k1 = interval * d_v;
    l1 = -interval * (d_r * d_v + d_b * d_x - value);

    k2 = interval * (d_v + 0.5*l1);
    l2 = -interval * (d_r * (d_v + 0.5*l1) + d_b * (d_x + 0.5*k1) - value);

    k3 = interval * (d_v + 0.5*l2);
    l3 = -interval * (d_r * (d_v + 0.5*l2) + d_b * (d_x + 0.5*k2) - value);

    k4 = interval * (d_v + l3);
    l4 = -interval * (d_r * (d_v + l3) + d_b * (d_x + k3) - value);

    d_x = d_x + k1/6 + k2/3 + k3/3 + k4/6;
    d_v = d_v + l1/6 + l2/3 + l3/3 + l4/6; 
    
  }

  /*
k1 = interval * d_v;
l1 = -interval * (d_r * d_v + d_b * d_x) ;

k2 = interval * (d_v + l1 / 2.0);
l2 = -interval * (d_r * (d_v + l1 / 2.0) + d_b * (d_x + k1 / 2.0)) ;

k3 = interval * (d_v + l2 / 2.0);
l3 = -interval * (d_r * (d_v + l2 / 2.0) + d_b * (d_x + k2 / 2.0)) ;

k4 = interval * (d_v + l3);
l4 = -interval * (d_r * (d_v + l3) + d_b * (d_x + k3)) ;
  */

  void reset()
  {
    d_v=d_x=0;
  }

  T get() const
  {
    return 0.5 * 1.0 * d_v * d_v + 0.5 * d_b * d_x * d_x;
  }
};


/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

struct FrequencyGenerator
{
  explicit FrequencyGenerator(double freq, double amplitude) : 
    d_omega(2*M_PI*freq),
    d_amplitude(amplitude)
  {}

  double operator()(double t) const
  {
    return d_amplitude*sin(d_omega*t);
  }
  
  double d_omega;
  double d_amplitude;
};

template<typename SignalGenerator>
class HarmonicOscillatorFunctor
{
public:
  HarmonicOscillatorFunctor(double freq, double q, const SignalGenerator& sg)
    : d_generator(&sg)
  {
    d_r = 2*M_PI*freq / q;  // oscillator frequency
    d_b = 2*M_PI*2*M_PI*freq * freq / (1-(0.5/(q*q)));
    d_freq = freq;
    d_q = q;
  }
  
  void operator()(const state_type& x, state_type &dxdt, const double& t) const
  {
    dxdt[0] = x[1];                          // x[0] = position
    dxdt[1] = -d_b*x[0] - d_r*x[1] + (*d_generator)(t);   // x[1] = speed, dxdt = force
  }


  const SignalGenerator* d_generator;
  double d_b;
  double d_r;
  double d_freq;
  double d_q;
};

using namespace boost::numeric::odeint;

class OscillatorBank
{
public:
  OscillatorBank(double startFreq, double endFreq, double gap);

  void feed(double val, double interval);
  std::vector<std::pair<double, double> > get();
private:
  double d_startFreq;
  double d_endFreq;
  double d_gap;
  std::vector<Oscillator<double>> d_bank;
};


#endif
