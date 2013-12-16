#ifndef OSCIL_HH
#define OSCIL_HH
#include <cmath>
#include <iostream>
#include <utility>
using namespace std;

template<class T> 
class Oscillator
{
public:
  typedef T value_type;

  T d_r;     // friction, newton/meter/second
  T d_b;     // restoring force, newton/meter
  //  T d_m;
  T d_rate, d_irate;  // configured feeding rate

  T d_x;
  T d_v;

  T d_freq;
  T d_brake;
  //  Averager d_pavg, d_vavg, d_xavg;
  unsigned int d_ticks;

public:
  Oscillator() 
  {}

  Oscillator(T freq, T q, T rate)
  {
    init(freq, q, rate);
  }

  void init(T freq, T q, T rate)
  {
    /* Math:
       Linear oscillator with linear friction.
    */
    d_freq = freq;
    freq *= 2 * M_PI;
    //    d_m = 1.0;
    d_r = freq * 1.0 / q;                        // 1.0 was d_m
    d_b = freq * freq * 1.0 / (1-(0.5/(q*q)));   // 1.0 was d_m
    d_rate = rate;
    d_irate = 1/rate;
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
  inline void feed(T value, T interval)
  {
    /*
    T force = 0;
    force -= d_r * d_v;  // friction
    force -= d_b * d_x;  // restoring force
    force += value;      // our input
    */

    const T force = -d_r * d_brake * d_v - d_b * d_x + value;

    d_v += force * d_irate;       // account speed    

    d_x += d_v * d_irate;

    /*
    d_pavg.feed(get());

    d_vavg.feed(fabs(d_v));
    d_xavg.feed(fabs(d_x));
    */

  }

  void reset()
  {
    d_v=d_x=0;
    //    d_pavg.reset();
//    amp=0;
  }

  T get() const
  {
    return 0.5 * 1.0 * d_v * d_v + 0.5 * d_b * d_x * d_x;
  }

  T getL() const
  {
    return sqrt(0.5 * 1.0 * d_v * d_v + 0.5 * d_b * d_x * d_x);
  }
};
#endif
