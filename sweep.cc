#include "oscil.hh"
#include <fenv.h>
#include <fstream>

int main(int argc, char** argv)
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW); 

  typedef runge_kutta_cash_karp54< state_type > error_stepper_type;

  double ffreq=atoi(argv[1]);
  double step=0.1/ffreq;
  double q=atoi(argv[2]);

  double tau = 2*q/(2*M_PI*ffreq);
  cerr << "tau:"<< tau << endl;
  ofstream indiv("indiv");
  for(double freq  = 0.95*ffreq ; freq < 1.05*ffreq; freq+=0.25) {
    runge_kutta4< state_type > stepper;
    
    state_type x(2);
    x[0] = 0.0; // start at x=1.0, p=0.0
    x[1] = 0.0;
   
    FrequencyGenerator fg(freq, 200);
    HarmonicOscillatorFunctor<FrequencyGenerator> hof(ffreq, q, fg);

    integrate_adaptive(make_controlled<error_stepper_type>(1.0e-7, 1.0e-5), 
    		       hof, x, 0.0, 20*tau, step);

    Oscillator<double> o1(ffreq, q);
    /*
    for(double t = 0 ; t < 10; t+=step) {
      o1.feed(sin(2*M_PI*freq*t), step);
      if(!((counter++)%10) && fabs(ffreq-freq) < 3)
	indiv << freq << '\t' << t << '\t' << o1.get()<<'\n';
    }
    */
    cout << freq <<'\t'<<ffreq*ffreq*o1.get()/q<<'\t'<<
      sqrt(77.7895*ffreq*ffreq*(0.5*x[1]*x[1] + 0.5*x[0]*x[0]*hof.d_b)/(q*q))<<endl;
  }
}
