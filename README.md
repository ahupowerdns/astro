## astro
Hi! This code takes FITS Kepler data, and feeds the electron counts into a
simulated oscillator bank. These oscillators are tuned to the frequency
range we expect the electron counts to vary with.

Some simple pre-processing is performed to remove discontinuities from the
data, but otherwise the data used is 'raw'.

No matter what the frequency of input data is, for each oscillator we output
per day the energy it contains. If a powerful frequency component is
present, the oscillator will contain more energy.

The behaviour of each oscillator is output to 'perfreq/frequency'. The
fields are:
	
 * number of kiloseconds since start
 * current energy
 * averaged energy (150 days surrounding window, -75 - +75)
 * variance of energy (150 days)
 * sigmas distance from 0 based on 150 day data

There is also a '.power' file which contains a fourier spectrum of the
energy timeseries of the oscillator.

The tool can be run as:

	$ ./hw kplr009145955-20*

Or alternatively:

	$ ./hw allBert/pollux_ts.txt

NOTE: hw currently scans from 0.09 milliHertz to 0.18 milliHertz, to change,
edit this line in hw.cc:

	for(double f = 0.09; f< 0.18; f+=0.00001) {


## Further output
A PPM graph is generated called 'plot' which can be opened with many image
viewers (gimp works, chrome too). This shows a waterfall, with t=0 on top of
the image and the end of the data below.

In addition, all data is saved to a file called 'unlikelies'.

This file has one (very long) line per oscillator:

 * Frequency in millihertz
 * Static 'unlikely' score
 * Mean of signal strength over all observations
 * Mean of signal strength of 'surrounding' oscillators
 * Mean of standard deviation of signals strength of 'surrounding' oscillators
 * All recorded signal strengths (1 entry per day)

The static 'unlikely' score is the result of a built-in model, based on how
often (and how strongly) the signal contradicts a null-hypothesis that there
is no signal.

## The model
Each oscillator is studied individually. The challenge is to distinguish
between oscillators excited by random noise and those excited by actual
stellar modes. To do this, several measures are available. One of these is
the 'unlikely score'.

### The unlikely score
The daily energy values of the oscillator are averaged over a sliding
window, typically several months long.  Over this window, we calculate the
mean and standard deviation of the energy. We then calculate how many
sigmas this signal is away from 0. This can be seen as a sort of
'self-significance test'. 

For each day, if the 'sigma' is more than a high value, we award several
points. If it is higher than a medium value, we award less points. Finally,
if significance is below a certain value, we deduct points.

### Comparison with neighbouring oscillators
A window covering 20% of all oscillators, 10% with lower frequencies, 10%
with higher is used to determine the mean and standard deviation of the
energy of 'surrounding' oscillators.

We then compare the mean energy of an oscillator to that if the neighbours.
We do the same to the standard deviation of the oscillator.

### The decision model
We declare a mode to be real if:

 * The unlikely score is higher than a threshold (unl)
 * The energy is more than 'meanfact' times higher than the surroundings
 * The standard deviation of the energy is more than 'stdfact' times higher
   than the surroundings
 
In code:

    if(unlikely > unl && p.mean/p.smoothedMean > meanfact && p.stddev/p.smoothedStddev > stdfact)

## All model parameters
Some typical values for the parameters are:

 * unl = 0
 * meanfact = 4
 * stdfact = 3

However, there are further parameters hiding in the unlikely score:
 * highsig, for which sigmas are considered truely significant
   (leading to 3*fact points)
 * sig, which sigma is considered significant (leading to 1*fact
   points)
 * nosig, below which sigma we deduct 1*fact point
        * fact, number of points
 * sigdays, how many days we average over to determin mean, stddev
   and thus sigma

The relevance of 'fact' is not just scaling of the unlikely score: there is
a zone in which no points are awarded

### Picking the parameters
Although the parameters can be picked 'by eye', the 8 dimensional space can
be scanned using the tool 'optim'. Optim can be fed with Hare and Hound
synthetic data (so called 'True' files)'. Optim then evolves all 8
parameters to minimize:

 * Falsely detecting a mode
 * Not detecting a mode that is there

And optimize:

 * The number of modes covered

'optim' prints the values of the 8 parameters found, and also creates lots
of files called 'fits.n' which describe which modes were matched using the
model values.

These can then be plotted using (for example) gnuplot:

	plot 'fit.153' using 1:(100*$2) with impulses, 'allBert/Pollux_True.txt' using ($1/1000):(-0.1*$4) with impulses

Optim needs to be passed the name of the True file as a parameter:

	$ ./optim allBert/Pollux_True.txt

