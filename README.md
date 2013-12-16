astro
=====
Hi! This code depends on NFFTLS for Lomb-Scargle periodograms. NFFTLS can be
found on http://corot-be.obspm.fr/code.php

NFFTLS in turn depends on libfftw3-dev, which is available in most repositories.

NFFTLS also needs the NFFT library (http://www-user.tu-chemnitz.de/~potts/nfft/download.php)
version 3.2.0 or later. 

NFFTLS expects to find NFFT in <nfft/somewhere>, but that's nto where NFFT
installs.  Edit nfft.c to remove 'nfft/' in the #include, and things will
work.
