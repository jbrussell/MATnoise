# MATnoise
### MATLAB package to perform ambient noise tomography.

This package consists of two parts: (1) Calculation of ambient noise cross-spectra and measuring interstation phase velocities and (2) inverting interstation velocities for 1D or 2D isotropic phase velocity maps and 1D or 2D azimuthal anisotropy.

> :warning: **WARNING**: This code has only been tested with data that are in 24 hour segments and have been downsampled to 1 Hz!

### 1. Ambient noise & phase velocities (a1-a7)
- **a1** - Calculate ambient noise cross correlations (vertical, radial, transverse) in the frequency domain.

  **Processing options:**
  - One-bit normalization and spectral whitening after Bensen et al. (2007) GJI; [DOI:10.1111/j.1365-246X.2007.03374.x](https://academic.oup.com/gji/article/169/3/1239/626431). 
  - Frequency-time normalization (FTN) after Ekström et al. (2009) GRL; [DOI:10.1029/2009GL039131](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2009GL039131) and Shen et al. (2012) BSSA; [DOI:10.1785/0120120023](https://pubs.geoscienceworld.org/ssa/bssa/article/102/4/1872-1877/325525) rather than typical one-bit normalization and whitening. 
  - Simple prefiltering of seismograms prior to cross-correlating following Ekström (2011) EPSL [DOI: 10.1016/j.epsl.2013.11.022](https://www.sciencedirect.com/science/article/pii/S0012821X13006560?via%3Dihub)

  Assumed directory structure for data:
  ````{datadirectory}/{station}/{station}.{yyyy}.{jday}.{hh}.{mm}.{SS}.{COMP}.sac````\
  e.g.: ````mydata/CC05/CC05.2018.112.00.00.00.BDH.sac````

- **a2** - Plot empirical Green's functions in the time domain
- **a3** - Apply group velocity window around the surface wave arrivals of interest. This as original developed in order to remove the water column arrival from ocean bottom datasets but should also clean up land spectra.
- **a4** - Plot power spectral density of the cross spectra in order to determine the frequency content of the dominant signal.
- **a5** - Plot the real and imaginary parts of the cross spectra to evaluate where signal is best and degree of bias from inhomogeneous noise sources.
- **a6** - Fit J_0 Bessel functions to real part of cross-spectra to extract interstation phase velocity dispersion.
  
  Menke & Jin (2015) BSSA; [DOI:10.1785/0120140245](https://pubs.geoscienceworld.org/ssa/bssa/article/105/3/1619/332118/Waveform-Fitting-of-Cross-Spectra-to-Determine)
- **a7** - Try fitting the phase velocities assuming 1D isotropic velocity variations and 1D azimuthal anisotropy. Azimuthal anisotropy includes options for both 2-theta and 4-theta variations.

### 2. ray_tomo: Ray tomography for 2D structure (b0-b1)
- **b0** - Plot ray paths between station pairs colored by velocity.
- **b1** - Use simple travel-time ray tomography to invert for 2D maps of isotropic phase velocity and 1D or 2D azimuthal anisotropy.
