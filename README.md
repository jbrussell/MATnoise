# MATnoise
### MATLAB package to perform ambient noise tomography.

This package consists of two parts: (1) Calculation of ambient noise cross-spectra and measuring interstation phase velocities and (2) inverting interstation velocities for 1D or 2D isotropic phase velocity maps and 1D or 2D azimuthal anisotropy.

### 1. Ambient noise & phase velocities (a1-a7)
- **a1** - Calculate ambient noise cross correlations (vertical, radial, transverse) in the frequency domain following Bensen et al. (2007) GJI; [DOI:10.1111/j.1365-246X.2007.03374.x](https://academic.oup.com/gji/article/169/3/1239/626431). This code assumes the data are in 24 hour segments and have been downsampled to 1 Hz.

  Assumed directory structure for data: <br/><br/>
  ````{datadirectory}/{station}/{station}.{yyyy}.{jday}.{hh}.{mm}.{SS}.{COMP}.sac````

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
