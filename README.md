# WR_simpleGCM
Further improvement of GFDL's idealized General Circulation Model

We employ actual orography with land and ocean into the idealized model. We use an iterative prodedure after Chang (2006) to determine equilibrium temperatures (Teq) that lead to reanalysis-like temperatures and diabatic heating.  We also vary the magnitude of the Rayleigh damping coefficients. The detailed discription of the improved model and the corresponding outcomes are in 

Wu, Z. and T. Reichler (2018): Towards a More Earth-like Circulation in an Idealized Model, Geophys. Res. Lett. (submitted).

We further use this modified model to study the surface controls of the frequency of stratospheric sudden warming events, as described in 

Wu, Z. and T. Reichler (2018): Surface Control of the Frequency of Stratospheric Sudden Warming Events, J. Climate (in preparation).

We base on Held and Suarez (1994) benchmark, and start our experiments using an analytically determined Teq and Tau profile given in the appendix of Jucker et al. (2014). The main references are

Held, I. I. M., and M. M. J. Suarez, 1994: A proposal for the intercomparison of the dynamical cores of atmospheric general circulation models. Bull. Amer. Meteor. Soc., 75, 1825–1830.

Jucker, M., S. Fueglistaler, and G. K. Vallis, 2014: Stratospheric sudden warmings in an idealized GCM. J. Geophys. Res. Atmos., 119, 11 054–11 064, doi:10.1002/2014JD022170.

Chang, E. K. M., 2006: An idealized nonlinear model of the Northern Hemisphere winter storm tracks. J. Atmos. Sci., 63, 1818–1839.

This code uses "read from file" to get the Teq profile. In order to make the code more efficient, this modified code reduces the times for interpolation. It is better to replace the whole "src" folder to replace your "src" folder if you already have the experiment from the Flexible Modeling System (FMS). 
