# WR_simpleGCM
#### Further improvement of GFDL's idealized General Circulation Model

We employ actual orography with land and ocean in the idealized model. The Teq we provide corresponds to iteration 31 of DRAG experimeent D3 as described in Wu and Reichler (2018a). It is the result of an iterative procedure after Chang (2006) and leads to reanalysis-like temperatures and diabatic heating. We also vary the magnitude of the Rayleigh damping coefficients. A description of the improved model and corresponding outcomes are in 

Wu, Z. and T. Reichler (2018a): Towards a More Earth-like Circulation in an Idealized Model, Geophys. Res. Lett. (submitted).

We further use this modified model to study the surface control of the frequency of stratospheric sudden warming events, as described in 

Wu, Z. and T. Reichler (2018b): Surface Control of the Frequency of Stratospheric Sudden Warming Events, J. Climate (in preparation).

We base on Held and Suarez (1994) benchmark, and start our experiments using an analytically determined Teq and Tau profile given in the appendix of Jucker et al. (2014). The main references are

Held, I. I. M., and M. M. J. Suarez, 1994: A proposal for the intercomparison of the dynamical cores of atmospheric general circulation models. Bull. Amer. Meteor. Soc., 75, 1825–1830.

Jucker, M., S. Fueglistaler, and G. K. Vallis, 2014: Stratospheric sudden warmings in an idealized GCM. J. Geophys. Res. Atmos., 119, 11 054–11 064, doi:10.1002/2014JD022170.

Chang, E. K. M., 2006: An idealized nonlinear model of the Northern Hemisphere winter storm tracks. J. Atmos. Sci., 63, 1818–1839.

This code is based on JFV-strat (https://github.com/mjucker/JFV-strat). The following modules were modified: ... The new code reads in twelve monthly Teq profiles during initialization and then calculates daily Teq fields by linear interpolation. The entire interpolation is done only once during initialization to reduce the run time. Additional modifications concern ... Since several modules of the original JFV-strat code were modified, it is probably best to replace the entire "src" folder. 
 
