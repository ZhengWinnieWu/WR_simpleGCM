# WR_simpleGCM
#### Improving GFDL's Idealized General Circulation Model

Our experiments are similar to the Held and Suarez (1994) benchmark. We start our experiments using an analytically determined equilibrium temperature (Teq) and Newtonian relaxation time scale (Tau) profile after Jucker et al. (2014). We then gradually optimize Teq by introducing zonal asymmetries into it using the iterative procedure following Chang (2006). The Teq we provide here corresponds to iteration 31 of DRAG experiment D3 with the surface drag being 0.9 1/day. The file contains 12 monthly Teq fields, each varying in the three spatial dimensions. We also employ an actual orography with land and ocean. Given these modifications, the GCM produces temperatures and diabatic heating that are very similar to that of the reanalysis. A detailed description of the improved model and the outcomes of the various drag experiments are described in Wu and Reichler (2018a, 2018b). 

This code is based on **JFV-strat** (https://github.com/mjucker/JFV-strat). The following program and modules were modified:<br />
•	atmos_model<br />
•	hs_forcing_mod<br />
•	atmosphere_mod<br />
•	time_manager_mod	

The new code reads in twelve monthly Teq and Tau profiles during initialization and then uses linear interpolation to calculate daily fields. The interpolation is done only once during initialization to reduce the actual run time. Additional modification concerns the inclusion of a more flexible surface drag. The Rayleigh damping (surface drag) in the original code is defined in module “hs_forcing_mod”, which is a fixed number. In our modification, we introduce the Rayleigh damping as a two-dimensional matrix (longitude-latitude), and allow the code to read in the Rayleigh damping data files. Certain new variables are read in through the name-list: <br />
•	equilibrium_option = 'from_file', <br />
•	drag_file_name = 'name of drag data file.txt', <br />
•	temp_file = 'name of Teq data file.nc', <br />
•	tau_file = 'name of Tau data file.nc'. <br />

The main references are:<br />
Chang, E. K. M., 2006: An idealized nonlinear model of the Northern Hemisphere winter storm tracks. J. Atmos. Sci., 63, 1818–1839.<br />
Held, I. I. M., and M. M. J. Suarez, 1994: A proposal for the intercomparison of the dynamical cores of atmospheric general circulation models. Bull. Amer. Meteor. Soc., 75, 1825–1830.<br />
Jucker, M., S. Fueglistaler, and G. K. Vallis, 2014: Stratospheric sudden warmings in an idealized GCM. J. Geophys. Res. Atmos., 119, 11 054–11 064, doi:10.1002/2014JD022170.<br />
Wu, Z. and T. Reichler (2018a): Towards a More Earth-like Circulation in an Idealized Model, Geophys. Res. Lett. (submitted).<br />
Wu, Z. and T. Reichler (2018b): Surface Control of the Frequency of Stratospheric Sudden Warming Events, J. Climate (in preparation).
