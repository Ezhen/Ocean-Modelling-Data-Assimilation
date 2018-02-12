# Ocean-Modelling-Data-Assimilation
Here are the scripts built in frames of the following courses of Alex Barth (http://modb.oce.ulg.ac.be/mediawiki/index.php/User:Alex): Ocean Data Assimilation & Structure and Application of Numerical Ocean Models.


* 3D_var.py / Aplication of 3D-Var approach to assimilate on an abstract case
* Covariance.py / Covariance approach to restore the correct temperature based on speeds of currents
* Ensemble_Karman_filter.py / Application of the ensemble Kalman filter to correct a path of a random buoy
* Karman_filter.py / Application of the Kalman filter to correct a path of a random buoy
* SSH_to_geostrophic_currents.py / Recalculation of sea surface elevations to currents
* Wave_energy_integration_shallow_water.py / The examination exercise for the Structure and Application of Numerical Ocean Models. Studying, how a change in domain discretization will affect the integrated wave energy, behaving accordig to the 1D system of shallow water equation.
* Shallow_water_Karman_filter.py / The examination exercise for the Ocean Data Assimilation. Application of the Kalman filter on the 1D shallow water equation system for restoration of the correct shape of a wave trapped inside of a domain in order to explore its efficiency as a function of the time step, availability of data and a random error in the observations.

Addinitonal netcdf files used by some of the scripts above:
* ocean_his_surf.nc / 
* ssh_20071127.nc / 
