![logo](https://github.com/user-attachments/assets/e6d43934-838c-4137-9df0-3e163166e27c)

# Thermal Inertia for Mars Analogs (TIMA) Model

 TIMA is a finite difference numerical model of soil heat transfer that that considers time-variability in thermal conductivity as a function of temperature, moisture, and composition. The detailed parametrization of thermal conductivity is meant to enable direct translation between soil thermophysical properties on Earth and Mars. This repository is maintained by [Ari Koeppel](https://earthsciences.dartmouth.edu/people/ari-koeppel). Please cite Koeppel et al., (2024) when using the model.

 Note: the 2D implementation of the model using raster data is requiring some reconstruction and will be available here soon!

# References:
Koeppel, A.H., Edwards, C.S., Edgar, L.A., Nowicki, S., Bennett, K.A., Gullikson, A., Piqueux, S., Eifert, H., Chapline, D. and Rogers, A.D., 2024. A novel surface energy balance method for thermal inertia studies of terrestrial analogs. Earth and Space Science, 11(9), p.e2023EA003259.

# Overview: 
The original implementation of the model was designed to derive a characteristic dry soil thermal conductivity at 300K. The procedure achieves this by inputing micrometeorological data into a surface energy balance forward model and adjusting thermal conductivity (along with 5 other modifying parameters) to fit surface skin temperature data (typically obtained through radiometer observation). The preferred fitting protocol involves Markov Chain Monte Carlo simulations to identify the most probable set of constants and associated uncertainties. Model updates have included multilayer parametrizations, including consideration of subsurface ice. Thus, the model can be used to derive soil physical properties, layering and the depth to subsurface material transitions based entirely on surface temperature and weather observations.

The required fitting parameters are:
1.      k_dry_std_upper: [W/mK] bulk dry thermal conductivty of upper layer at T_std (scalar)
2.      m: [unitless] pore-network-connectivity-parameter, typically 0-1.3 (scalar)
3.      CH: [Unitless] resistance to sensible heat flux coefficient,
          similar to the aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2)
          1/CH should be 0.0028-0.0075 (or CH~100-400) for smooth to
          roughly open soils on Davenport Scale, CH is larger (more
          resistence) with rougher topography or larger vegetation (scalar)
4.      CH: [Unitless] resistance to sensible heat flux coefficient,
          similar to the aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2)
          1/CH should be 0.0028-0.0075 (or CH~100-400) for smooth to
          roughly open soils on Davenport Scale, CH is larger (more
          resistence) with rougher topography or larger vegetation (scalar)
5.       theta_k: [0-1, fraction by volume] conductivity soil moisture inflection
          point - in theory this should be similar to saturation value or porosity (scalar) 
6.       theta_E: [0-1, fraction by volume] latent heat soil moisture inflection
          point (scalar)
     
and optional parameters:

7.       'depth_transition': [m] Depth of major transition between upper
          material and lower material. (scalar)
8.       'k_dry_std_lower': [W/mK] bulk dry thermal conductivty of lower layer at T_std (scalar)


# Forward Model:
    [T_surf_C,T_sub_C,q_latent,k_eff_dt,q_conv,q_rad,q_G] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,  windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,pressure_air_pa,varargin)

## Input Parameters ('single quotes' means optional varargin)
  ### Timeseries data
      air_Temp_C: [C] near surface air temperature, typically at 3 m AGL.
          (1D vector)
      'albedo': [0-1] time variant albedo (e.g., due to wetting) for
          fitting surface. (1D vector)
      evap_depth: [m] Depth of the evaporation front (e.g., as function
          of VWC). (1D vector)
      'f_diff': [0-1] Fraction of  Global Horizontal Irradiance (GHI)
          or r_short_upper that is diffuse. (1D vector)
      RH: [0-1] array of near surface relative humidity values,
          typically at 3m AGL. (1D vector)
      pressure_air_pa: [Pa] station pressure, typically at 3m AGL. (1D vector)
      r_long_upper: [W/m^2] Integrated longwave radiation  (4.5 to 42 μm) incident on flat
          surface.(1D vector)
      r_short_lower: [W/m^2] Integrated upwelling shortwave radiation (305 to 2800 nm) from flat
          surface. (1D vector)
      r_short_upper: [W/m^2] Integrated shortwave radiation (305 to 2800 nm) incident on flat
          surface. (1D vector)
      'shadow_time_ind': Index of shadow_data corresponding to
          timeseries. (1D vector)
      'solar_azimuth_cwfromS': [degrees] Solar azimuth in degrees
          clockwise from South, typically -180:180 or 0:360 (1D vector)
      'solar_zenith_apparent': [degrees] Solar zenith in degrees,
          corrected for atmospheric refraction. (1D vector)
      VWC_column: [0-1, decimal fraction by volume] array of volumetric water content
          for each model layer with each time step, typically
          interpolated. (2D vector)
      windspeed_horiz: [m/s] Near surface horizontal wind speed,
          typically at 3m AGL. (1D vector)

  ### Constants:
      'aspect_cwfromS': [degrees] slope aspect clockwise from South. (scalar)
      CE: [Unitless] resistance to latent heat flux coefficient, similar
          to the aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2) (scalar)
      CH: [Unitless] resistance to sensible heat flux coefficient,
          similar to the aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2)
          1/CH should be 0.0028-0.0075 (or CH~100-400) for smooth to
          roughly open soils on Davenport Scale, CH is larger (more
          resistence) with rougher topography or larger vegetation (scalar)
      'depth_transition': [m] Depth of major transition between upper
          material and lower material. (scalar)
      dt: [s] Time step (scalar)
      emissivity: [0-1] Weighted thermal emissivity over wavelength
          range of sensor. (scalar)
      'e_fxn': spectral emissivity function that determines emissivity as
          function of wavelength in um. (function)
      initial_temps: [K] Initial condition, list of center temperatures
          for each layer at start of simulation (1D vector)
      'k_dry_std_lower': [W/mK] bulk dry thermal conductivty of lower layer at T_std (scalar)
      'k_dry_std_mantle': [W/mK] bulk dry thermal conductivty of topmost
          mantling layer at T_std (scalar)
      k_dry_std_upper: [W/mK] bulk dry thermal conductivty of upper layer at T_std (scalar)
      layer_size: [m] List of vertical thickness for each layer
          from top to bottom. (1D vector)
      m: [unitless] pore-network-connectivity-parameter, typically 0-1.3 (scalar)
      'mantle_thickness': [m] thickness of topmost mantling layer. (scalar)
      'material': ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at the surface (char)
      'material_lower':  ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at depth (char)
      NDAYS: [unitless] # of days to run equilib model for: More = closer to equilib, fewer = faster (scalar)
      'rho_dry_lower': [kg/m^3] Value for density of lower layer soil.  (scalar)
      rho_dry_upper: [kg/m^3] Value for density of soil.  (scalar)
      'shadow_data': [0-1] Fraction of ROI/pixel in shadow corresponding to time. (1D Vector)
      'slope_angle': [degrees] angle of slope. (scalar)
      theta_k: [0-1, fraction by volume] conductivity soil moisture inflection
          point - in theory this should be similar to saturation value or porosity (scalar) 
      theta_E: [0-1, fraction by volume] latent heat soil moisture inflection
          point (scalar)
      T_deep: [K] Lower boundary condition, fixed temperature (scalar)
      T_std: [K] Standard temperature; typically 300 (scalar)
      'T_adj1': [index, temperature K] pair used to force column
          temperature change at a given time point due to wetting
          (1D vector)
      'T_adj2': [index, temperature K] pair used to force a second
          column temperature change at a given time point due to wetting
          (1D vector)

## Output Parameters:
    T_Surf_C = [C] Surface temperature time series (1D vector)
    T_sub_C = [C] Suburface temperature time series for all layers (2D vector)
    q_latent = [W/m^2] Latent heat flux time series for all layers (2D vector)
    k_eff_dt = [W/mK] Effective thermal conductivity time series for all layers (2D vector)
    q_conv = [W/m^2] Surface sensible heat flux time series (1D vector)
    q_rad = [W/m^2] Surface radiative heat flux time series (1D vector)
    q_G = [W/m^2] Heat flux between top and second subsurface layer time series (1D vector)

# Fitting Process Input/Output Setup:
  ### Time data - struct of timeseries data variables
      TData.air_Temp_C: [C] near surface air temperature, typically at 3
          m AGL (1D vector)
      TData.DF: [decimal fraction] Fraction of  Global Horizontal Irradiance (GHI)
          or r_short_upper that is diffuse (1D vector)
      TData.evap_depth: [m] Depth of the evaporation front. (1D vector)
      TData.temps_to_fit_interp: [C] Surface temperature values to be used for
          fitting with any data gaps filled in. (1D vector)
      TData.timed_albedo: [decimal fraction] time variant albedo (e.g.,
          due to wetting) for fitting surface. (1D vector)
      TData.TIMESTAMP: datetime array associated with each table row
          spaced by Mdata.dt [s] (1D datetime)
      TData.humidity: [decimal fraction] array of near surface relative humidity values,
          typically at 3m AGL TData.r_long_upper: [W/m^2] Integrated longwave radiation (4.5 to 42 μm) incident on flat
          surface (1D vector)
      TData.pressure_air_Pa: [Pa] station pressure, typically at 3m AGL (1D vector)
      TData.r_long_upper: [W/m^2] Integrated longwave radiation (4500 to 42000 nm) incident on flat
          surface (1D vector)
      TData.r_short_lower: [W/m^2] Integrated upwelling shortwave radiation (305 to 2800 nm) from flat
          surface (1D vector)
      TData.r_short_upper: [W/m^2] Integrated shortwave radiation (305 to 2800 nm) incident on flat
          surface (1D vector)
      TData.solarazimuth_cwfromS: [degrees] Solar azimuth in degrees
          clockwise from South, typically -180:180 or 0:360 (1D vector)
      TData.solarzenith_apparent: [degrees] Solar zenith in degrees,
          corrected for atmospheric refraction. (1D vector)
      TData.VWC_column: [decimal fraction by volume] array of volumetric water content
          for each model layer with each time step, typically
          interpolated. (2D vector)
      TData.windspeed_horiz_ms: [m/s] Near surface horizontal wind speed,
          typically at 3m AGL. (1D vector)

   ### Model Data - Struct of static and model format variables
      MData.burnin_fit: initial time length (s) to ignore in fitting (vector, default=0)
      MData.burnin_mcmc: fraction of the chain that should be removed.
          (scalar, default=0)
      MData.density: [kg/m^3] Value for density of soil beneath tower.
          (scalar)
      MData.dt: [s] Time step (scalar)
      MData.emissivity: [0-1] Weighted thermal emissivity over wavelength
          range of sensor. (scalar)
      MData.erf: Uncertainty as function of observed temperature (function_handle)
      MData.fit_ind: Indecies of temps_to_fit_interp in which to apply fitting
          to (scalar)
      MData.k_dry_std_mantle: [W/mK] bulk dry thermal conductivty of topmost
          mantling layer at T_std (scalar)
      MData.layer_size: [m] List of vertical thickness for each layer
          from top to bottom. (1D vector)
      MData.mantle_thickness: [m] thickness of topmost mantling layer. (scalar)
      MData.material: ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at the surface (char)
      MData.material_lower:  ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at depth (char)
      MData.minit: Vector of initialized test variables with as many
          randomized samples as desired for fitting, 50 is good such that
          vector is nvarx50 (2D vector)
      MData.notes: Details to record in data structure (string)
      MData.nstep: Number of iterations for curve fitting, 250 is good
          (scalar)
      MData.nvars: Number of variables being fit for (scalar)
      MData.nwalkers: Number of walkers in MCMC ensemble (scalar)
      MData.ThinChain: MCMC data reduction by thinning all output chains by only storing every N'th step (scalar, default=10)
      MData.T_adj1: [index, temperature K] pair used to force column
          temperature change at a given time point due to wetting (1D vector,
          optional)
      MData.T_adj2: [index, temperature K] pair used to force a second
          column temperature change at a given time point due to wetting (1D vector,
          optional)
      MData.T_deep: [K] Lower boundary condition, fixed temperature (scalar)
      MData.T_start: [K] Initial condition, list of center temperatures
          for each layer at start of simulation (scalar)
      MData.T_std: [K] Standard temperature; typically 300 (scalar)
      MData.vars_init: [k-upper [W/mK], Pore network con. par. (mk) [unitless],...
          Surf. ex. coef. (CH) [unitless], Surf. ex. coef. (CE) [unitless], Soil Moist. Infl. (thetak) [% by volume],...
          Soil Moist. Infl. (thetaE) [% by volume], (Transition Depth [m]), (k-lower [W/mK])]
          List of 6-8 inputs for variables to serve as either initial or fixed values. (1D vector)

Output directory: 

      out_DIR: Full path to directory for outputs (string)
    
## Input data:
Example input files are located in the [`Example Data`](hhttps://github.com/NAU-PIXEL/Thermal-Inertia-for-Mars-Analogs/tree/main/Example%20Data) folder.

# Run Fitting Process
    [models,names] = tima_TI_Earth_Tower(TData,MData,outDIR,'Mode','1layer','Parallel',true,'SaveModels',true);

# Visualize Results
    tima_plot_results(TData,MData,models,names,varargin)
    
![TempFit](https://github.com/user-attachments/assets/8ec806ed-61ec-43e2-b552-0ca5e52da66d)
![Corner](https://github.com/user-attachments/assets/4aa489ee-407d-4f55-9cf2-73198073981d)
![SubTemps](https://github.com/user-attachments/assets/550d6d10-7ccc-4662-83e2-4105e22bee93)
![keff](https://github.com/user-attachments/assets/3205714b-81e2-4aff-8482-569936f775ad)
![GroundHeat](https://github.com/user-attachments/assets/0aa85abe-f935-4d81-bfc0-a8eac5020e55)
![Rad](https://github.com/user-attachments/assets/ac3addf0-a153-46ed-bff8-1c209747a36b)
![SensibleHeat](https://github.com/user-attachments/assets/8da21e97-13c3-45cd-9b39-40866374c9c9)
![LatentHeat](https://github.com/user-attachments/assets/e7b38248-414b-465b-9bb4-a355ad493011)
