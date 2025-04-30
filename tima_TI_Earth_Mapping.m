function [] = tima_TI_Earth_Mapping(TData,MData,inDIR,outDIR,row,varargin)
%% TIMA_TI_EARTH_MAPPING (Thermal Inertia for Mars Analogs)
%   Surface energy balance model for deriving thermal inertia in terrestrial sediments using diurnal
%   observations taken in the field to fit 1D multi-parameter model to each pixel. Justification for
%   approach: https://www.mathworks.com/help/gads/table-for-choosing-a-solver.html
%
% Description
%   This model uses a surrogate optimization solver with 1-3 paramters to fit IR-image-observed
%   surface temperatures using meteorological data and surface parameters derived from
%   aerial/satellite imagery. Fitting Parameters can include any or all of: top layer bulk dry thermal
%   conductivity, Depth of thermal conductivity transition, bottom layer bulk dry thermal conductivity. The tool is
%   optimized for raster files to run one row at a time (or in parallel per compute node)
%   and then run columns in parallel on the same node. (With respect to
%   tower script, this script requires new values for MData.col_min,
%   MData.col_max, MData.lbound, MData.UAV_flight_times, MData.ubound, MData.vars_init,MData.erf)
%
% Syntax
%   tima_TI_Earth_Mapping(TData,MData,inDIR,outDIR,row)
%   tima_TI_Earth_Mapping(TData,MData,inDIR,outDIR,row,varargin)
%   tima_TI_Earth_Mapping(TData,MData,inDIR,outDIR,row,'Mode','2layer')
%
% varargin options:
%   'IncreaseSampling': Option to increase sampling rate if model does not
%       converge. Increases compute time, but avoids nonconvergence. (default=false) (logical)
%   'Initialize': Option to reInitialize subsurface temperatures for each run. Increases compute time, 
%       but avoids dramatic shifts with parameter optimization leading to instability. (default=false) (logical)
%   'Mode': Fitting mode (options: '1layer', '2layer', '2layer_fixed_lower','2layer_fixed_depth' 
%       -- 2Layer includes fitting of bulk dry thermal conductivity of lower 
%       layer and the transition depth) (default='1layer') (string)
%   'Parallel': Run surrogate optimizer in parallel. (default=false)
%       (logical) NOTE: The script is already double parallelized for
%       raster analysis (row and column). This parameter is an additional
%       option for a third level of parallelization, which could cause
%       issues.
%
% Input Parameters
%   TData: Time data - struct of timeseries data variables (all vectors)
%       TData.air_Temp_C: [C] near surface air temperature, typically at 3 m AGL
%       TData.DF: [decimal fraction] Fraction of  Global Horizontal Irradiance (GHI)
%           or r_short_upper that is diffuse (see
%           https://github.com/sandialabs/MATLAB_PV_LIB)
%       TData.evap_depth: [m] Depth of the evaporation front.
%       TData.TIMESTAMP: datetime array associated with each table row
%           spaced by Mdata.dt [s] (datetime)
%       TData.humidity: [decimal fraction] array of near surface relative humidity values,
%           typically at 3m AGL TData.r_long_upper: [W/m^2] Integrated longwave radiation (4.5 to 42 μm) incident on flat
%           surface
%       TData.pressure_air_Pa: [Pa] station pressure, typically at 3m AGL
%       TData.r_long_upper: [W/m^2] Integrated longwave radiation (4500 to 42000 nm) incident on flat
%           surface
%       TData.r_short_lower: [W/m^2] Integrated upwelling shortwave radiation (305 to 2800 nm) from flat
%           surface
%       TData.r_short_upper: [W/m^2] Integrated shortwave radiation (305 to 2800 nm) incident on flat
%           surface
%       TData.solarazimuth_cwfromS: [degrees] Solar azimuth in degrees
%           clockwise from South, typically -180:180 or 0:360
%       TData.solarzenith_apparent: [degrees] Solar zenith in degrees,
%           corrected for atmospheric refraction.
%       TData.VWC_column: [decimal fraction by volume] array of volumetric water content
%           for each model layer with each time step, typically
%           interpolated.
%       TData.windspeed_horiz_ms: [m/s] Near surface horizontal wind speed,
%           typically at 3m AGL.
%   
%   MData: Model Data - Struct of static and model format variables
%       MData.col_min: For reduced column range, minimum (vector)
%       MData.col_max: For reduced column range, maximum (vector)  
%       MData.density: [kg/m^3] Value for density of soil beneath tower.  (vector)
%       MData.dt: [s] Time step (vector)
%       MData.emissivity: [0-1] Weighted thermal emissivity over wavelength
%           range of sensor. (vector)
%       MData.erf: Uncertainty as function of observed temperature (function_handle)
%       MData.k_dry_std_mantle: [W/mK] bulk dry thermal conductivty of topmost
%           mantling layer at T_std (scalar)
%       MData.layer_size: [m] List of vertical thickness for each layer
%           from top to bottom. (vector)
%       MData.lbound: List of lower limits on variables being fit for, in
%           same order as MData.vars_init (vector, size MData.nvars)
%       MData.mantle_thickness: [m] thickness of topmost mantling layer. (scalar)
%       MData.material: ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at the surface (char)
%       MData.material_lower:  ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at depth (char)
%       MData.minit: Vector of initialized test variables with as many
%           randomized samples as desired for fitting, 50 is good such that
%           vector is nvarx50 (vector)
%       MData.ndays: Number of days to run initialization script for
%           (scalar)
%       MData.notes: Details to record in data structure (string)
%       MData.nstep: Number of iterations for curve fitting, 250 is good
%           (vector)
%       MData.nvars: Number of variables being fit for (vector)
%       MData.T_deep: [K] Lower boundary condition, fixed temperature (vector)
%       MData.T_start: [K] Initial condition, list of center temperatures
%           for each layer at start of simulation (vector)
%       MData.T_std: [K] Standard temperature; typically 300 (vector)
%       MData.UAV_flight_times: list of capture times of thermal mosaics of field region (datetime)
%       MData.ubound: List of upper limits on variables being fit for, in
%           same order as MData.vars_init (vector, size MData.nvars)
%       MData.vars_init: [k-upper [W/mK], Pore network con. par. (mk) [unitless],...
%           Surf. ex. coef. (CH) [unitless], Surf. ex. coef. (CE) [unitless], Soil Moist. Infl. (thetak) [% by volume],...
%           Soil Moist. Infl. (thetaE) [% by volume], (Transition Depth [m]), (k-lower [W/mK])]
%           List of 6-8 inputs for variables to serve as either initial or fixed values. (vector)
%
%   inDIR: Full path to directory of inputs (string), must contain
%       Albedo.csv [values 0-1]
%       Aspect_cwfromS.csv [values in degrees]
%       Slope.csv [values in degrees]
%       Shadows: A subdirectory containing modeled shadows throughout the
%           day in format 'Shadow_yyyyMMdd_HHmmss.csv' [values 0-1, with 0 being full shadow]
%       For each thermal map: TempC_#.csv [values in degrees C, where # is the map identified in
%           chronological order]
%
%   outDIR: Full path to directory for outputs (string)
%
%   row: row number of raster to run fitting on (vector)
%
% Outputs:
%   TK_Line_#.txt: [W/mK] Text file containing derived bulk dry thermal
%       conductivity values for row #.
%   Depth_Line_#.txt: [m] Text file containing derived depth to lower layer for row #.
%   fval_Line_#.txt: [chi_v^2] Text file containing derived goodness of fit
%       values for row #.
%   **Row files can be recombined with tima_mapping_recombine_rows.m**
%
% Author
%    Ari Koeppel -- Copyright 2023
%
% Sources
%   Optimization tool: https://www.mathworks.com/help/gads/table-for-choosing-a-solver.html
%   Subsurface heat flux: Kieffer et al., 2013
%   Irradiance Calculation: https://github.com/sandialabs/MATLAB_PV_LIB
%   Sensible heat flux:
%   Latent heat flux: Daamen and Simmonds 1996 + Mahfouf and Noilhan (1991)+ Kondo1990
%   Thermal conductivity mixing:
%   
% See also 
%   time_heat_transfer.m tima_initialize.m tima_ln_prior.m
%   tima_mapping_recombine_rows.m

format shortG
p = inputParser;
p.addRequired('TData',@isstruct);
p.addRequired('MData',@isstruct);
p.addRequired('inDIR',@ischar);
p.addRequired('outDIR',@ischar);
p.addParameter('Mode','1layer',@ischar);%'2layer','2layer_fixed_depth','2layer_fixed_lower'
p.addParameter('IncreaseSampling',false,@islogical);
p.addParameter('Initialize',false,@islogical);
p.addParameter('Parallel',false,@islogical);
p.parse(TData, MData, inDIR, outDIR, varargin{:});
p=p.Results;
Initialize = p.Initialize;
IncreaseSampling = p.IncreaseSampling;

imptopts = detectImportOptions([inDIR,'Slope.csv']);
imptopts.DataLines = [row row];
Data_Slope_X = readmatrix([inDIR,'Slope.csv'],imptopts);
Data_Aspect_X = readmatrix([inDIR,'Aspect_cwfromS.csv'],imptopts);
Data_Albedo_X = readmatrix([inDIR,'Albedo.csv'],imptopts); 
ShadowDataDir = [inDIR,'Shadows/'];
ShadowFiles = dir(fullfile(ShadowDataDir,'*.csv')); %gets all files with yyyyMMdd_HHmmss.csv suffix
Data_Shadows_X = NaN([size(Data_Albedo_X,2) length(ShadowFiles)]); 
Shadow_Times = datetime.empty(length(ShadowFiles),0);
for k = 1:length(ShadowFiles)
    FileName = fullfile(ShadowDataDir, ShadowFiles(k).name);
    Data_Shadows_X(:,k) = readmatrix(FileName,imptopts);
    Shadow_Times(k) = datetime(ShadowFiles(k).name(end-18:end-4),'InputFormat',"yyyyMMdd_HHmmss");
end
shadow_time_ind = NaN([length(TData.TIMESTAMP) 1]);
for k = 1:length(TData.TIMESTAMP)
    [~, shadow_time_ind(k)] = min(abs(timeofday(Shadow_Times) - timeofday(TData.TIMESTAMP(k))));
end
UAV_flight_ind = NaN([length(MData.UAV_flight_times) 1]);
for k = 1:length(MData.UAV_flight_times)
    [~, UAV_flight_ind(k)] = min(abs(MData.UAV_flight_times(k) - TData.TIMESTAMP));
end
if size(shadow_time_ind,2)>size(shadow_time_ind,1), shadow_time_ind=shadow_time_ind';end
Data_UAV_X = NaN([size(Data_Albedo_X,2) size(UAV_flight_ind,2)]);
for t = 1:size(UAV_flight_ind,2)
    Data_UAV_X(:,t) = readmatrix([inDIR,sprintf('TempC_%u.csv',t)],imptopts);
end
if strcmp(p.Mode,'1layer')
        MData.minit = MData.minit(1,:);
elseif strcmp(p.Mode,'2layer')
        MData.minit = MData.minit([1 7:8],:);
elseif strcmp(p.Mode,'2layer_fixed_depth')
        MData.minit = MData.minit([1 8],:);
elseif strcmp(p.Mode,'2layer_fixed_lower')
        MData.minit = MData.minit([1 7],:);
end
poolobj = gcp('nocreate');
delete(poolobj);
%parpool('Processes',MData.workers)
%% Run Model in 1 point mode for each pixel
%row = defined in function call
for col = MData.col_min:MData.col_max
        RESULTS = NaN([1 MData.nvars]);
        fval = NaN;
        if any(isnan(Data_UAV_X(col,1:end))) || isnan(Data_Slope_X(col)) || isnan(Data_Aspect_X(col)) || isnan(Data_Albedo_X(col)) || isnan(single(Data_Shadows_X(col,1)))
            continue
        end

        if strcmp(p.Mode,'1layer')
            formod = @(theta) tima_full_model(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,TData.temps_to_fit_interp,MData.ndays,'albedo',Data_Albedo_X(col),'slope_angle',Data_Slope_X(col),...
                'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',shadow_time_ind,'material',MData.material,'mantle_thickness',MData.mantle_thickness,...
                'k_dry_std_mantle',MData.k_dry_std_mantle,IncreaseSampling',IncreaseSampling,'Initialize',Initialize);
        elseif strcmp(p.Mode,'2layer')
            formod = @(theta) tima_full_model(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,TData.temps_to_fit_interp,MData.ndays,'albedo',Data_Albedo_X(col),...
                'slope_angle',Data_Slope_X(col),'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',shadow_time_ind,'material',MData.material,'depth_transition',...
                 theta(2),'k_dry_std_lower',theta(3),'material_lower',MData.material_lower,'mantle_thickness',MData.mantle_thickness,...
                'k_dry_std_mantle',MData.k_dry_std_mantle,'IncreaseSampling',IncreaseSampling,'Initialize',Initialize);
        elseif strcmp(p.Mode,'2layer_fixed_depth')
            formod = @(theta) tima_full_model(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,TData.temps_to_fit_interp,MData.ndays,'albedo',Data_Albedo_X(col),'slope_angle',Data_Slope_X(col),...
                'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',shadow_time_ind,'material',MData.material,'depth_transition',...
                 MData.vars_init(7),'k_dry_std_lower',theta(2),'material_lower',MData.material_lower,'mantle_thickness',MData.mantle_thickness,...
                'k_dry_std_mantle',MData.k_dry_std_mantle,'IncreaseSampling',IncreaseSampling,'Initialize',Initialize);
        elseif strcmp(p.Mode,'2layer_fixed_lower')
            formod = @(theta) tima_full_model(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.evap_depth,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,TData.temps_to_fit_interp,MData.ndays,'albedo',Data_Albedo_X(col),...
                'slope_angle',Data_Slope_X(col),'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',shadow_time_ind,'material',MData.material,...
                'depth_transition',theta(2),'k_dry_std_lower',MData.vars_init(8),'material_lower',MData.material_lower,'mantle_thickness',MData.mantle_thickness,...
                'k_dry_std_mantle',MData.k_dry_std_mantle,'IncreaseSampling',IncreaseSampling,'Initialize',Initialize);
        end
        Temps_Obs = Data_UAV_X(col,:);
        Temps_Obs = Temps_Obs(:);
        opts = optimoptions('surrogateopt','InitialPoints',MData.minit,'UseParallel',p.Parallel,'MaxFunctionEvaluations',MData.nstep);
        Obj = @(theta) tima_fval_chi2v(Temps_Obs,tima_formod_subset(theta,UAV_flight_ind,formod),MData.erf(Temps_Obs),MData.nvars);
        problem = struct('solver','surrogateopt','lb',MData.lbound,'ub',MData.ubound,'objective',Obj,'options',opts,'PlotFcn',[]) ; 
        [RESULTS_holder,fval_holder] = surrogateopt(problem);
        if isempty(RESULTS_holder) || isempty(fval_holder)
            continue
        else
            RESULTS(:) = RESULTS_holder';
            fval = fval_holder;
            writematrix(round(RESULTS(1),3),[outDIR,sprintf('TK_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            writematrix(round(fval,3),[outDIR,sprintf('fval_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            if strcmp(p.Mode,'2layer')
                writematrix(round(RESULTS(2),3),[outDIR,sprintf('Depth_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
                writematrix(round(RESULTS(3),3),[outDIR,sprintf('TK_lower_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            elseif strcmp(p.Mode,'2layer_fixed_depth')
                writematrix(round(RESULTS(2),3),[outDIR,sprintf('TK_lower_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            elseif strcmp(p.Mode,'2layer_fixed_lower')
                writematrix(round(RESULTS(2),3),[outDIR,sprintf('Depth_Row_%u_Col_%u.txt',row,col)],'Delimiter',',')
            end
        end
end
poolobj = gcp('nocreate');
delete(poolobj);
%% Combine data into rows
LineTKData = NaN([1 MData.col_max]);
LinefvalData = NaN([1 MData.col_max]);
if strcmp(p.Mode,'2layer')
    LineTKlowerData = NaN([1 MData.col_max]);
    LineDepthData = NaN([1 MData.col_max]);
elseif strcmp(p.Mode,'2layer_fixed_depth')
    LineTKlowerData = NaN([1 MData.col_max]);
elseif strcmp(p.Mode,'2layer_fixed_lower')
    LineDepthData = NaN([1 MData.col_max]);
end
for col = 1:MData.col_max
    LineTKFile = [out_DIR,sprintf('TK_Row_%u_Col_%u.txt',row,col)];
    LineTKlowerFile = [out_DIR,sprintf('TK_lower_Row_%u_Col_%u.txt',row,col)];
    LineDepthFile = [out_DIR,sprintf('Depth_Row_%u_Col_%u.txt',row,col)];
    LinefvalFile = [out_DIR,sprintf('fval_Row_%u_Col_%u.txt',row,col)];
    if isfile(LineTKFile)
        LineTKData(1,col) = readmatrix(LineTKFile);
        LinefvalData(1,col) = readmatrix(LinefvalFile);
    end
    if isfile(LineTKlowerFile)
        LineTKlowerData(1,col) = readmatrix(LineTKlowerFile);
    end
    if isfile(LineDepthData)
	    LineDepthData(1,col) = readmatrix(LineDepthFile);
    end
end
writematrix(LineTKData,[outDIR,sprintf('TK_Line_%u.txt',row)],'Delimiter',',')
writematrix(LinefvalData,[outDIR,sprintf('fval_Line_%u.txt',row)],'Delimiter',',')
if strcmp(p.Mode,'2layer')
    writematrix(LineTKlowerData,[outDIR,sprintf('TK_lower_Line_%u.txt',row)],'Delimiter',',')
    writematrix(LineDepthData,[outDIR,sprintf('Depth_Line_%u.txt',row)],'Delimiter',',')
elseif strcmp(p.Mode,'2layer_fixed_depth')
    writematrix(LineTKlowerData,[outDIR,sprintf('TK_lower_Line_%u.txt',row)],'Delimiter',',')
elseif strcmp(p.Mode,'2layer_fixed_lower')
    writematrix(LineDepthData,[outDIR,sprintf('Depth_Line_%u.txt',row)],'Delimiter',',')
end
delete([outDIR,sprintf('*Row_%u_Col*.txt',row)])
end