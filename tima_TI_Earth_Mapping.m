function [] = tima_TI_Earth_Mapping(TData,MData,inDIR,outDIR,row,varargin)
%% TIMA_TI_EARTH_MAPPING_GENERAL (TI Mars Analogs)
%   Surface energy balance model for deriving thermal inertia in terrestrial sediments using diurnal
%   observations taken in the field to fit 1D multi-parameter model to each pixel. Justification for
%   approach: https://www.mathworks.com/help/gads/table-for-choosing-a-solver.html
%
% Description
%   This model uses a surrogate optimization solver with 1-7 paramters to fit IR-image-observed
%   surface temperatures using meteorological data and surface parameters derived from
%   aerial/satellite imagery. Fitting Parameters can include any or all of: top layer thermal
%   conductivity at 300K, bottome layer thermal conductivity at 300K, Depth of thermal conductivity
%   transition, pore-network-connectivity, Surf. ex. coef. (sensible), Surf. ex. coef. (latent), Soil Moist. Inflection (conductivity) (%), Soil Moist. Infl. (latent heat) (%)
%
% Input Parameters
%   TData: Time data - struct of timeseries data variables (all vectors)
%       TData.air_Temp_C: [C] near surface air temperature, typically at 3 m height
%       TData.r_short_upper: [W/m^2] Integrated shortwave radiation (305 to 2800 nm) incident on flat
%           surface
%       TData.r_short_lower: [W/m^2] Integrated upwelling shortwave radiation (305 to 2800 nm) from flat
%           surface
%       TData.r_long_upper: [W/m^2] Integrated longwave radiation (4.5 to 42 μm) incident on flat
%           surface
%       TData.windspeed_horiz_ms: [m/s] Near surface horizontal wind speed,
%           typically at 3m height.
%       TData.VWC_column: [% by volume] array of volumetric water content
%           for each model layer with each time step, typically
%           interpolated.

%       TData.DF:
%       TData.dug_VWC_smooth: 
%       TData.evap_depth 
%       TData.humidity: as a fraction
%       TData.pressure_air_Pa: 
%       TData.r_long_upper:
%       TData.r_short_upper: 
%       TData.r_short_lower: 
%       TData.solarazimuth_cwfromS:
%       TData.solarzenith_apparent:
%       TData.timed_albedo: 
%       TData.TIMESTAMP: 
%       TData.temps_to_fit: 
%       TData.windspeed_horiz_ms: 
%
%   
%   MData: Model Data - Struct of static and model format variables
%       MData.vars_init: [k-upper [W/mK], Pore network con. par. (mk) [unitless],...
%           Surf. ex. coef. (CH) [unitless], Surf. ex. coef. (CE) [unitless], Soil Moist. Infl. (thetak) [% by volume],...
%           Soil Moist. Infl. (thetaE) [% by volume], (Transition Depth [m]), (k-lower [W/mK])]
%           List of 6-8 inputs for variables to serve as either initial or fixed values. (vector)
%       MData.density: [kg/m^3] Value for density of soil beneath tower.  (vector)
%       MData.dt: [s] Time step (vector)
%       MData.T_std: [K] Standard temperature; typically 300 (vector)
%       MData.T_deep: [K] Lower boundary condition, fixed temperature (vector)
%       MData.T_start: [K] Initial condition, list of center temperatures
%           for each layer at start of simulation (vector)
%       MData.layer_size: [m] List of vertical thickness for each layer
%           from top to bottom. (vector)
%       MData.emissivity: [0-1] Weighted thermal emissivity over wavelength
%           range of sensor. (vector)
%       MData.material: ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at the surface (char)
%       MData.material_lower:  ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at depth (char)
%       MData.UAV_flight_times: list of capture times of thermal mosaics of field region (datetime)
%       MData.col_min: For reduced column range, minimum (vector)
%       MData.col_max: For reduced column range, maximum (vector)      
%       Mdata.parallel: true or false whether to run fitting tool in parallel  (logical, default false)
%       MData.nvars: Number of variables being fit for (vector)
%       MData.lbound: List of lower limits on variables being fit for, in
%           same order as MData.vars_init (vector, size MData.nvars)
%       MData.ubound: List of upper limits on variables being fit for, in
%           same order as MData.vars_init (vector, size MData.nvars)
%       MData.notes: Details to record in data structure (string)
%       MData.minit: Vector of initialized test variables with as many
%           randomized samples as desired for fitting, 50 is good such that
%           vector is nvarx50 (vector)
%       MData.nstep: Number of iterations for curve fitting, 250 is good
%           (vector)
%       MData.erf: Uncertainty as function of observed temperature (function_handle)

%
%   inDIR: Full path to directory of inputs (string), must contain
%       Slope.csv [values in degrees]
%       Aspect_cwfromS.csv [values in degrees]
%       Albedo.csv [values 0-1]
%       Fore each thermal map: TempC_#.csv [values in degrees C, where # is the map identified in
%           chronological order]
%       Shadows: A subdirectory containing modeled shadows throughout the
%           day in format 'Shadow_yyyyMMdd_HHmmss.csv' [values 0-1, with 0 being full shadow]
%   out_DIR: Full path to directory for outputs (string)

% Outputs:
%   TK_Line_%u.txt
%   Depth_Line_%u.txt
%   fval_Line_%u.txt
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
%   TIMA_HEAT_TRANSFER TIMA_INITIALIZE TIMA_LATENT_HEAT_MODEL TIMA_LN_PRIOR TIMA_SENSIBLE_HEAT_MODEL TIMA_GWMCMC TIMA_COMBINE_ROWS

format shortG
p = inputParser;
p.addRequired('TData',@isstruct);
p.addRequired('MData',@isstruct);
p.addRequired('inDIR',@ischar);
p.addRequired('outDIR',@ischar);
p.addParameter('Mode','1layer',@ischar);%'2layer','2layer_fixed_depth','2layer_fixed_lower'
p.parse(TData, MData, inDIR, outDIR, varargin{:});
p=p.Results;
if ~isfield(MData, 'parallel')
    Mdata.parallel = false;
end
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

poolobj = gcp('nocreate');
delete(poolobj);
parpool('Processes',10)
%% Run Model in 1 point mode for each pixel
%row = defined in function call
parfor col = MData.col_min:MData.col_max
        RESULTS = NaN([1 MData.nvars]);
        fval = NaN;
        if any(isnan(Data_UAV_X(col,1:end))) || isnan(Data_Slope_X(col)) || isnan(Data_Aspect_X(col)) || isnan(Data_Albedo_X(col)) || isnan(single(Data_Shadows_X(col,1)))
            continue
        end

        if strcmp(p.Mode,'1layer')
            formod = @(theta) tima_heat_transfer(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,'albedo',Data_Albedo_X(col),'slope_angle',Data_Slope_X(col),...
                'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',shadow_time_ind,'MappingMode',true,'material',MData.material);
        elseif strcmp(p.Mode,'2layer')
            formod = @(theta) tima_heat_transfer(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,'albedo',Data_Albedo_X(col),...
                'slope_angle',Data_Slope_X(col),'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',shadow_time_ind,'MappingMode',true,'material',MData.material,'depth_transition',...
                 theta(2),'k_dry_std_lower',theta(3),'material_lower',MData.material_lower);
        elseif strcmp(p.Mode,'2layer_fixed_depth')
            formod = @(theta) tima_heat_transfer(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,'albedo',Data_Albedo_X(col),'slope_angle',Data_Slope_X(col),...
                'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',shadow_time_ind,'MappingMode',true,'material',MData.material,'depth_transition',...
                 MData.vars_init(7),'k_dry_std_lower',theta(2),'material_lower',MData.material_lower);
        elseif strcmp(p.Mode,'2layer_fixed_lower')
            formod = @(theta) tima_heat_transfer(theta(1),MData.vars_init(2),MData.vars_init(3),...
                MData.vars_init(4),MData.vars_init(5),MData.vars_init(6),MData.density,MData.dt,MData.T_std,TData.air_Temp_C,TData.r_short_upper,...
                TData.r_short_lower,TData.r_long_upper,TData.windspeed_horiz_ms,MData.T_deep,MData.T_start,MData.layer_size,...
                TData.VWC_column,TData.humidity,MData.emissivity,...
                TData.pressure_air_Pa,'albedo',Data_Albedo_X(col),...
                'slope_angle',Data_Slope_X(col),'aspect_cwfromS',Data_Aspect_X(col),'solar_azimuth_cwfromS',...
                TData.solarazimuth_cwfromS,'solar_zenith_apparent',TData.solarzenith_apparent,...
                'f_diff',TData.DF,'shadow_data',single(Data_Shadows_X(col,:)),...
                'shadow_time_ind',shadow_time_ind,'MappingMode',true,'material',MData.material,...
                'depth_transition',theta(2),'k_dry_std_lower',vars_init(8),'material_lower',MData.material_lower);
        end
        Temps_Obs = Data_UAV_X(col,:);
        Temps_Obs = Temps_Obs(:);
        opts = optimoptions('surrogateopt','InitialPoints',MData.minit,'UseParallel',Mdata.parallel,'MaxFunctionEvaluations',MData.nstep);
        Obj = @(theta) sum((Temps_Obs-tima_formod_subset(theta,UAV_flight_ind,formod)).^2./MData.erf(Temps_Obs).^2)/(length(Temps_Obs)-MData.nvars); %Reduced Chi_v         
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