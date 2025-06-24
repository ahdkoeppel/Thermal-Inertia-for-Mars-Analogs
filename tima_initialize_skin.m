function [temperature_column_K] = tima_initialize_skin(k_dry_std,rho_dry,m,theta_k,T_std,T_deep,surface_temperature_C,dt,layer_size,VWC_column,RH,NDAYS,material,varargin)
%***************
% TIMA_INITIALIZE_SKIN
%   Simple version of Heat tansfer model used to estimate a realistic subsurface 
%   temperatures at start of model. The skin approach ultimately solves the
%   upper boundary condition and balances energy fluxes at each time step
%   by adjusting the skin temperature to control upwelling radiative flux.
%
% Syntax
%   [temperature_column_C] = tima_initialize_skin(k_dry_std,rho_dry,m,theta_k,T_deep,surface_temperature_C,dt,layer_size,VWC_column,RH,NDAYS)
%
% Description
%   Script uses the measured surface temperature for 1 day on repeat as the forcing and then distrubutes heat
%   This function uses observational data and assigned thermophysical properties 
%   to estimate the the temperature profile of the subsurface at the start of the simulation.
%   Assumes homogenous bulk density and bulk dry thermal conductivity with depth
%
% Input Parameters
%       k_dry_std: [W/mK] mean dry soil bulk thermal conductivity (vector)
%       rho_dry: [kg/m^3] mean dry soil bulk density (scalar)
%       m: [unitless] pore-network-connectivity-parameter (scalar)
%       theta_k: [0-1, fraction by volume] conductivity soil moisture
%           inflection point - in theory this should be similar to saturation value or porosity (scalar)
%       T_std: [deg K] Standard temperature; typically 300 (scalar)
%       T_deep: [deg C] bottom layer constant temperature (scalar)
%       surface_temperature_C: [deg C] Soil surface temperature (vector)
%       layer_size: [m] array of thickness of subsurface grid layers (vector)
%       dt: [s] time step (scalar)
%       VWC_column: [fraction by volume] volumetric water content measured at multiple depths (2D vector)
%       RH: [fraction] relative humidity of air (vector)
%       NDAYS: [unitless] # of days to run equilib model for: More = closer to equilib, fewer = faster (scalar)
%       material: ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at the surface (char)
%
% Outputs:
%       temperature_column_K = [deg K] subsurface temperatures through NDAYS days of equilibration (vector size layer_size x time)
%
% Author
%    Ari Koeppel -- Copyright 2023
%   
% See also 
%   Hanks 1992: Good tmperature approximations can be made...even for many nonuniform soils by assuming a uniform thermal diffusivity. 
%   tima_conductivity_model_lu2007.m tima_specific_heat_model_hillel.m
    % ***************    
    p = inputParser;
    p.KeepUnmatched=true;
    p.addParameter('material_lower',"basalt",@ischar);
    p.addParameter('depth_transition',sum(layer_size),@isnumeric);
    p.parse(varargin{:});
    p=p.Results;
    material_lower = p.material_lower;
    depth_transition = p.depth_transition;
    [~,D_z] = min(abs(depth_transition-cumsum(layer_size))); %index of depth for k transition


    % ***************
    % Inputs and constants:
    Soil_Temp_K = surface_temperature_C+273.15;
    Day_Dur = 1440/(dt/60); %#of mins/day div by min per interval
    Depth_Max = sum(layer_size);
    NLAY = length(layer_size);
    k = k_dry_std.*ones(1,NLAY);
    k(NLAY) = k_dry_std;
    rho_H2O = 997; %kg/m^3, does not vary much in T range so assume constant density
    %k_H2O = 0.6096; %from Haigh 2012

    % ***************
    % Initialize temperatures of subsurface layers linearly at start
    temperature_column_K = NaN(NDAYS,Day_Dur,NLAY); %Set up Matrix of Temperatures for each day, each minute, and each layer
    temperature_column_K(:,:,NLAY) = T_deep; %Constant lower boundary temp
    for ind = 1:NLAY-1
        temperature_column_K(1,1,ind) = Soil_Temp_K(1)+(T_deep-Soil_Temp_K(1))*(sum(layer_size(1:ind))-layer_size(ind)/2)/Depth_Max;
        k(ind) = tima_conductivity_model_lu2007(k_dry_std,temperature_column_K(1,1,ind),T_std,VWC_column(1,end),theta_k,m,RH(1),material);
    end
    % ***************
    % Cycle through day 1 NDAYS times using measured surf temp and VWCs to estimate subsurface forcing
    for day = 1:NDAYS %Day loop
        for t = 1:Day_Dur %Minute loop
            %k(1) = 1/(VWC_column(t,1)/k_H2O+(k_dry_std*sqrt(temperature_column_K(day,t,1)/T_deep))^-1);
            for z = 1:NLAY-1 %layer loop
                if t == 1 && day > 1 %First step in day (requires pulling from previous day)
                    rho = rho_dry + rho_H2O*VWC_column(t,z); %H2O dep
                    Cp = tima_specific_heat_model_hillel(rho_dry,rho);
                    k(z) = tima_conductivity_model_lu2007(k_dry_std,temperature_column_K(day-1,end,z),T_std,VWC_column(t,z),theta_k,m,RH(t),material);
                    if z == 1  % explicit layer 1 heat flow 
                        Hdn = 2*k(1)*(Soil_Temp_K(t)-temperature_column_K(day-1,end,1))/layer_size(1);
                        Hup = 2*(temperature_column_K(day-1,end,1)-temperature_column_K(day-1,end,2))/(layer_size(1)/k(1)+layer_size(2)/k(2));
                        dT = dt/(Cp*rho*layer_size(1))*(Hdn - Hup);
                        temperature_column_K(day,t,1) = temperature_column_K(day-1,end,1) + dT;
                    elseif strcmp(material_lower,"ice") && z > D_z
                        k(z:end) = 632./temperature_column_K(day-1,end,z)+0.38-0.00197.*temperature_column_K(day-1,end,z); %Wood 2020/Andersson and Inaba 2005;
                        temperature_column_K(day,t,z) = 273.15;
                    elseif strcmp(material_lower,"ice") && z == D_z
                        k(z) = tima_conductivity_model_lu2007(k_dry_std,temperature_column_K(day-1,end,z),T_std,1,theta_k,m,RH(Day_Dur),material);
                        temperature_column_K(day,t,z)  = T_deep;
                    else
                        %Subsurface Multip Factors (Kieffer, 2013)
                        F1 = 2*dt*k(z)/((Cp*rho_dry*layer_size(z)^2)*(1+layer_size(z+1)/layer_size(z)*k(z)/k(z+1)));
                        F3 = (1+(layer_size(z+1)*k(z))/(layer_size(z)*k(z+1)))/(1+(layer_size(z-1)/layer_size(z)*k(z)/k(z-1)));
                        F2 = -(1 + F3);
                        %Temperature response
                        dT = F1*(temperature_column_K(day-1,end,z+1)+F2*temperature_column_K(day-1,end,z)+F3*temperature_column_K(day-1,end,z-1));
                        temperature_column_K(day,t,z)  = temperature_column_K(day-1,end,z) + dT;
                    end
                elseif t >= 2 %All other steps
                    rho = rho_dry + rho_H2O*VWC_column(t,z); %H2O dep
                    Cp = tima_specific_heat_model_hillel(rho_dry,rho);
                    k(z) = tima_conductivity_model_lu2007(k_dry_std,temperature_column_K(day,t-1,z),T_std,VWC_column(t,z),theta_k,m,RH(t),material);
                    if z == 1 % explicit layer 1 heat flow 
                        Hdn = 2*k(1)*(Soil_Temp_K(t)-temperature_column_K(day,t-1,1))/layer_size(1);
                        Hup = 2*(temperature_column_K(day,t-1,1)-temperature_column_K(day,t-1,2))/(layer_size(1)/k(1)+layer_size(2)/k(2));
                        dT = dt/(Cp*rho*layer_size(1))*(Hdn - Hup);
                        temperature_column_K(day,t,1) = temperature_column_K(day,t-1,1) + dT;
                    elseif strcmp(material_lower,"ice") && z > D_z
                        k(z:end) = 632./temperature_column_K(day,t-1,z)+0.38-0.00197.*temperature_column_K(day,t-1,z); %Wood 2020/Andersson and Inaba 2005;
                        temperature_column_K(day,t,z) = 273.15;
                    elseif strcmp(material_lower,"ice") && z == D_z
                        k(z) = tima_conductivity_model_lu2007(k_dry_std,temperature_column_K(day,t-1,z),T_std,1,theta_k,m,RH(t),material);
                        temperature_column_K(day,t,z)  = T_deep;
                    else
                        %Subsurface Multip Factors (Kieffer, 2013)
                        F1 = 2*dt*k(z)/((Cp*rho_dry*layer_size(z)^2)*(1+layer_size(z+1)/layer_size(z)*k(z)/k(z+1)));
                        F3 = (1+(layer_size(z+1)*k(z))/(layer_size(z)*k(z+1)))/(1+(layer_size(z-1)/layer_size(z)*k(z)/k(z-1)));
                        F2 = -(1 + F3);
                        %Temperature response
                        dT = F1*(temperature_column_K(day,t-1,z+1)+F2*temperature_column_K(day,t-1,z)+F3*temperature_column_K(day,t-1,z-1));
                        temperature_column_K(day,t,z) = temperature_column_K(day,t-1,z)+dT;
                    end
                end                 
            end
        end
    end
end