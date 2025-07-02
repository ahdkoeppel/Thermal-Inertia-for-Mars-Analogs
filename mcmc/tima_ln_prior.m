function [logp] = tima_ln_prior(FitParams,varargin)
%***************
% TIMA_LN_PRIOR
%    Set up resonable ranges for each test variable. 
%
% SYNTAX
%   [logp] = tima_ln_prior(FitParams)

% varargin to adjust values from defaults as needed:
%   'k_upper_min','k_upper_max','k_lower_min','k_lower_max','depth_min','depth_max',
%   'm_min','m_max','CH_min','CH_max','CE_min','CE_max','theta_k_min','theta_k_max','theta_E_min','theta_E_max')
%   k (upper and lower) = thermal conductivity (W/mK);
%   depth = subsurface depth transition between k_upper and k_lower (m)
%   m = pore-network-connectivity (unitless)
%   CH = sensible heat resistence coefficent (unitless)
%   CE = latent heat resistence coefficient (unitless)
%   theta_k = soil moisture inflection point for thermal conductivity (% by volume)
%   theta_E = soil moisture inflection point for latent heat (% by volume)

% INPUTS
%   FitParams: parameters for fitting (vector)
% 
% OUTPUT 
%   logp = If the simulation is within the ranges, the log-probability is
%   0. Otherwise, it's -infinity (vector)
%
% AUTHOR
%   Ari Koeppel, 2023

    p = inputParser;
    p.addRequired('FitParams',@(x) isvector(x) && isnumeric(x));
    p.addOptional('k_upper_min',0.024, @(x) isvector(x) && isnumeric(x) && (x>0)); 
    p.addOptional('k_upper_max',3.7, @(x) isvector(x) && isnumeric(x) && (x<10) && (x>0)); 
    p.addOptional('k_lower_min',0.024, @(x) isvector(x) && isnumeric(x) && (x>0)); 
    p.addOptional('k_lower_max',3.7, @(x) isvector(x) && isnumeric(x) && (x<10) && (x>0)); 
    p.addOptional('depth_min',0, @(x) isvector(x) && isnumeric(x) && (x>=0));
    p.addOptional('depth_max',5, @(x) isvector(x) && isnumeric(x) && (x<=10) && (x>0));
    p.addOptional('m_min',0.01, @(x) isvector(x) && isnumeric(x) && (x>=0));
    p.addOptional('m_max',1.3, @(x) isvector(x) && isnumeric(x) && (x<=1.3) && (x>0));
    p.addOptional('CH_min',1, @(x) isvector(x) && isnumeric(x) && (x>=0));
    p.addOptional('CH_max',1000, @(x) isvector(x) && isnumeric(x) && (x>0));
    p.addOptional('CE_min',1, @(x) isvector(x) && isnumeric(x) && (x>=0));
    p.addOptional('CE_max',10000, @(x) isvector(x) && isnumeric(x) && (x>0));
    p.addOptional('theta_k_min',0.05, @(x) isvector(x) && isnumeric(x) && (x>=0) && (x<1));
    p.addOptional('theta_k_max',0.9, @(x) isvector(x) && isnumeric(x) && (x<=1) && (x>0));
    p.addOptional('theta_E_min',0.01, @(x) isvector(x) && isnumeric(x) && (x>=0) && (x<1));
    p.addOptional('theta_E_max',0.75, @(x) isvector(x) && isnumeric(x) && (x<=1) && (x>0));
    p.parse(FitParams,varargin{:});
    k_upper_min = p.Results.k_upper_min(:);
    k_upper_max = p.Results.k_upper_max(:);
    k_lower_min = p.Results.k_lower_min(:);
    k_lower_max = p.Results.k_lower_max(:);    
    depth_min = p.Results.depth_min(:);
    depth_max = p.Results.depth_max(:);
    m_min = p.Results.m_min(:);
    m_max = p.Results.m_max(:);
    CH_min = p.Results.CH_min(:);
    CH_max = p.Results.CH_max(:);
    CE_min = p.Results.CE_min(:);
    CE_max = p.Results.CE_max(:);
    theta_k_min = p.Results.theta_k_min(:);
    theta_k_max = p.Results.theta_k_max(:);
    theta_E_min = p.Results.theta_E_min(:);
    theta_E_max = p.Results.theta_E_max(:);
    if length(FitParams)==6
        if (FitParams(1)>k_upper_min) && (FitParams(1)<k_upper_max)...
                && (FitParams(2)>m_min) && (FitParams(2)<m_max)...
                && (FitParams(3)>CH_min) && (FitParams(3)<CH_max) && (FitParams(4)>CE_min)...
                && (FitParams(4)<CE_max) && (FitParams(5)>=theta_k_min) && (FitParams(5)<theta_k_max)...
                && (FitParams(6)>theta_E_min) && (FitParams(6)<theta_E_max)
            logp = 0;
        else
            logp = -Inf;
        end
    elseif length(FitParams)==7
        if (FitParams(1)>k_upper_min) && (FitParams(1)<k_upper_max)...
                && (FitParams(2)>m_min) && (FitParams(2)<m_max)...
                && (FitParams(3)>CH_min) && (FitParams(3)<CH_max) && (FitParams(4)>CE_min)...
                && (FitParams(4)<CE_max) && (FitParams(5)>=theta_k_min) && (FitParams(5)<theta_k_max)...
                && (FitParams(6)>theta_E_min) && (FitParams(6)<theta_E_max)...
                && (FitParams(7)>depth_min) && (FitParams(7)<depth_max)
            logp = 0;
        else
            logp = -Inf;
        end
    else
        if (FitParams(1)>k_upper_min) && (FitParams(1)<k_upper_max)...
                && (FitParams(2)>m_min) && (FitParams(2)<m_max)...
                && (FitParams(3)>CH_min) && (FitParams(3)<CH_max) && (FitParams(4)>CE_min)...
                && (FitParams(4)<CE_max) && (FitParams(5)>=theta_k_min) && (FitParams(5)<theta_k_max)...
                && (FitParams(6)>theta_E_min) && (FitParams(6)<theta_E_max)...
                && (FitParams(8)>k_lower_min) && (FitParams(8)<k_lower_max)...
                && (FitParams(7)>depth_min) && (FitParams(7)<depth_max)
            logp = 0;
        else
            logp = -Inf;
        end
    end
end