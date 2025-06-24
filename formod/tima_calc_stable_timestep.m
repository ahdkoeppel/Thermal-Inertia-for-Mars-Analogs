function [dt] = tima_calc_stable_timestep(C,dz,diffusivity)   
    dt = C.*dz.^2./diffusivity;
end