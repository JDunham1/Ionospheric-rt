%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% 3D Ionogram script. Implements the 
%
% Author: Jake Dunham
% Date: 2/10/2026 (time management is not my strongest quality)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
clear raytrace_3d


% Input parameters
UT = [2000 9 27 8 0];
R12 = 168.8; % taken from WDC-SILSO

%location and heading
origin_lat = 32.609081;                                          % latitude of the start point of ray
origin_long = -85.481728;                                        % longitude of the start point of ray
origin_ht = 0;
end_lat = 27.796419;                                             % latitude of the receiver of ray
end_long = -97.404129;                                            % longitude of the receiver of ray
[target_ground_distance_m,ray_bear] = latlon2raz(end_lat,end_long,origin_lat,origin_long);

target_ground_distance_km = target_ground_distance_m * 1e-3;

% Doppler Parameters
doppler_flag = 1;            % generate ionosphere 5 minutes later so that
                             % Doppler shift can be calculated
                             
% irregularity parameters
irregs_flag = 0;             % no irregularities - not interested in 
                             % Doppler spread or field aligned irregularities
kp = 0;                      % kp not used as irregs_flag = 0. Set it to a 
                            % dummy value 

% simulation tolerances
tol = [1e-7 .01 10];         % ODE tolerance and min/max step sizes

%% generate ionospheric, geomagnetic and irregularity grids

ht_start = 0;          % start height for ionospheric grid (km)
ht_inc = 2;             % height increment (km)
num_ht = 201;           
lat_start = 20;
lat_inc = 0.3;
num_lat = 101.0;
lon_start= -90;
lon_inc = 1.0;
num_lon = 20.0;
iono_grid_params = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = 10;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_lat_start = lat_start;
B_lat_inc = 1.0;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_lon_inc = 1.0;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc); 
geomag_grid_params = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

tic
fprintf('Generating ionospheric grid... ')
[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_params, geomag_grid_params, doppler_flag);
toc

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

%% Ray Tracing Loop

% ray tracing params
nhops = 1;
elevs = 1:0.2:90;               % initial elevation of rays
freqs = 1:0.1:10; % Ray frequency in MHz (recommended by prof)

%ray filtering 
distance_from_receiver_threshold = 25; % km, allowed delta from target_ground_distance (guess)
landed_threshold = 5; % km, allowed distance from ground (0) to be considered landed (guess)

%Data Collection Helpers
valid_ray_heights_o = [];
valid_ray_heights_x = [];
valid_ray_heights_iso = [];

freq_of_valid_rays_o = [];
freq_of_valid_rays_x = [];
freq_of_valid_rays_iso = [];

ray_bears = ones(size(elevs))*ray_bear; % initial bearing of rays

% iono_grid_params = double(reshape(iono_grid_params, 1, 9));
% iono_grid_params([3 6 9]) = round(iono_grid_params([3 6 9]));
% geomag_grid_params = double(reshape(geomag_grid_params, 1, 9));
% geomag_grid_params([3 6 9]) = round(geomag_grid_params([3 6 9]));

for freq = freqs
    
    tracing_freqs = ones(size(elevs))*freq;

    % O-mode 
    OX_mode = 1;
    args = { ...
            origin_lat, origin_long, origin_ht, ...     % 1-3
            elevs, ray_bears, tracing_freqs, ...        % 4-6
            OX_mode, nhops, tol, ...                    % 7-9
            iono_en_grid, iono_en_grid_5, collision_freq, ...  % 10-12
            iono_grid_params, Bx, By, Bz, geomag_grid_params ... % 13-17
           };

    [ray_data_O, ray_O, ray_state_vec_O] = raytrace_3d(args{:});
    
    % Acceptance of O-mode ray data based on tolerances 
    ray_ranges = [ray_data_O.ground_range];
    end_heights = arrayfun(@(p) p.height(end), ray_O); %for some reason the final ray_dath_data entry gives you all of its height values while the other elevations dont...
   
    close_enough = abs(ray_ranges - target_ground_distance_km) <= distance_from_receiver_threshold;
    if any(close_enough)
        fprintf('\nclose enough O-mode ray at frequency = %f',freq)
    end
    landed = abs(end_heights) <= landed_threshold;
    if any(landed)
        fprintf('\nlanded O-mode ray at frequency = %f',freq)
    end

    keep = close_enough & landed;
    if any(keep)
        fprintf('\nkept O-mode ray at frequency = %f',freq)
    end

    kept_O_rays = ray_data_O(keep);
    kept_O_ray_path_data = ray_O(keep);

    if isempty(kept_O_rays)
        continue
    end

    % virtual height is not returned in the 3D ray structure.
    % Instead, use apogee. This is the physical maximum height the ray
    % reaches. A similar height related metric, though not the same.

    valid_ray_heights_o = [valid_ray_heights_o,kept_O_rays.apogee];
    freq_of_valid_rays_o = [freq_of_valid_rays_o,freq*ones(size(kept_O_rays))];


    % Generate the X mode rays - note in the raytrace_3d call the ionosphere does
    % not need to be passed in again as it is already in memory
    OX_mode = -1;
 
    [ray_data_X, ray_X, ray_sv_X] = ...
        raytrace_3d(origin_lat, origin_long, ht_start, elevs, ray_bears, tracing_freqs, ...
                    OX_mode, nhops, tol);

     % Acceptance of X-mode ray data based on tolerances 
    ray_ranges = [ray_data_X.ground_range];
    end_heights = arrayfun(@(p) p.height(end), ray_X); %for some reason the final ray_dath_data entry gives you all of its height values while the other elevations dont...
   
    close_enough = abs(ray_ranges - target_ground_distance_km) <= distance_from_receiver_threshold;
    if any(close_enough)
        fprintf('\nclose enough X-mode ray at frequency = %f',freq)
    end
    landed = abs(end_heights) <= landed_threshold;
    if any(landed)
        fprintf('\nlanded X-mode ray at frequency = %f',freq)
    end

    keep = close_enough & landed;
    if any(keep)
        fprintf('\nkept X-mode ray at frequency = %f',freq)
    end

    kept_X_rays = ray_data_X(keep);
    kept_X_ray_path_data = ray_X(keep);

    if isempty(kept_X_rays)
        continue
    end
    
    % virtual height is not returned in the 3D ray structure.
    % Instead, use apogee. This is the physical maximum height the ray
    % reaches. A similar height related metric, though not the same.

    valid_ray_heights_x = [valid_ray_heights_x,kept_X_rays.apogee];
    freq_of_valid_rays_x = [freq_of_valid_rays_x,freq*ones(size(kept_X_rays))];

    % Generate the rays for the case where the magnetic field is ignored
    % This will be compared to 2D case as a sanity check for the 3D setup.
    OX_mode = 0;
 
    [ray_data_N, ray_N, ray_sv_N] = ...
    raytrace_3d(origin_lat, origin_long, ht_start, elevs, ray_bears, tracing_freqs, ...
                  OX_mode, nhops, tol);

    % Acceptance of X-mode ray data based on tolerances 
    ray_ranges = [ray_data_N.ground_range];
    end_heights = arrayfun(@(p) p.height(end), ray_N); %for some reason the final ray_dath_data entry gives you all of its height values while the other elevations dont...
   
    close_enough = abs(ray_ranges - target_ground_distance_km) <= distance_from_receiver_threshold;
    if any(close_enough)
        fprintf('\nclose enough isotropic ray at frequency = %f',freq)
    end
    landed = abs(end_heights) <= landed_threshold;
    if any(landed)
        fprintf('\nlanded isotropic ray at frequency = %f',freq)
    end

    keep = close_enough & landed;
    if any(keep)
        fprintf('\nkept isotropic ray at frequency = %f',freq)
    end

    kept_iso_rays = ray_data_N(keep);
    kept_iso_ray_path_data = ray_N(keep);

    if isempty(kept_iso_rays)
        continue
    end
    
    % virtual height is not returned in the 3D ray structure.
    % Instead, use apogee. This is the physical maximum height the ray
    % reaches. A similar height related metric, though not the same.

    valid_ray_heights_iso = [valid_ray_heights_iso,kept_iso_rays.apogee];
    freq_of_valid_rays_iso = [freq_of_valid_rays_iso,freq*ones(size(kept_iso_rays))];
end

%% Plot ionograms
% valids
figure
scatter(freq_of_valid_rays_o,valid_ray_heights_o,'DisplayName','O-Mode')
hold on
scatter(freq_of_valid_rays_x,valid_ray_heights_x,'DisplayName','X-Mode')
legend
hold off
xlabel('Frequency (MHz)')
ylabel('Apogee (km)')
title(sprintf('Auburn -> Chesapeake | 9/27/2000 08:00 UTC'))

%sanity check
figure
scatter(freq_of_valid_rays_iso,valid_ray_heights_iso)
xlabel('Frequency (MHz)')
ylabel('Apogee (km)')
title(sprintf('Auburn -> Chesapeake | 9/27/2000 08:00 UTC\nIsotropic 3D Ray for Comparison to 2D'))

