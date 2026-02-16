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
origin_lat = 36.768208;                                          % latitude of the start point of ray
origin_long = -76.287491;                                        % longitude of the start point of ray
end_lat = 32.609081;                                             % latitude of the receiver of ray
end_long = -85.481728;                                            % longitude of the receiver of ray
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

%lat and long bounds finding and num setting
% ionospheric grid parameters (taken from example)
latlong_padding = 5; % degrees
lat_min = min(origin_lat, end_lat) - latlong_padding;
lat_max = max(origin_lat, end_lat) + latlong_padding;

lon_min = min(origin_long, end_long) - latlong_padding;
lon_max = max(origin_long, end_long) + latlong_padding;

lat_inc = 0.3;
lon_inc = 1.0;

lat_start = lat_min;
lon_start = lon_min;

num_lat = floor((lat_max - lat_min)/lat_inc) + 1;
num_lon = floor((lon_max - lon_min)/lon_inc) + 1;

iono_grid_params = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht];

% geomagnetic grid parameters (only change grid steps (inc) here)
B_lat_inc = 1.0;
B_lon_inc = 1.0;
B_ht_inc = 10;

lat_span = (num_lat - 1) * lat_inc;
lon_span = (num_lon - 1) * lon_inc;
ht_span  = (num_ht - 1)  * ht_inc;

B_num_lat = ceil(lat_span / B_lat_inc) + 1;
B_num_lon = ceil(lon_span / B_lon_inc) + 1;
B_num_ht = ceil(ht_span / B_ht_inc) + 1;

B_lat_start = lat_start;
B_lon_start = lon_start;
B_ht_start = ht_start;

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
elevs = 1:0.25:90;               % initial elevation of rays
freqs = 1; % Ray frequency in MHz (recommended by prof)

%ray filtering 
distance_from_receiver_threshold = 50; % km, allowed delta from target_ground_distance (guess)
landed_threshold = 5; % km, allowed distance from ground (0) to be considered landed (guess)

%Data Collection Helpers
valid_ray_heights_o = [];
valid_ray_heights_x = [];
valid_ray_heights_iso = [];

freq_of_valid_rays_o = [];
freq_of_valid_rays_x = [];
freq_of_valid_rays_iso = [];

gr_of_valid_rays_o = [];
gr_of_valid_rays_x = [];
gr_of_valid_rays_iso = [];

ray_bears = ones(size(elevs))*ray_bear; % initial bearing of rays

iono_grid_params = double(reshape(iono_grid_params, 9, 1));
iono_grid_params([3 6 9]) = round(iono_grid_params([3 6 9]));
geomag_grid_params = double(reshape(geomag_grid_params, 9, 1));
geomag_grid_params([3 6 9]) = round(geomag_grid_params([3 6 9]));

for freq = freqs
    
    tracing_freqs = ones(size(elevs))*freq;

    % O-mode 
    OX_mode = 1;
    
   
    fprintf('Generating %d O-mode rays ...', length(elevs));
    tic
    [ray_data_O, ray_O, ray_state_vec_O] = ...
        raytrace_3d(origin_lat, origin_long, ht_start, elevs, ray_bears, tracing_freqs, ...
                    OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	                collision_freq, iono_grid_params, ...
                    Bx, By, Bz, geomag_grid_params);
	      
    NRT_total_time = toc;
    fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
            [ray_data_O.NRT_elapsed_time], NRT_total_time)
    
    % Acceptance of O-mode ray data based on tolerances
    ray_ranges = [ray_data_O.ground_range];
    startend_heights = [ray_O.height];
    receive_heights = startend_heights(1:2:end); % first element is start height 0, removed

    close_enough = abs(ray_ranges - target_ground_distance_km) <= distance_from_receiver_threshold;
    landed = abs(receive_heights) <= landed_threshold;

    keep = close_enough & landed;

    kept_O_rays = ray_data_O(keep);
    kept_O_ray_path_data = ray_O(keep);

    if isempty(kept_O_rays)
        continue
    end

    valid_ray_heights_o = [valid_ray_heights_o,kept_O_rays.virtual_height];


    % Generate the X mode rays - note in the raytrace_3d call the ionosphere does
    % not need to be passed in again as it is already in memory
    OX_mode = -1;
 
    fprintf('Generating %d X-mode rays ...', length(elevs));
    tic
    [ray_data_X, ray_X, ray_sv_X] = ...
        raytrace_3d(origin_lat, origin_long, ht_start, elevs, ray_bears, tracing_freqs, ...
                    OX_mode, nhops, tol);
    NRT_total_time = toc;
    fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
            [ray_data_X.NRT_elapsed_time], NRT_total_time)

    % Generate the rays for the case where the magnetic field is ignored
    % This will be compared to 2D case as a sanity check for the 3D setup.
    OX_mode = 0;
 
    fprintf('Generating %d ''no-field'' rays ...', length(elevs));
    tic
    [ray_data_N, ray_N, ray_sv_N] = ...
    raytrace_3d(origin_lat, origin_long, ht_start, elevs, ray_bears, tracing_freqs, ...
                  OX_mode, nhops, tol);
    NRT_total_time = toc;
    fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
            [ray_data_N.NRT_elapsed_time], NRT_total_time)
end