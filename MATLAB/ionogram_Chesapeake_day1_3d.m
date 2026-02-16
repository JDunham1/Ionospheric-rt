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

% Input parameters
UT = [2000 9 27 8 0];
R12 = 168.8; % taken from WDC-SILSO

%location and heading
origin_lat = 36.768208;                                          % latitude of the start point of ray
origin_long = -76.287491;                                        % longitude of the start point of ray
origin_ht = 0;                                                   % height at start point of ray (assumed ground level like 2D)
end_lat = 32.609081;                                             % latitude of the receiver of ray
end_long = -85.481728;                                            % longitude of the receiver of ray
ray_bear = get_bearing(origin_lat,origin_long,end_lat,end_long); % bearing of rays

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

% ionospheric grid parameters (taken from example)
ht_start = 60;          % start height for ionospheric grid (km)
ht_inc = 2;             % height increment (km)
num_ht = 201;           
lat_start = -20.0;
lat_inc = 0.3;
num_lat = 101.0;
lon_start= 128.0;
lon_inc = 1.0;
num_lon = 5.0;
iono_grid_params = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];

% geomagnetic grid parameters
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

%% Ray Tracing Loop

% O-mode 
OX_mode = 1;
elevs = [3:1:81];               % initial elevation of rays
freqs = ones(size(elevs))*15;   % frequency (MHz)
ray_bears = zeros(size(elevs)); % initial bearing of rays
fprintf('Generating %d O-mode rays ...', num_elevs);
tic
[ray_data_O, ray_O, ray_state_vec_O] = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	          collision_freq, iono_grid_parms, Bx, By, Bz, ...
	          geomag_grid_parms);
	      
NRT_total_time = toc;
fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
        [ray_data_O.NRT_elapsed_time], NRT_total_time)
