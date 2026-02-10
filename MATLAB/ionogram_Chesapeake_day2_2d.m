%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% 2D Ionogram script. Implements the Built-in function from PHaRLAP
% Use as comparison to a custom built one.
%
% Author: Jake Dunham
% Date: 2/9/2026 (time management is not my strongest quality)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear workspace

% Useful constants
speed_of_light = 2.99792458e8;

% Input parameters
UT = [2000 9 28 10 0];
R12 = 168.8; % taken from WDC-SILSO

%location and heading
origin_lat = 36.768208;                                          % latitude of the start point of ray
origin_long = -76.287491;                                        % longitude of the start point of ray
end_lat = 32.609081;                                             % latitude of the receiver of ray
end_long = -85.481728;                                            % longitude of the receiver of ray
[target_ground_distance_m,ray_bear] = latlon2raz(end_lat,end_long,origin_lat,origin_long);

target_ground_distance = 1e-3*target_ground_distance_m;


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

% grid parameters (taken from example)
max_range = 10000;      % maximum range for sampling the ionosphere (km)
num_range = 201;        % number of ranges (must be < 2000)
range_inc = max_range ./ (num_range - 1);  % range cell size (km)

start_height = 0 ;      % start height for ionospheric grid (km)
height_inc = 3;         % height increment (km)
num_heights = 500;      % number of  heights (must be < 2000)

clear iri_options
% iri_options.Ne_B0B1_model = 'Bil-2000'; % this is a non-standard setting for 
%                                         % IRI but is used as an example

tic
fprintf('Generating ionospheric grid... ')
[iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
    gen_iono_grid_2d(origin_lat, origin_long, R12, UT, ray_bear, ...
                     max_range, num_range, range_inc, start_height, ...
		             height_inc, num_heights, kp, doppler_flag);
toc

% convert plasma frequency grid to electron density in electrons/cm^3
iono_en_grid   = iono_pf_grid.^2   / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

%% Ray Tracing Loop
% ray tracing parameters
nhops = 1;                   % number of hops to raytrace
elevs = 2:2:80;            % initial ray elevation, taken from example
num_elevs = length(elevs);

freqs = 1:0.1:10; % Ray frequency in MHz (recommended by prof)

%ray filtering 
% target_ground_distance = get_ground_distance(origin_lat,origin_long,end_lat,end_long);
distance_from_receiver_threshold = 50; % km, allowed delta from target_ground_distance (guess)
landed_threshold = 5; % km, allowed distance from ground (0) to be considered landed (guess)

%lists of saved results
% collects all rays close to target's ground distance
valid_ray_heights = [];
freq_of_valid_rays = [];
gr_of_valid_rays = [];
% collects the single valid ray with least group range 
min_ray_heights = [];
freq_of_min_rays = [];
gr_of_min_rays = [];

% main loop
for freq = freqs
    tracing_freqs = freq.*ones(size(elevs)); %create an array of freqs for every tested elevation

    [ray_data, ray_path_data] = raytrace_2d(origin_lat, origin_long, elevs, ...
                                            ray_bear, tracing_freqs, nhops, ...
                                            tol, irregs_flag, iono_en_grid, iono_en_grid_5,  ...
                                            collision_freq, start_height, height_inc, range_inc, irreg);

    % access the right part of ray_data and ray_path_data to have matching
    % mask sizes | ray_data.ground_range 1xsize(elevs) |
    % ray_path_data.height 1xsize(heights along path)
    ray_ranges = arrayfun(@(p) p.ground_range(end), ray_path_data);
    ending_heights = arrayfun(@(p) p.height(end), ray_path_data);

    close_enough = abs(ray_ranges - target_ground_distance) <= distance_from_receiver_threshold;
    landed = abs(ending_heights) <= landed_threshold;

    keep = close_enough & landed;

    kept_rays = ray_data(keep);
    kept_ray_path_data = ray_path_data(keep);

    if isempty(kept_rays)
        continue
    end
    
    valid_ray_heights = [valid_ray_heights, kept_rays.virtual_height];
    freq_of_valid_rays = [freq_of_valid_rays, freq*ones(size(kept_rays))];
    gr_of_valid_rays = [gr_of_valid_rays, kept_rays.group_range];

    [minimum_group_range, idx] = min([kept_rays.group_range]);
    min_ray_heights = [min_ray_heights, kept_rays(idx).virtual_height];
    freq_of_min_rays(end+1) = freq;
    gr_of_min_rays(end+1) = minimum_group_range;
end

fprintf('Target Ground Distance: %.2f km\n', target_ground_distance)

%% Ray Path Sanity check
% plots all rays
figure(1) 
start_range = 0;
start_range_idx = fix(start_range ./ range_inc) + 1;
end_range = target_ground_distance+2*distance_from_receiver_threshold;
end_range_idx = fix(end_range ./ range_inc) + 1;
start_ht = start_height;
start_ht_idx = 1;
end_ht = 597;
end_ht_idx = fix(end_ht ./ height_inc) + 1;
iono_pf_subgrid = iono_pf_grid(start_ht_idx:end_ht_idx, ...
    start_range_idx:end_range_idx);
plot_ray_iono_slice(iono_pf_subgrid, start_range, end_range, range_inc, ...
    start_ht, end_ht, height_inc, ray_path_data, 'color', 'w', 'linewidth', 2);

% plots accepted rays
figure(2)
plot_ray_iono_slice(iono_pf_subgrid, start_range, end_range, range_inc, ...
    start_ht, end_ht, height_inc, kept_ray_path_data, 'color', 'w', 'linewidth', 2);

%% Plot ionograms
% valids
figure(3)
scatter(freq_of_valid_rays,valid_ray_heights)
xlabel('Frequency (MHz)')
ylabel('Virtual Height (km)')

% min group distance
figure(4)
plot(freq_of_min_rays,min_ray_heights)
xlabel('Frequency (MHz)')
ylabel('Virtual Height (km)')

