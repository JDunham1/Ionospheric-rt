function ground_distance = get_ground_distance(start_lat,start_long,end_lat,end_long)
%GET_GROUND_DISTANCE Haversine calculation for ground distance between two
%locations
%   Used to filter landing sites
%   ref: R. Bullock, “Great circle distances and bearings between two locations.” 2007.

R = 6378.14; % Earth's Radius [km] see ref

% convert deg to rad
phi_a = deg2rad(start_lat);
L_a = deg2rad(start_long);
phi_b = deg2rad(end_lat);
L_b = deg2rad(end_long);

del_phi = phi_a - phi_b;
del_L = L_a - L_b;

arg = sin(del_phi/2)^2 + cos(phi_a)*cos(phi_b)*sin(del_L/2)^2;


theta = 2 * asin(sqrt(arg));

ground_distance = R*theta;
end