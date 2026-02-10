function bearing = get_bearing(start_lat,start_long,end_lat,end_long)
%GET_BEARING converts given latitudes and longitudes for starting and
%ending point into a azimuth bearing for use in PHaRLAP

% Bearing is found using great-circle bearing
% ref: R. Bullock, “Great circle distances and bearings between two locations.” 2007.
  
%convert to radians
phi_A = deg2rad(start_lat);
L_A = deg2rad(start_long);
phi_B = deg2rad(end_lat);
L_B = deg2rad(end_long);

del_L = L_B - L_A;

S = cos(phi_B)*sin(del_L);
C = cos(phi_A)*sin(phi_B) - sin(phi_A)*cos(phi_B)*cos(del_L);

beta = atan2(S,C); % bearing in rad
bearing = mod(rad2deg(beta) + 360, 360); % places in range 0-360 degrees for PHaRLAP
end

