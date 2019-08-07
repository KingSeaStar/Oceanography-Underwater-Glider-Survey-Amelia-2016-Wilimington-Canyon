% Getting raw,  original,  pre-interpolated CTD measurements
% Daniel Haixing Wang
% June 2019

%% Code that is based on Donglai's code: dg_glider_ctd_processing.m

% load('amelia-20160226-maracoos-post.mat')

% initialize variables to load from EBD data structure
time_ebd_raw = []; pressure_raw = []; temperature_raw = []; conductivity_raw = [];

% initialize variables to load from DBD data structure
time_dbd_raw = []; lon_raw = []; lat_raw = []; gps_lon_raw = []; gps_lat_raw = [];

% load glider time,  pressure,  temperature,  and conductivity data from EBD
% Glider pressure in unit of bar,  conductivity in unit of S/m.
for iter = 1:length(gdataebd)
  if ~isempty(gdataebd{iter})
    time_ebd_raw = [time_ebd_raw;gdataebd{iter}(:, slebd.sci_m_present_time)];
    pressure_raw = [pressure_raw;gdataebd{iter}(:, slebd.sci_water_pressure)];
    temperature_raw = [temperature_raw;gdataebd{iter}(:, slebd.sci_water_temp)];
    conductivity_raw = [conductivity_raw;gdataebd{iter}(:, slebd.sci_water_cond)];
  end %if
end %for

% sort the data by time
[time_ebd_raw_sorted, sorting_indices_time_ebd_raw] = sort(time_ebd_raw);
pressure_raw = pressure_raw(sorting_indices_time_ebd_raw);
temperature_raw = temperature_raw(sorting_indices_time_ebd_raw);
conductivity_raw = conductivity_raw(sorting_indices_time_ebd_raw);

% load glider time,  lon,  lat data from DBD
for iter = 1:length(gdatadbd)
  if ~isempty(gdatadbd{iter})
    time_dbd_raw = [time_dbd_raw;gdatadbd{iter}(:, sldbd.m_present_time)];
    lon_raw = [lon_raw;gdatadbd{iter}(:, sldbd.m_lon)];
    lat_raw = [lat_raw;gdatadbd{iter}(:, sldbd.m_lat)];
    gps_lon_raw = [gps_lon_raw;gdatadbd{iter}(:, sldbd.m_gps_lon)];
    gps_lat_raw = [gps_lat_raw;gdatadbd{iter}(:, sldbd.m_gps_lat)];
  end %if
end %for

% sort the data by time
[time_dbd_raw_sorted, sorting_indices_time_dbd_raw] = sort(time_dbd_raw);
lon_raw = lon_raw(sorting_indices_time_dbd_raw);
lat_raw = lat_raw(sorting_indices_time_dbd_raw);
gps_lon_raw = gps_lon_raw(sorting_indices_time_dbd_raw);
gps_lat_raw = gps_lat_raw(sorting_indices_time_dbd_raw);

% convert from glider lon/lat format to decimal degrees lon/lat
[lonDR, latDR] = dg_llg2lld(lon_raw, lat_raw);

% extract not NaN GPS location fixes
nnan_gps = find(~isnan(gps_lon_raw));
[unique_t_gps, unique_t_gps_indices] = unique(time_dbd_raw_sorted);

% extracting the good GPS lon/lat indices from the DBD variables.
good_gps_indices = intersect(nnan_gps, unique_t_gps_indices);

% obtain time of GPS fixes
gps_time = time_dbd_raw_sorted(good_gps_indices);

% convert from glider lon/lat format to decimal degrees lon/lat
[gps_lon,  gps_lat] = dg_llg2lld(gps_lon_raw(good_gps_indices),  gps_lat_raw(good_gps_indices));

% interplate GPS fixes to obtain glider position
gps_lon_lat_complex = gps_lon + gps_lat*i; % use complex numbers to represent GPS lon/lat
lonlat_ebd = interp1(gps_time, gps_lon_lat_complex, time_ebd_raw_sorted); % interpolate to EBD time stamps
lonlat_dbd = interp1(gps_time, gps_lon_lat_complex, time_dbd_raw_sorted); % interplate to DBD time stamps
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% below are used for later calculation, 
...such as divide_raw_data_to_downcasts_and_upcasts.m
lon_ebd = real(lonlat_ebd); % lon matched to EBD time stamps
lat_ebd = imag(lonlat_ebd); % lat matched to EBD time stamps
lon_dbd = real(lonlat_dbd); % lon matched to DBD time stamps
lat_dbd = imag(lonlat_dbd); % lat matched to DBD time stamps
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CTD data QC
bad_ctd_indices = find(conductivity_raw < 0.1); % need to make sure why we use this cretiron
conductivity_raw(bad_ctd_indices) = NaN;
temperature_raw(bad_ctd_indices) = NaN;

% use GSW toolbox to calculate potential temperature,  density,  and salinity
% note GSW toolbox need pressure in unit of dbar and conductivity in unit of mS/cm

%
z_raw = gsw_z_from_p(pressure_raw*10, lat_ebd); % depth,  converting pressure from bar to dbar
salinity_raw = gsw_SP_from_C(conductivity_raw*10, temperature_raw, pressure_raw*10); % salinity,  converting pressure from bar to dbar and conductivity from S/m to mS/cm.
absolute_salinity_raw = gsw_SA_from_SP(salinity_raw, pressure_raw*10, lon_ebd, lat_ebd); % absolute salinity,  converting pressure from bar to dbar
conservative_temperature_raw = gsw_CT_from_t(absolute_salinity_raw, temperature_raw, pressure_raw*10); % conservative temperature,  converting pressure from bar to dbar
potential_temperature_raw = gsw_pt_from_CT(absolute_salinity_raw, conservative_temperature_raw); % potential temperature
rho_raw = gsw_rho(absolute_salinity_raw, conservative_temperature_raw, pressure_raw); % in-situ density
sigma0_raw = gsw_sigma0(absolute_salinity_raw, conservative_temperature_raw); % potential density anomaly

save amelia_2016_glider_CTD_data_correction_raw_data.mat



%% houskeeping. store relevant data into structure for the entire glider deployment


whole_deployment.ebd_time = time_ebd_raw_sorted;
whole_deployment.dbd_time = time_dbd_raw_sorted;
whole_deployment.pressure_raw = pressure_raw;
whole_deployment.z_raw = z_raw;
whole_deployment.ebd_lon = lon_ebd;
whole_deployment.ebd_lat = lat_ebd;
whole_deployment.dbd_lon = lon_dbd;
whole_deployment.dbd_lat = lat_dbd;
whole_deployment.temperature_raw = temperature_raw;
whole_deployment.conductivity_raw = conductivity_raw;
whole_deployment.salinity_raw = salinity_raw;
whole_deployment.absolute_salinity_raw = absolute_salinity_raw;
whole_deployment.conservative_temperature_raw = conservative_temperature_raw;
whole_deployment.potential_temperature_raw = potential_temperature_raw;
whole_deployment.rho_raw = rho_raw;
whole_deployment.sigma0_raw = sigma0_raw;

%%
save('amelia_2016_whole_deployment_CTD_raw_data.mat','whole_deployment') ;





