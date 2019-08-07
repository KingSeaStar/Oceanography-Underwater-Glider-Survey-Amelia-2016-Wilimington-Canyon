%
% Divide raw measurements into downcasts and upcasts
% based on haixing_down_up_cast.m
% Daniel Haixing Wang
% June 2019


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('amelia_2016_whole_deployment_CTD_raw_data.mat');

%%
z_m_ind = find(~isnan(whole_deployment.z_raw));
z_m = whole_deployment.z_raw(z_m_ind); % in meters
pres_m = whole_deployment.pressure_raw(z_m_ind);
t_m = whole_deployment.ebd_time(z_m_ind); % in seconds
lon_m = whole_deployment.ebd_lon(z_m_ind);
lat_m = whole_deployment.ebd_lat(z_m_ind);

% unaligned temperature and conductivity data. used for further correction
% and calculation
temp_m = whole_deployment.temperature_raw(z_m_ind);
cond_m = whole_deployment.conductivity_raw(z_m_ind);

% data calculated from un-corrected temperature and conductivity data in
% another script. used for comparison with corrected data in the end.
ptemp_m = whole_deployment.potential_temperature_raw(z_m_ind);
ctemp_m = whole_deployment.conservative_temperature_raw(z_m_ind);
salt_m = whole_deployment.salinity_raw(z_m_ind);
saltA_m = whole_deployment.absolute_salinity_raw(z_m_ind);
sigma0_m = whole_deployment.sigma0_raw(z_m_ind);

dz_m = z_m(2:end)-z_m(1:end-1);
dt_m = t_m(2:end)-t_m(1:end-1);
% dt_m = 1;
vz_m = dz_m./dt_m; % one size smaller than z_m

%% thermistor response correction using Fofonoff (1974) method. 2019.07. Haixing
dtemp_m = temp_m(2:end) - temp_m(1:end-1);
dT_dt = dtemp_m./dt_m;
% Johnson et al. 2007: 0.53 s is the value of the nominal time constant for the SBE-41CP thermistor
% John Kerfoot 2019 poster also uses 0.53 s.
% see Johnson et al. 2007 paper, Kerfoot et al. 2019 poster
tao_T = 0.53; % in seconds.
temp_Ture = temp_m;
temp_Ture(2:end) = temp_m(2:end) + tao_T.*dT_dt;




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data from downcast/upcast measurements

% creteria: downward/upward vertical velocity>0.1 m/s, measurement time interval <
% 10 seconds
% downcast
dn_cast_ind1 = intersect(find(vz_m<-0.1), find(dt_m<10));
% upcast
up_cast_ind1 = intersect(find(vz_m>0.1), find(dt_m<10));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find indices of downcast measurements

% method that considers shift in data caused by using v_z and dt_m
% downcast
t_dn1 = 0.5*(t_m(dn_cast_ind1) + t_m(dn_cast_ind1 + 1));
lon_dn1 = 0.5*(lon_m(dn_cast_ind1) + lon_m(dn_cast_ind1 + 1));
lat_dn1 = 0.5*(lat_m(dn_cast_ind1) + lat_m(dn_cast_ind1 + 1));
z_dn1 = 0.5*(z_m(dn_cast_ind1) + z_m(dn_cast_ind1 + 1));
pres_dn1 = 0.5*(pres_m(dn_cast_ind1) + pres_m(dn_cast_ind1 + 1));
temp_dn1 = 0.5*(temp_m(dn_cast_ind1) + temp_m(dn_cast_ind1 + 1));
temp_True_dn1 = 0.5*(temp_Ture(dn_cast_ind1) + temp_Ture(dn_cast_ind1 + 1));
cond_dn1 = 0.5*(cond_m(dn_cast_ind1) + cond_m(dn_cast_ind1 + 1));
ctemp_dn1 = 0.5*(ctemp_m(dn_cast_ind1) + ctemp_m(dn_cast_ind1 + 1));
ptemp_dn1 = 0.5*(ptemp_m(dn_cast_ind1) + ptemp_m(dn_cast_ind1 + 1));
salt_dn1 = 0.5*(salt_m(dn_cast_ind1) + salt_m(dn_cast_ind1 + 1));
saltA_dn1 = 0.5*(saltA_m(dn_cast_ind1) + saltA_m(dn_cast_ind1 + 1));
sigma0_dn1 = 0.5*(sigma0_m(dn_cast_ind1) + sigma0_m(dn_cast_ind1 + 1));

% upcast
t_up1 = 0.5*(t_m(up_cast_ind1) + t_m(up_cast_ind1 + 1));
z_up1 = 0.5*(z_m(up_cast_ind1) + z_m(up_cast_ind1 + 1));
lon_up1 = 0.5*(lon_m(up_cast_ind1) + lon_m(up_cast_ind1 + 1));
lat_up1 = 0.5*(lat_m(up_cast_ind1) + lat_m(up_cast_ind1 + 1));
pres_up1 = 0.5*(pres_m(up_cast_ind1) + pres_m(up_cast_ind1 + 1));
temp_up1 = 0.5*(temp_m(up_cast_ind1) + temp_m(up_cast_ind1 + 1));
temp_True_up1 = 0.5*(temp_Ture(up_cast_ind1) + temp_Ture(up_cast_ind1 + 1));
cond_up1 = 0.5*(cond_m(up_cast_ind1) + cond_m(up_cast_ind1 + 1));
ctemp_up1 = 0.5*(ctemp_m(up_cast_ind1) + ctemp_m(up_cast_ind1 + 1));
ptemp_up1 = 0.5*(ptemp_m(up_cast_ind1) + ptemp_m(up_cast_ind1 + 1));
salt_up1 = 0.5*(salt_m(up_cast_ind1) + salt_m(up_cast_ind1 + 1));
saltA_up1 = 0.5*(saltA_m(up_cast_ind1) + saltA_m(up_cast_ind1 + 1));
sigma0_up1 = 0.5*(sigma0_m(up_cast_ind1) + sigma0_m(up_cast_ind1 + 1));


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find local low points

zind_low1 = intersect(find(dz_m(1:end-1)<0),find(dz_m(2:end)>0))+1;
t_low1 = t_m(zind_low1);
z_low1=z_m(zind_low1);

% zind_low2 = find(z_low1<-30&z_low1>-350);
% choose -20 m by eyeball (plot(z_m)), then double check plot(filtered_upcast(ii).t,filtered_upcast(ii).z)
zind_low2 = find(z_low1<-20 & z_low1>-360);
t_low2 = t_low1(zind_low2);
z_low2=z_low1(zind_low2);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% store data of each downcast

n_iter = size(t_low2,1);
for ii=1:n_iter
    if ii==1
        raw_downcast(ii).tind = find(t_dn1<=t_low2(ii));
    else
        raw_downcast(ii).tind = find(t_dn1<=t_low2(ii) & t_dn1>t_low2(ii-1));
    end
    
    raw_downcast(ii).t = t_dn1(raw_downcast(ii).tind);
    raw_downcast(ii).lon = lon_dn1(raw_downcast(ii).tind);
    raw_downcast(ii).lat = lat_dn1(raw_downcast(ii).tind);
    raw_downcast(ii).z = z_dn1(raw_downcast(ii).tind);
    raw_downcast(ii).pres = pres_dn1(raw_downcast(ii).tind);
    raw_downcast(ii).temp = temp_dn1(raw_downcast(ii).tind);
    raw_downcast(ii).temp_True = temp_True_dn1(raw_downcast(ii).tind);
    raw_downcast(ii).cond = cond_dn1(raw_downcast(ii).tind);
    
    raw_downcast(ii).ctemp = ctemp_dn1(raw_downcast(ii).tind);
    raw_downcast(ii).ptemp = ptemp_dn1(raw_downcast(ii).tind);
    raw_downcast(ii).saltA = saltA_dn1(raw_downcast(ii).tind);
    raw_downcast(ii).salt = salt_dn1(raw_downcast(ii).tind);
    raw_downcast(ii).sigma0 = sigma0_dn1(raw_downcast(ii).tind);   
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% store data of each upcast
% I made a mistake in indicing in the March-May version in the haixing_down_up_casts.m
% indices for upcasts were off by 1; for example, 
...Mistake:
...raw_upcast(ii).tind = find(t_up1<=t_low2(ii) & t_up1>t_low2(ii-1));
...Correction:
...raw_upcast(ii).tind = find(t_up1<=t_low2(ii+1) & t_up1>t_low2(ii));

n_iter = size(t_low2,1);
for ii=1:n_iter
    if ii == n_iter
        raw_upcast(ii).tind = find(t_up1>t_low2(ii));
    else
        raw_upcast(ii).tind = find(t_up1<=t_low2(ii+1) & t_up1>t_low2(ii));
    end    
    
    raw_upcast(ii).t = t_up1(raw_upcast(ii).tind);
    raw_upcast(ii).lon = lon_up1(raw_upcast(ii).tind);
    raw_upcast(ii).lat = lat_up1(raw_upcast(ii).tind);
    raw_upcast(ii).z = z_up1(raw_upcast(ii).tind);
    raw_upcast(ii).pres = pres_up1(raw_upcast(ii).tind);
    raw_upcast(ii).temp = temp_up1(raw_upcast(ii).tind);
    raw_upcast(ii).temp_True = temp_True_up1(raw_upcast(ii).tind);
    raw_upcast(ii).cond = cond_up1(raw_upcast(ii).tind);
    
    raw_upcast(ii).ctemp = ctemp_up1(raw_upcast(ii).tind);
    raw_upcast(ii).ptemp = ptemp_up1(raw_upcast(ii).tind);
    raw_upcast(ii).saltA = saltA_up1(raw_upcast(ii).tind);
    raw_upcast(ii).salt = salt_up1(raw_upcast(ii).tind);
    raw_upcast(ii).sigma0 = sigma0_up1(raw_upcast(ii).tind);
    
end


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtering raw data to real downcasts and upcasts 
% filtering conditions: casts with actual measurements > 30 counts 
...for both downcast and upcast during the same yo
    
filtered_index_down = [];
for ii = 1:size(raw_downcast,2)
    if size(raw_downcast(ii).t,1)>30 % select casts with over 30 measurements
        filtered_index_down = [filtered_index_down; ii];
    end
end

%
filtered_index_up = [];
for ii = 1:size(raw_upcast,2)
    if size(raw_upcast(ii).t,1)>30 % select casts with over 30 measurements
        filtered_index_up = [filtered_index_up; ii];
    end
end

filtered_index_pair = intersect(filtered_index_down, filtered_index_up);
filtered_raw_downcast = raw_downcast(filtered_index_pair);
filtered_raw_upcast = raw_upcast(filtered_index_pair);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
strictly_filtered_down_up_pair_indices = [];
n_iter = size(raw_downcast,2);
% n_iter = size(raw_upcast,2); % same size as above

% select pairs of downcasts and upcasts subject to condition:
... at least 30 measurements in both downcast and upcast;
... ratio between the numbers of measurements in downcast and upcast is between (3/4, 4/3);
for ii = 1:n_iter
    ratio_downcast_upcast_number_of_measurements = size(raw_downcast(ii).t,1)/size(raw_upcast(ii).t,1);
    if ratio_downcast_upcast_number_of_measurements >= (2/3) ...
            && ratio_downcast_upcast_number_of_measurements <= (3/2)...
            && size(raw_downcast(ii).t,1)>30 ...
            && size(raw_upcast(ii).t,1)>30 ...
        strictly_filtered_down_up_pair_indices = [strictly_filtered_down_up_pair_indices; ii];
    end
end

strictly_filtered_raw_downcast = raw_downcast(strictly_filtered_down_up_pair_indices);
strictly_filtered_raw_upcast = raw_upcast(strictly_filtered_down_up_pair_indices);
size(strictly_filtered_down_up_pair_indices)

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diagnosis plots

plot(strictly_filtered_down_up_pair_indices,'*b');
hold on
plot(filtered_index_pair, '.r');
hold off

%%
% dignosis_index = 1;
% dignosis_index = 3;
% dignosis_index = 10;
% dignosis_index = 50;
% dignosis_index = 97;
% dignosis_index = 500;

plot(raw_downcast(dignosis_index).t, raw_downcast(dignosis_index).z, '.b')
hold on
plot(raw_upcast(dignosis_index).t, raw_upcast(dignosis_index).z, '.r')
hold off