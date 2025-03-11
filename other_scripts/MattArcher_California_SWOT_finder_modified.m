%% Data Processing
%% Load and Subset SWOT data
warning('Obsolete as of Aug 15, 2024')
clear;clc
% Get file list
% files = dir('/mnt/flow/swot/KaRIn/SWOT_L2_LR_SSH_1.0/*Basic_00*.nc'); % list all netCDF files in SWOT directory
files = dir('/mnt/flow/swot/KaRIn/SWOT_L2_LR_SSH_2.0_california_xover/*Basic*.nc');
%
lenfiles = length(files); % get # of files
% Define California Current System (CCS) area
ln = 50; ls = 30; le = -116; lw = -133; % lat/lon limits for CCS
count=0;
for i = 1:lenfiles % foreach file (i.e. for each SWOT pass)
  disp([num2str(i) ' out of ' num2str(lenfiles)]) % Display progress
  fname = [files(i).folder '/' files(i).name]; % Construct filename
  % Load in Lon/Lat data
  lon1 = ncread(fname,'longitude');lon1(lon1>180) = lon1(lon1>180)-360;
  lat1 = ncread(fname,'latitude');
  % Conditional indexing - find only CCS region
  cswath = 35; % get center track of swath
  lonc = lon1(cswath,:); latc = lat1(cswath,:); % reduce lon/lat to a 1D track for indexing purposes
  indc = find(lonc >= lw & lonc <= le & latc >= ls & latc <= ln); % index to retain only swot data in CCS
  if length(indc) > 1 % if there's data in region
    count = count + 1;
    % Subset to region (creates a cell array)
    name{count} = files(i).name;
    time{count} = ncread(fname,'time',indc(1),length(indc))/60/60/24 + datenum(2000,1,1);
    lon{count} = lon1(:,indc);
    lat{count} = lat1(:,indc);
    % Load in and save only key variables from SWOT files
    ssha{count} = ncread(fname,'ssha_karin_2',[1 indc(1)],[size(lon1,1) length(indc)]);
    flag{count} = ncread(fname,'ssha_karin_2_qual',[1 indc(1)],[size(lon1,1) length(indc)]);
    mss{count} = ncread(fname,'mean_sea_surface_cnescls',[1 indc(1)],[size(lon1,1) length(indc)]);
    ith{count} = ncread(fname,'internal_tide_hret',[1 indc(1)],[size(lon1,1) length(indc)]);
    d2c{count} = ncread(fname,'distance_to_coast',[1 indc(1)],[size(lon1,1) length(indc)]);

    % Added by LK
    lon_full{count} = lon1;
    lat_full{count} = lat1;
    time_full{count} = ncread(fname,'time')/60/60/24 + datenum(2000,1,1);
    ssha_full{count} = ncread(fname,'ssha_karin_2');
    flag_full{count} = ncread(fname,'ssha_karin_2_qual');
  end
end
save('/home/kachelein/MATLAB/08.2024/SWOTdata_CCS.mat', 'name', 'time', 'lon', 'lat', 'ssha', 'flag', 'mss', 'ith', 'd2c', 'ssha_full', 'lat_full', 'lon_full', 'flag_full', 'time_full');
%%
error('Forced error.')
%% Load in subsetted data, perform ad-hoc quality control, and plot it up
clear;clc
warning off
% Load subsetted SWOT data
load SWOTdata_CCS
% Load coast data
load coastlines
% Define California Current System (CCS) area
ln = 45; ls = 25; le = -110; lw = -140; % lat/lon limits for CCS
figure % initialize figure window
%
for i = 1:length(time) % for each pass
  disp([num2str(i) ' out of ' num2str(length(time))]) % Display progress
  % Grab variables at time i
  sshai = ssha{i};
  d2ca = d2c{i};
  flaga = flag{i};
  lona = lon{i};
  lata = lat{i};
  % Ad-hoc Quality Control
  sshai(d2ca > 10000000 | d2ca <= 1000 | isnan(d2ca) ) = NaN; % Remove data near coast
  sshai(sshai > 3 | sshai < -3) = NaN; % Remove data with unrealistic amplitudes
  sshai(flaga > 200) = NaN; % ssha(flaga~=0) = NaN; % Use flag variable to remove bad data
  sshf = detrend(sshai','omitnan')'; % detrend to remove long wavelength error
  % Plot coastlines
  plot(coastlon,coastlat),hold on
  % Plot SWOT swath pass
  pcolor(lona,lata,sshf), shading flat, hold on
  caxis([-.35 .35])
  % Clean figure
  set(gcf,'color','w'),set(gca,'Layer','Top','fontsize',20), box on, grid on
  xlim([lw le]); ylim([ls ln])
  Y=get(gca,'ylim');set(gca,'dataaspectratio',[1 cos((Y(2)+Y(1))/2/180*pi) 1])
end
