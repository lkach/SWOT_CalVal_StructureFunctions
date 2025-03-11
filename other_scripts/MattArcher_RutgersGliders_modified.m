%% Data Analysis
%% Revisiting this glider processing work in 2024 now that Luke is using the data
%% I will try this code on the new glider data to see how it looks.
%% LK: modifed further (Jun7,2024) to enable examining multiple integration depths
% Download all glider data, and produce a single netCDF trajectory file

% Load in glider datasets
% Cut to time period of interest
% Quality control procedures
% Compute steric height
% Save new netCDF with variables labelled as per SIO mooring format

clear;clc
warning off

% --------------------------------------------> Preparations for RUTGERS

% Path to file
% path2glider = '/Users/archer/Documents/SWOT/glider/2024_Luke/DATA/'; % laptop
% path2glider = '/mnt/flow/swot/calval/Data/calval_Rutgers_gliders/'; % EDDY
path2glider = '/mnt/flow/swot/Analysis_Luke/Gliders/DATA/'; % EDDY redo with new data

% Glider File
ru32 = 'ru32-20230330T1626-profile-sci-delayed.nc';
ru38 = 'ru38-20230420T1602-profile-sci-delayed.nc';

gx = [ru32; ru38];

% Pick an index time to cut off
% >>> manually/visually chosen: figure;plot(plon,plat,'.') - and - figure;plot(plat,'.')
% OLD DATASET indstart = [37216+574; 34278]; % these indices are for when the glider arrives at it's first station-keeping location
indstart = [433165; 733907]; % these indices are for when the glider arrives at it's first station-keeping location
indend = [5531020; 8283180];
% NOTE: may want to start earlier since there's more data that could be of interest within CalVal region

% Load in SIO time mean RHO profile (provided by Luke K.)
% OLD DATASET SIOdata = '/Users/archer/Documents/SWOT/moorings/JPL_QC/SWOTPOSTLAUNCH_L2_JPLQC_MOORING-P1_CTD-PROFILER_START20230220_END20230713_RT_VER002_5s4s_alltimes.nc';
% SIOdata = '/Users/archer/Documents/SWOT/glider/2024_Luke/SWOTPOSTLAUNCH_L2_JPLQC_MOORING-P1_CTD-PROFILER-BOTH_START20230218_END20230918_DM_VER01_5s4s_alltimes.nc';
SIOdata = ['/mnt/flow/swot/Analysis_Luke/Moorings2/DATA/PROF/' ...
    'SWOTPOSTLAUNCH_L2_JPLQC_MOORING-P1_CTD-PROFILER-BOTH_START20230218_END20230918_DM_VER01_5s4s_alltimes.nc']; % new EDDY
sio_rho_depth = ncread(SIOdata,'RHO_TIMEMEAN_DEPTH');
sio_rho_mean = ncread(SIOdata,'RHO_TIMEMEAN');
% Add deeper average profiles to enable deeper glider dives:
SIOdata2 = ['/mnt/flow/swot/Analysis_Luke/Moorings2/DATA/FIXED/' ...
    'SWOTPOSTLAUNCH_L2_JPLQC_MOORING-S1_CTD-FIXED_START20230303_END20230917_DM_VER01.nc']; % new EDDY
sio_rho_depth2 = ncread(SIOdata2,'RHO_TIMEMEAN_DEPTH');
sio_rho_mean2 = ncread(SIOdata2,'RHO_TIMEMEAN');
% Combine and eliminate NaNs:
sio_rho_depth = [sio_rho_depth; sio_rho_depth2(sio_rho_depth2>500)];
sio_rho_mean  = [sio_rho_mean;  sio_rho_mean2( sio_rho_depth2>500)];
sio_rho_depth = sio_rho_depth(isfinite(sio_rho_mean));
sio_rho_mean = sio_rho_mean(isfinite(sio_rho_mean));
% figure;plot(sio_rho_mean,-sio_rho_depth); % Check plot

% --------------------------------------------> Load in glider, prepare, and save
for k = 1:size(gx,1)

    eval(['fname = [''' path2glider gx(k,:) '''];'])
    % ncdisp(fname)

    % ------------------------------------------------------> Load in variables
    time = ncread(fname,'time')./86400 + datenum(1970,1,1);
    lon = ncread(fname,'longitude');
    lat = ncread(fname,'latitude');
    depth = ncread(fname,'depth');
    temp = ncread(fname,'temperature');
    sal = ncread(fname,'salinity');
    dens = ncread(fname,'density');
    pres = ncread(fname,'pressure');
    cond = ncread(fname,'conductivity');
    %
    ptime = ncread(fname,'profile_time')./86400 + datenum(1970,1,1);
    plon = ncread(fname,'profile_lon');
    plat = ncread(fname,'profile_lat');
    %
    dive = ptime; %repmat(ptime,1,length(time))';

    % ----------------------> Remove indexed times (see above for selection)
    index = [1:indstart(k) indend(k):length(time)];
    time(index) = [];
    lon(index) = [];
    lat(index) = [];
    depth(index) = [];
    pres(index) = [];
    temp(index) = [];
    sal(index) = [];
    dens(index) = [];
    dive(index) = [];
    cond(index) = [];
    %
    ptime(index) = [];
    plon(index) = [];
    plat(index) = [];

    % ----------------------------------------> Catch repeated time records
    [~,IA] = unique(time); % when I checked values, repeated time also repeated P,S,T values (no NaNs)
    %
    depth = depth(IA);
    pres = pres(IA);
    temp = temp(IA);
    cond = cond(IA);
    time = time(IA);
    lon = lon(IA);
    lat = lat(IA);
    plon = plon(IA);
    plat = plat(IA);
    dive = dive(IA);
    whos depth pres temp cond time
    sal = sal(IA);

    % ----------------------------------------> Catch known record-delay errors
    % (when glider pressure, temperature and salinity sensors saved data at
    %  different pressures)

    % -> 1-D interpolation

    % Get index where any of - pres/temp/sal - have values
    indexi = zeros(size(time));
    indexi(~isnan(pres) | ~isnan(temp) | ~isnan(sal)) = 1; %indexi = indexi+1;
    % Testing via David Aragon suggestions (did not work)
    %     indexi(~isnan(pres) | ~isnan(temp)) = 1; %indexi = indexi+1;
    %     indexi(~isnan(pres)) = 1; %indexi = indexi+1;

    % To smooth at locations where there's partial data (i.e. where interpolated, because it's noisier)
    indexs = zeros(size(time));
    sumnnan = ~isnan(pres) + ~isnan(temp) + ~isnan(sal);
    indexs(sumnnan < 3 & sumnnan > 0) = 1; % i.e. where at least one value
    indexs(sumnnan == 0) = [];

    % Piecewise Cubic Hermite Interpolating Polynomial
    mchoice = 'pchip';  % mchoice = 'linear';
    presi = interp1(time(~isnan(pres)),pres(~isnan(pres)),time(indexi==1),mchoice);
    depthi = interp1(time(~isnan(depth)),depth(~isnan(depth)),time(indexi==1),mchoice);
    tempi = interp1(time(~isnan(temp)),temp(~isnan(temp)),time(indexi==1),mchoice);
    condi = interp1(time(~isnan(cond)),cond(~isnan(cond)),time(indexi==1),mchoice);
    %
    loni = interp1(time(~isnan(lon)),lon(~isnan(lon)),time(indexi==1),mchoice);
    lati = interp1(time(~isnan(lat)),lat(~isnan(lat)),time(indexi==1),mchoice);
    ploni = interp1(time(~isnan(plon)),plon(~isnan(plon)),time(indexi==1),mchoice);
    plati = interp1(time(~isnan(plat)),plat(~isnan(plat)),time(indexi==1),mchoice);
    %
    timei = time(indexi==1);
    divei = dive(indexi==1);

    % LK remove everything where p<-1.5 for gsw functions to work:
    pthresh = -1.5;
    depthi = depthi(presi>pthresh);
    tempi = tempi(presi>pthresh);
    condi = condi(presi>pthresh);
    loni = loni(presi>pthresh);
    lati = lati(presi>pthresh);
    ploni = ploni(presi>pthresh);
    plati = plati(presi>pthresh);
    timei = timei(presi>pthresh);
    divei = divei(presi>pthresh);
    indexs = indexs(presi>pthresh);
    presi = presi(presi>pthresh);

    % Re-compute Salinity
    sali = gsw_SP_from_C((condi.*10),tempi,presi); % x10 for correct units % MA

    % Filter
    salis = vfilt(sali,5);
    salis(indexs==0) = sali(indexs==0);
    %     keyboard
    %     % Plot to Check
    %     figure;plot(timei,sali,'.'),hold on,plot(timei,salis,'.');hold on,plot(time,sal,'.')

    % -------------------------------------------------->  Quality Control Data

    qcFlag = zeros(size(timei)); % 0 = good data / 99 = bad data (so starting position is all GOOD)

    % >>>>---------------> Ginput remove any blankly bad data periods
    sal1 = salis;

    %     % First time only:
    %     figure('Position',[441 379 2060 937])
    %     scatter(1:length(timei),-depthi,10,sal1,'filled')
    %     %
    %     figure('Position',[441 379 2060 937])
    %     scatter(1:length(timei),-depthi,10,tempi,'filled')
    %
    % [x,y]=ginput(2); round(x)

    % Remove bad period
    if str2double(gx(k,3:4)) == 32
        %         sal1(265102:276292)=NaN;
        %         qcFlag(265102:276292)=99;
        %         %
        sal1(186707:197017)=NaN;
        qcFlag(186707:197017)=99;
    elseif str2double(gx(k,3:4)) == 38
    end

    % >>>>---------------> Full sample STD-QC
    % Create depth vector
    dvec = 0:5:max(depthi);
    sal2 = sal1;
    temp2 = tempi;

    if str2double(gx(k,3:4)) == 32
        fac = 3; % number of STDs to set as threshold
        facsaveGLOBAL = fac;
    elseif str2double(gx(k,3:4)) == 38
        fac = 3.5; % number of STDs to set as threshold
        facsaveGLOBAL = fac;
    end

    for dd = 1:length(dvec)-1

        indd = find(depthi >= dvec(dd) & depthi < dvec(dd+1));

        % Salinity
        data = detrend(sal1(indd),'omitnan');
        salm = nanmedian(data);
        sals = nanstd(data);
        indqc = find(data < salm-fac.*sals | data > salm+fac.*sals);
        %     [~,indqc] = rmoutliers(data);
        sal2(indd(indqc)) = NaN;
        qcFlag(indd(indqc)) = 99;
        clear data indqc salm sals

        % Temperature
        datat = detrend(tempi(indd),'omitnan');
        tmed = nanmedian(datat);
        tstd = nanstd(datat);
        indqct = find(datat < tmed-fac.*tstd | datat > tmed+fac.*tstd);
        temp2(indd(indqct)) = NaN;
        qcFlag(indd(indqct)) = 99;
        clear datat indqct tmed tstd

    end

    % >>>>---------------> For dt-period STD-QC
    sal3 = sal2; % initialize new data variable
    temp3 = temp2;
    dt = 2; % dt days window
    datett = timei(1); % starting time; % datett = datenum(2023,4,17);

    if str2double(gx(k,3:4)) == 32
        fac = 3.75; % fac*STD to threshold outliers
        facsaveWINDOW = fac;
    elseif str2double(gx(k,3:4)) == 38
        fac = 4; % number of STDs to set as threshold
        facsaveWINDOW = fac;
    end

    while datett < timei(end)

        % Window in time
        indtime = (find(timei >= datett & timei < datett + dt));
        datett = datett + dt/4; % step forward by dt/x
        %
        tempt = temp3(indtime);
        salt = sal3(indtime);
        deptht = depthi(indtime);
        qct = qcFlag(indtime);

        for dd = 1:length(dvec)-1 % ie each depth bin
            indd = find(deptht >= dvec(dd) & deptht < dvec(dd+1));

            % Salinity
            %         data = detrend(salt(indd),'omitnan');
            data = salt(indd);
            salm = nanmedian(data);
            sals = nanstd(data);
            indqc = find(data <= salm-(fac.*sals) | data >= salm+(fac.*sals));
            %         [~,indqc] = rmoutliers(data);
            salt(indd(indqc)) = NaN;
            qct(indd(indqc)) = 99;
            clear data indqc salm sals

            % Temperature
            datat = tempt(indd);
            tmed = nanmedian(datat);
            tstd = nanstd(datat);
            indqct = find(datat <= tmed-(fac.*tstd) | datat >= tmed+(fac.*tstd));
            tempt(indd(indqct)) = NaN;
            qct(indd(indqct)) = 99;
            clear datat indqct tmed tstd

        end

        sal3(indtime) = salt; clear salt
        temp3(indtime) = tempt; clear tempt deptht
        qcFlag(indtime) = qct; clear qct

    end

    salinity = sal3; clear sal1 sal2 sal3
    temperature = temp3; clear temp2 temp3

    % -------------------------------------------> Get density

    % Prepare T/S - Convert using TEOS-10 (McDougall & Barker Toolbox)
    SA = gsw_SA_from_SP(salinity,presi,loni,lati);
    CT = gsw_CT_from_t(SA,temperature,presi);

    % Compute density
    rho = gsw_rho(SA,CT,presi);

    %     figure;scatter(timei,-depthi,10,RHO)

    % -------------------------------------------> Interpolate to uniform z + compute steric height

    rho0 = 1029; % Update: Luke used 1029 % Same as Jinbo used % gsw_rho(35,0,10);
    pXXX = sio_rho_depth; % it is 0:1:500m
    CUTD = [500 600 1000]; % cutd = 500; % cut off depth in meters
    for cutd = CUTD %$
        dz = mode(diff(pXXX));
        udive = unique(divei); % unique dive

        % Initialize Variables
        sh=nan(size(udive)); maxPD=sh;minPD=sh;goodd=sh;
        indsave=sh;midproftime=sh;profdur=sh;
        udcast=sh;mlon=sh;mlat=sh;presmid=sh;
        prof_num=sh;

        errcount=0;threshfail=0;
        for i = 1:length(udive)

            disp([num2str(i) ' out of ' num2str(length(udive))])

            clear divedepth timedepth odh oSA oCT

            index = find(divei == udive(i));
            prof_num(index) = i*ones(size(index));
            indsave(i) = sum(~isnan(rho(index))); % includes full depth profile

            % Threshold (ensure enough data)
            thresh = length(index)*0.8;

            if length(index) > 50 && min(depthi(index)) < 500
                if indsave(i) > thresh

                    divedepth = depthi(index);
                    timedepth = timei(index);
                    p = presi(index);

                    % Save Profile Variables
                    midproftime(i) = nanmean(timedepth(divedepth<=500)); % time at mid of 500m profile
                    profdur(i) = range(timedepth(divedepth<=500))*24*60; % profile duration
                    udcast(i) = mode(diff(divedepth))/abs(mode(diff(divedepth))); % up=-1/down=1 (depth +ve down)
                    mlon(i) = nanmean(loni(index)); % lon of profile
                    mlat(i) = nanmean(lati(index)); % lat of profile
                    presmid(i) = nanmean(p);

                    oRHO = rho(index); % Get density of profile

                    % Interpolate sio_rho_mean to oRHO
                    sio_rhomi = interp1(sio_rho_depth,sio_rho_mean,divedepth);
                    %             figure;plot(sio_rho_mean,-sio_rho_depth),hold on;plot(sio_rhomi,-divedepth,'o')

                    % Compute rho_prime
                    oRHOP = oRHO - sio_rhomi;

                    % Compute steric height from irregular depth
                    if udcast(i) == -1 % if UPcast, need to flip profile

                        %                 [~,IA,~] = unique(p);
                        %                 p2 = p(IA);
                        %                 OS2 = oSA(IA);
                        %                 OC2 = oCT(IA);
                        %                 OR2 = oRHO(IA);
                        %                 ORP2 = oRHOP(IA);

                        index0 = find(divedepth<=cutd); % max depth = cutd meters
                        d2i1 = divedepth(index0); % cut for max depth
                        [d2i,IA,~] = unique(d2i1);
                        r2i= oRHOP(index0(IA)); % cut for max depth
                        ri = interp1(d2i(~isnan(r2i)),r2i(~isnan(r2i)),d2i(~isnan(d2i))); % interpolate small gaps

                        % Compute steric height anomaly
                        %                 sh(i) = squeeze(-1/rho0 .* num_int_trap(flipud(d2i),flipud(ri)));
                        sh(i) = squeeze(-1/rho0 .* num_int_trap((d2i),(ri)));

                        % Save
                        indfind = find(~isnan(r2i));
                        maxPD(i) = max(d2i([indfind(1) indfind(end)])); % max profile depth
                        minPD(i) = min(d2i([indfind(1) indfind(end)])); % min profile depth
                        goodd(i) = length(indfind);   % good data in interpolated profile

                        clear index0 d2i r2i ri indfind

                    elseif udcast(i) == 1

                        index0 = find(divedepth<=cutd); % max depth = cutd meters
                        d2i1 = divedepth(index0); % cut for max depth
                        [d2i,IA,~] = unique(d2i1);
                        r2i= oRHOP(index0(IA)); % cut for max depth
                        ri = interp1(d2i(~isnan(r2i)),r2i(~isnan(r2i)),d2i(~isnan(d2i))); % interpolate small gaps

                        % Compute steric height anomaly
                        sh(i) = squeeze(-1/rho0 .* num_int_trap((d2i),(ri)));

                        % Save
                        indfind = find(~isnan(r2i));
                        maxPD(i) = max(d2i([indfind(1) indfind(end)])); % max profile depth
                        minPD(i) = min(d2i([indfind(1) indfind(end)])); % min profile depth
                        goodd(i) = length(indfind);   % good data in interpolated profile

                        clear index0 d2i r2i ri indfind
                    end
                else
                    threshfail = threshfail+1;
                end
            end
        end

        %     thresh = 0.105;
        %     runber = 8;
        %     %
        %     [ODH,dfo] = despike(sh(udcast==1), thresh, runber);
        %     figure;plot(udive(udcast==1),ODH);
        %     %
        %     [ODH,dfo] = despike(sh(udcast==-1), thresh, runber);
        %     hold on;plot(udive(udcast==-1),ODH)

        % --------------------------------------------------------------> Re-naming

        % Trajectory Data
        TIME = timei;
        LONGITUDE = loni;
        LATITUDE = lati;
        PRES = presi;
        DEPTH = depthi;
        TEMP = tempi;
        CNDC = condi;
        PSAL = salis;
        RHO = rho;
        RHO_0 = rho0;
        QC_FLAG = qcFlag;
        PROFILE_NUM = prof_num;

        % From Luke
        RHO_TIMEMEAN = sio_rho_mean;
        RHO_TIMEMEAN_DEPTH = sio_rho_depth;

        % Profile  data
        TIME_MIDPROFILE = midproftime;
        DURATION_PROFILE = profdur;
        LONGITUDE_PROFILE = mlon;
        LATITUDE_PROFILE = mlat;
        PRESSURE_MIDPROFILE = presmid;
        UPDOWN_CAST = udcast;

        % Steric Height data
        TIME_STERIC_HEIGHT = midproftime;
        STERIC_HEIGHT_ANOMALY = sh;
        STERIC_HEIGHT_MAX_PROF_DEPTH = maxPD;
        STERIC_HEIGHT_MIN_PROF_DEPTH = minPD;
        PROF_GOOD_REALIZATIONS = goodd;
        %
        STERIC_HEIGHT_CUTOFF_DEPTH = cutd;
        STERIC_HEIGHT_MEAN = squeeze(-1/RHO_0 .* num_int_trap(RHO_TIMEMEAN_DEPTH,RHO_TIMEMEAN-RHO_0));

        % -------------------------------------------> Save as a single netCDF file
        if cutd == CUTD(1)
            fileNC = [path2glider gx(k,1:4) '_L2_Processed_' ...
                datestr(TIME(1),'yyyymmddTHH') '_to_' datestr(TIME(end),'yyyymmddTHH') '.nc']; % filename

            % >>---------------Dimensions
            traj = 'trajectory'; trajn = length(TIME);
            trajprof = 'profiles'; trajprofn = length(TIME_MIDPROFILE);

            % >>---------------Initialize TRAJECTORY
            nccreate(fileNC,'TIME','Dimensions',{traj,trajn},'Datatype','double')
            nccreate(fileNC,'LONGITUDE','Dimensions',{traj,trajn},'Datatype','double')
            nccreate(fileNC,'LATITUDE','Dimensions',{traj,trajn},'Datatype','double')
            nccreate(fileNC,'PRES','Dimensions',{traj,trajn},'Datatype','double')
            nccreate(fileNC,'DEPTH','Dimensions',{traj,trajn},'Datatype','double')
            nccreate(fileNC,'TEMP','Dimensions',{traj,trajn},'Datatype','double')
            nccreate(fileNC,'CNDC','Dimensions',{traj,trajn},'Datatype','double')
            nccreate(fileNC,'PSAL','Dimensions',{traj,trajn},'Datatype','double')
            nccreate(fileNC,'RHO','Dimensions',{traj,trajn},'Datatype','double')
            nccreate(fileNC,'RHO_0','Dimensions',{'value',1},'Datatype','double')
            nccreate(fileNC,'QC_FLAG','Dimensions',{traj,trajn},'Datatype','double')
            nccreate(fileNC,'PROFILE_NUM','Dimensions',{traj,trajn},'Datatype','double')

            % >>---------------Initialize TRAJECTORY
            nccreate(fileNC,'TIME_MIDPROFILE','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'DURATION_PROFILE','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'LONGITUDE_PROFILE','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'LATITUDE_PROFILE','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'PRESSURE_MIDPROFILE','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'UPDOWN_CAST','Dimensions',{trajprof,trajprofn},'Datatype','double')

            % >>>---------------Initialize STERIC HEIGHT
            nccreate(fileNC,'TIME_STERIC_HEIGHT','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'STERIC_HEIGHT_ANOMALY','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'STERIC_HEIGHT_MAX_PROF_DEPTH','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'STERIC_HEIGHT_MIN_PROF_DEPTH','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'PROF_GOOD_REALIZATIONS','Dimensions',{trajprof,trajprofn},'Datatype','double')
            % Single Values
            nccreate(fileNC,'STERIC_HEIGHT_CUTOFF_DEPTH','Dimensions',{'value',1},'Datatype','double')

            % >>>---------------Write New Global Attributes
            ncwriteatt(fileNC,'/','experiment_name','SWOT CalVal');
            ncwriteatt(fileNC,'/','title',['AS PER CALVAL MOORING FILES: (1) Trajectory data with initial quality control (e.g. removing repeated time stamps, aligning sampled variables);',...
                newline,'(2) Quality Control Flag variable has full quality control data (removing outliers from detrended global median and sliding window median); and ',...
                newline,'(3) Steric height data (from irregular depth).']);
            ncwriteatt(fileNC,'/','notes',['Trajectory data has been cut in time (to select glider in CalVal region only), and ',...
                newline,'to align sampled variables the trajectory was pchip interpolated between gaps to solve for out-of-phase pressure/temperature/salinity sampling that created many large gaps in salinity. ',...
                newline,'The sigma-threshold for global outlier removal was: ' num2str(facsaveGLOBAL) ', and the sigma-threshold for window outlier removal was: ' num2str(facsaveWINDOW) '.']);
            ncwriteatt(fileNC,'/','creation_date',datestr(now));
            ncwriteatt(fileNC,'/','creator_name','Matthew R. Archer');
            ncwriteatt(fileNC,'/','creator_email','archer@jpl.nasa.gov');
            ncwriteatt(fileNC,'/','institution','Jet Propulsion Laboratory (JPL)');
            ncwriteatt(fileNC,'/','datasource_creator','John Kerfoot (kerfoot@marine.rutgers.edu)');
            ncwriteatt(fileNC,'/','datasource_institution','Rutgers: https://rucool.marine.rutgers.edu');
            ncwriteatt(fileNC,'/','datasource_JPL-EDDY_filename',gx(k,1:end-3));

            % >>>---------------Write Variables + Attributes
            % >>>>>>>>>>>>>>>>>>>>>>>>>>>TRAJECTORY<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ncwrite(fileNC,'TIME',TIME - datenum(1950,1,1,0,0,0))
            ncwriteatt(fileNC,'TIME','units','days since 1950-01-01 00:00:00.0');
            ncwriteatt(fileNC,'TIME','calendar','gregorian');
            ncwriteatt(fileNC,'TIME','long_name','Time of each measurement along trajectory');
            %
            ncwrite(fileNC,'LONGITUDE',LONGITUDE)
            ncwriteatt(fileNC,'LONGITUDE','units','degrees_east');
            ncwriteatt(fileNC,'LONGITUDE','comments','Longitude along trajectory: east positive, west negative, relative to Greenwich meridian')
            %
            ncwrite(fileNC,'LATITUDE',LATITUDE)
            ncwriteatt(fileNC,'LATITUDE','units','degrees_north');
            ncwriteatt(fileNC,'LATITUDE','comments','Latitude along trajectory: north is positive, south is negative, relative to equator')
            %
            ncwrite(fileNC,'DEPTH',DEPTH)
            ncwriteatt(fileNC,'DEPTH','units','meters');
            ncwriteatt(fileNC,'DEPTH','comments','Depth along trajectory: positive downwards')
            %
            ncwrite(fileNC,'PRES',PRES)
            ncwriteatt(fileNC,'PRES','units','dbar');
            ncwriteatt(fileNC,'PRES','comments','Pressure along trajectory: positive downwards')
            %
            ncwrite(fileNC,'CNDC',CNDC)
            ncwriteatt(fileNC,'CNDC','units','S m-1');
            ncwriteatt(fileNC,'CNDC','comments','sea_water_electrical_conductivity along trajectory with initial quality control')
            %
            ncwrite(fileNC,'TEMP',TEMP)
            ncwriteatt(fileNC,'TEMP','units','degree_C');
            ncwriteatt(fileNC,'TEMP','comments','sea_water_temperature along trajectory with initial quality control')
            %
            ncwrite(fileNC,'PSAL',PSAL)
            ncwriteatt(fileNC,'PSAL','units','Practical Salinity PSU');
            ncwriteatt(fileNC,'PSAL','comments','sea_water_practical_salinity along trajectory with initial quality control')
            %
            ncwrite(fileNC,'RHO',RHO)
            ncwriteatt(fileNC,'RHO','units','kg/m^3');
            ncwriteatt(fileNC,'RHO','comments','in-situ density from PSAL, TEMP, and PRES using gsw_toolbox')
            %
            ncwrite(fileNC,'QC_FLAG',QC_FLAG)
            ncwriteatt(fileNC,'QC_FLAG','long_name','Quality control flag');
            ncwriteatt(fileNC,'QC_FLAG','comments','0 = data passed all tests; 99 = data flagged as bad.')
            %
            ncwrite(fileNC,'PROFILE_NUM',PROFILE_NUM)
            ncwriteatt(fileNC,'PROFILE_NUM','comments','Unique number assigned to each vertical profile');
            %
            ncwrite(fileNC,'RHO_0',RHO_0)
            ncwriteatt(fileNC,'RHO_0','units','kg/m^3');
            ncwriteatt(fileNC,'RHO_0','comments','Reference density for calculating steric height');

            % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>PROFILE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ncwrite(fileNC,'TIME_MIDPROFILE',TIME_MIDPROFILE)
            ncwriteatt(fileNC,'TIME_MIDPROFILE','units','days since 1950-01-01 00:00:00.0');
            ncwriteatt(fileNC,'TIME_MIDPROFILE','comments','Mean time along each profile (from surface to STERIC_HEIGHT_CUTOFF_DEPTH)')
            %
            ncwrite(fileNC,'DURATION_PROFILE',DURATION_PROFILE)
            ncwriteatt(fileNC,'DURATION_PROFILE','units','minutes');
            ncwriteatt(fileNC,'DURATION_PROFILE','comments','Time taken to record profile (from surface to STERIC_HEIGHT_CUTOFF_DEPTH)')
            %
            ncwrite(fileNC,'LONGITUDE_PROFILE',LONGITUDE_PROFILE)
            ncwriteatt(fileNC,'LONGITUDE_PROFILE','units','degrees_east');
            ncwriteatt(fileNC,'LONGITUDE_PROFILE','comments','longitude at mid-point of a profile (from surface to STERIC_HEIGHT_CUTOFF_DEPTH)')
            %
            ncwrite(fileNC,'LATITUDE_PROFILE',LATITUDE_PROFILE)
            ncwriteatt(fileNC,'LATITUDE_PROFILE','units','degrees_north');
            ncwriteatt(fileNC,'LATITUDE_PROFILE','comments','latitude at mid-point of a profile (from surface to STERIC_HEIGHT_CUTOFF_DEPTH)')
            %
            ncwrite(fileNC,'PRESSURE_MIDPROFILE',PRESSURE_MIDPROFILE)
            ncwriteatt(fileNC,'PRESSURE_MIDPROFILE','units','dbar');
            %
            ncwrite(fileNC,'UPDOWN_CAST',UPDOWN_CAST)
            ncwriteatt(fileNC,'UPDOWN_CAST','comments','-1 is UPcast, +1 is DOWNcast');

            % >>>>>>>>>>>>>>>>>>>>>>>>>>>STERIC<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ncwrite(fileNC,'TIME_STERIC_HEIGHT',TIME_STERIC_HEIGHT - datenum(1950,1,1,0,0,0))
            ncwriteatt(fileNC,'TIME_STERIC_HEIGHT','units','days since 1950-01-01 00:00:00.0');
            ncwriteatt(fileNC,'TIME_STERIC_HEIGHT','calendar','gregorian');
            %
            ncwrite(fileNC,'STERIC_HEIGHT_ANOMALY',STERIC_HEIGHT_ANOMALY)
            ncwriteatt(fileNC,'STERIC_HEIGHT_ANOMALY','units','m');
            ncwriteatt(fileNC,'STERIC_HEIGHT_ANOMALY','comments','Steric height calculated from density anomaly RHO(trajectory)-RHO_TIMEMEAN(interpolated to trajectory). Trapezoidal integration, dz = variable');
            %
            ncwrite(fileNC,'STERIC_HEIGHT_MAX_PROF_DEPTH',STERIC_HEIGHT_MAX_PROF_DEPTH)
            ncwriteatt(fileNC,'STERIC_HEIGHT_MAX_PROF_DEPTH','units','m');
            ncwriteatt(fileNC,'STERIC_HEIGHT_MAX_PROF_DEPTH','comments','For each profile, the deepest depth with data used to calculate steric height');
            %
            ncwrite(fileNC,'STERIC_HEIGHT_MIN_PROF_DEPTH',STERIC_HEIGHT_MIN_PROF_DEPTH)
            ncwriteatt(fileNC,'STERIC_HEIGHT_MIN_PROF_DEPTH','units','m');
            ncwriteatt(fileNC,'STERIC_HEIGHT_MIN_PROF_DEPTH','comments','For each profile, the shallowest depth with data used to calculate steric height');
            %
            ncwrite(fileNC,'PROF_GOOD_REALIZATIONS',PROF_GOOD_REALIZATIONS)
            ncwriteatt(fileNC,'PROF_GOOD_REALIZATIONS','units','count');
            ncwriteatt(fileNC,'PROF_GOOD_REALIZATIONS','comments','The the number of good data in the profile');
            %
            ncwrite(fileNC,'STERIC_HEIGHT_CUTOFF_DEPTH',STERIC_HEIGHT_CUTOFF_DEPTH)
            ncwriteatt(fileNC,'STERIC_HEIGHT_CUTOFF_DEPTH','units','m');
            ncwriteatt(fileNC,'STERIC_HEIGHT_CUTOFF_DEPTH','comments','The deepest depth considered for this mooring for all steric height calculations.');

        else
            % Added by LK to enable saving just the steric height when using
            % different cutoff depths (600 m and 1000 m):
            fileNC = [path2glider gx(k,1:4) '_L2_Processed_' num2str(cutd) 'm_' ...
                datestr(TIME(1),'yyyymmddTHH') '_to_' datestr(TIME(end),'yyyymmddTHH') '.nc']; % filename

            % >>---------------Dimensions
            trajprof = 'profiles'; trajprofn = length(TIME_MIDPROFILE);

            % >>---------------Initialize TRAJECTORY
            nccreate(fileNC,'RHO_0','Dimensions',{'value',1},'Datatype','double')

            % >>---------------Initialize TRAJECTORY
            nccreate(fileNC,'TIME_MIDPROFILE','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'DURATION_PROFILE','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'LONGITUDE_PROFILE','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'LATITUDE_PROFILE','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'PRESSURE_MIDPROFILE','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'UPDOWN_CAST','Dimensions',{trajprof,trajprofn},'Datatype','double')

            % >>>---------------Initialize STERIC HEIGHT
            nccreate(fileNC,'TIME_STERIC_HEIGHT','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'STERIC_HEIGHT_ANOMALY','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'STERIC_HEIGHT_MAX_PROF_DEPTH','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'STERIC_HEIGHT_MIN_PROF_DEPTH','Dimensions',{trajprof,trajprofn},'Datatype','double')
            nccreate(fileNC,'PROF_GOOD_REALIZATIONS','Dimensions',{trajprof,trajprofn},'Datatype','double')
            % Single Values
            nccreate(fileNC,'STERIC_HEIGHT_CUTOFF_DEPTH','Dimensions',{'value',1},'Datatype','double')

            % >>>---------------Write New Global Attributes
            ncwriteatt(fileNC,'/','experiment_name','SWOT CalVal');
            ncwriteatt(fileNC,'/','title',['AS PER CALVAL MOORING FILES: (1) Trajectory data with initial quality control (e.g. removing repeated time stamps, aligning sampled variables);',...
                newline,'(2) Quality Control Flag variable has full quality control data (removing outliers from detrended global median and sliding window median); and ',...
                newline,'(3) Steric height data (from irregular depth).']);
            ncwriteatt(fileNC,'/','notes',['Trajectory data has been cut in time (to select glider in CalVal region only), and ',...
                newline,'to align sampled variables the trajectory was pchip interpolated between gaps to solve for out-of-phase pressure/temperature/salinity sampling that created many large gaps in salinity. ',...
                newline,'The sigma-threshold for global outlier removal was: ' num2str(facsaveGLOBAL) ', and the sigma-threshold for window outlier removal was: ' num2str(facsaveWINDOW) '.']);
            ncwriteatt(fileNC,'/','creation_date',datestr(now));
            ncwriteatt(fileNC,'/','creator_name','Matthew R. Archer');
            ncwriteatt(fileNC,'/','creator_email','archer@jpl.nasa.gov');
            ncwriteatt(fileNC,'/','institution','Jet Propulsion Laboratory (JPL)');
            ncwriteatt(fileNC,'/','datasource_creator','John Kerfoot (kerfoot@marine.rutgers.edu)');
            ncwriteatt(fileNC,'/','datasource_institution','Rutgers: https://rucool.marine.rutgers.edu');
            ncwriteatt(fileNC,'/','datasource_JPL-EDDY_filename',gx(k,1:end-3));

            % >>>---------------Write Variables + Attributes
            % >>>>>>>>>>>>>>>>>>>>>>>>>>>TRAJECTORY<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ncwrite(fileNC,'RHO_0',RHO_0)
            ncwriteatt(fileNC,'RHO_0','units','kg/m^3');
            ncwriteatt(fileNC,'RHO_0','comments','Reference density for calculating steric height');

            % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>PROFILE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ncwrite(fileNC,'TIME_MIDPROFILE',TIME_MIDPROFILE)
            ncwriteatt(fileNC,'TIME_MIDPROFILE','units','days since 1950-01-01 00:00:00.0');
            ncwriteatt(fileNC,'TIME_MIDPROFILE','comments','Mean time along each profile (from surface to STERIC_HEIGHT_CUTOFF_DEPTH)')
            %
            ncwrite(fileNC,'DURATION_PROFILE',DURATION_PROFILE)
            ncwriteatt(fileNC,'DURATION_PROFILE','units','minutes');
            ncwriteatt(fileNC,'DURATION_PROFILE','comments','Time taken to record profile (from surface to STERIC_HEIGHT_CUTOFF_DEPTH)')
            %
            ncwrite(fileNC,'LONGITUDE_PROFILE',LONGITUDE_PROFILE)
            ncwriteatt(fileNC,'LONGITUDE_PROFILE','units','degrees_east');
            ncwriteatt(fileNC,'LONGITUDE_PROFILE','comments','longitude at mid-point of a profile (from surface to STERIC_HEIGHT_CUTOFF_DEPTH)')
            %
            ncwrite(fileNC,'LATITUDE_PROFILE',LATITUDE_PROFILE)
            ncwriteatt(fileNC,'LATITUDE_PROFILE','units','degrees_north');
            ncwriteatt(fileNC,'LATITUDE_PROFILE','comments','latitude at mid-point of a profile (from surface to STERIC_HEIGHT_CUTOFF_DEPTH)')
            %
            ncwrite(fileNC,'PRESSURE_MIDPROFILE',PRESSURE_MIDPROFILE)
            ncwriteatt(fileNC,'PRESSURE_MIDPROFILE','units','dbar');
            %
            ncwrite(fileNC,'UPDOWN_CAST',UPDOWN_CAST)
            ncwriteatt(fileNC,'UPDOWN_CAST','comments','-1 is UPcast, +1 is DOWNcast');

            % >>>>>>>>>>>>>>>>>>>>>>>>>>>STERIC<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ncwrite(fileNC,'TIME_STERIC_HEIGHT',TIME_STERIC_HEIGHT - datenum(1950,1,1,0,0,0))
            ncwriteatt(fileNC,'TIME_STERIC_HEIGHT','units','days since 1950-01-01 00:00:00.0');
            ncwriteatt(fileNC,'TIME_STERIC_HEIGHT','calendar','gregorian');
            %
            ncwrite(fileNC,'STERIC_HEIGHT_ANOMALY',STERIC_HEIGHT_ANOMALY)
            ncwriteatt(fileNC,'STERIC_HEIGHT_ANOMALY','units','m');
            ncwriteatt(fileNC,'STERIC_HEIGHT_ANOMALY','comments','Steric height calculated from density anomaly RHO(trajectory)-RHO_TIMEMEAN(interpolated to trajectory). Trapezoidal integration, dz = variable');
            %
            ncwrite(fileNC,'STERIC_HEIGHT_MAX_PROF_DEPTH',STERIC_HEIGHT_MAX_PROF_DEPTH)
            ncwriteatt(fileNC,'STERIC_HEIGHT_MAX_PROF_DEPTH','units','m');
            ncwriteatt(fileNC,'STERIC_HEIGHT_MAX_PROF_DEPTH','comments','For each profile, the deepest depth with data used to calculate steric height');
            %
            ncwrite(fileNC,'STERIC_HEIGHT_MIN_PROF_DEPTH',STERIC_HEIGHT_MIN_PROF_DEPTH)
            ncwriteatt(fileNC,'STERIC_HEIGHT_MIN_PROF_DEPTH','units','m');
            ncwriteatt(fileNC,'STERIC_HEIGHT_MIN_PROF_DEPTH','comments','For each profile, the shallowest depth with data used to calculate steric height');
            %
            ncwrite(fileNC,'PROF_GOOD_REALIZATIONS',PROF_GOOD_REALIZATIONS)
            ncwriteatt(fileNC,'PROF_GOOD_REALIZATIONS','units','count');
            ncwriteatt(fileNC,'PROF_GOOD_REALIZATIONS','comments','The the number of good data in the profile');
            %
            ncwrite(fileNC,'STERIC_HEIGHT_CUTOFF_DEPTH',STERIC_HEIGHT_CUTOFF_DEPTH)
            ncwriteatt(fileNC,'STERIC_HEIGHT_CUTOFF_DEPTH','units','m');
            ncwriteatt(fileNC,'STERIC_HEIGHT_CUTOFF_DEPTH','comments','The deepest depth considered for this mooring for all steric height calculations.');


        end
    end
    clearvars -except path2glider gx facsaveGLOBAL facsaveWINDOW k indstart ...
        sio_rho_depth sio_rho_mean indend
end
error('Forced stop (only plotting happens after this line anyway)')
% ---------------------->

%% Checking against previous datasets

clear;clc

fname32_orig = '/Users/archer/Documents/SWOT/glider/ru_gliders/ru32_L2_Processed_20230413T06_to_20230713T15.nc';
fname38_orig = '/Users/archer/Documents/SWOT/glider/ru_gliders/ru38_L2_Processed_20230430T06_to_20230713T13.nc';

fname32_new = '/Users/archer/Documents/SWOT/glider/2024_Luke/DATA/ru32_L2_Processed_20230412T10_to_20230713T15.nc';
fname38_new = '/Users/archer/Documents/SWOT/glider/2024_Luke/DATA/ru38_L2_Processed_20230430T06_to_20230711T23.nc';

% -------------------------------------- Load
ncreadall(fname32_orig)
ncreadall(fname38_orig)

% List of variable names
vars = {'CNDC', 'PRES', 'RHO_0', 'TIME_MIDPROFILE', 'DEPTH', 'PRESSURE_MIDPROFILE', ...
    'STERIC_HEIGHT_ANOMALY', 'TIME_STERIC_HEIGHT', 'DURATION_PROFILE', ...
    'PROFILE_NUM', 'STERIC_HEIGHT_CUTOFF_DEPTH', 'UPDOWN_CAST', 'LATITUDE', ...
    'PROF_GOOD_REALIZATIONS', 'STERIC_HEIGHT_MAX_PROF_DEPTH', 'LATITUDE_PROFILE', ...
    'PSAL', 'STERIC_HEIGHT_MIN_PROF_DEPTH', 'LONGITUDE', 'QC_FLAG', 'TEMP', ...
    'LONGITUDE_PROFILE', 'RHO', 'TIME'};

% Loop through each variable name, create a new variable with '1' appended, and assign the value
for i = 1:length(vars)
    original_var = vars{i};
    new_var = [original_var '1'];
    eval([new_var ' = ' original_var ';']);
end

ncreadall(fname32_new)
ncreadall(fname38_new)

% -------------------------------------- Plotting

figure
plot(TIME_STERIC_HEIGHT + datenum(1950,1,1,0,0,0),STERIC_HEIGHT_ANOMALY), hold on
plot(TIME_STERIC_HEIGHT1 + datenum(1950,1,1,0,0,0),STERIC_HEIGHT_ANOMALY1), hold on
datetick
%
figure
plot(TIME_STERIC_HEIGHT(UPDOWN_CAST==-1) + datenum(1950,1,1,0,0,0),STERIC_HEIGHT_ANOMALY(UPDOWN_CAST==-1)), hold on
plot(TIME_STERIC_HEIGHT1 + datenum(1950,1,1,0,0,0),STERIC_HEIGHT_ANOMALY1), hold on
datetick
%
figure
plot(TIME_STERIC_HEIGHT(UPDOWN_CAST==1) + datenum(1950,1,1,0,0,0),STERIC_HEIGHT_ANOMALY(UPDOWN_CAST==1)), hold on
plot(TIME_STERIC_HEIGHT1 + datenum(1950,1,1,0,0,0),STERIC_HEIGHT_ANOMALY1), hold on
datetick

% Quick QC
thresh = 0.105;
runber = 8;
%
[ODH,dfo] = despike(STERIC_HEIGHT_ANOMALY(UPDOWN_CAST==1), thresh, runber);
%
figure;plot(TIME_STERIC_HEIGHT(UPDOWN_CAST==1) + datenum(1950,1,1,0,0,0),ODH);hold on
plot(TIME_STERIC_HEIGHT1 + datenum(1950,1,1,0,0,0),STERIC_HEIGHT_ANOMALY1), hold on
datetick

% Detrend and demean - downcast
[ODH,dfo] = despike(STERIC_HEIGHT_ANOMALY(UPDOWN_CAST==1), thresh, runber);
figure;plot(TIME_STERIC_HEIGHT(UPDOWN_CAST==1) + datenum(1950,1,1,0,0,0),detrend(ODH-nanmean(ODH),'omitnan'));hold on
plot(TIME_STERIC_HEIGHT1 + datenum(1950,1,1,0,0,0),STERIC_HEIGHT_ANOMALY1-nanmean(STERIC_HEIGHT_ANOMALY1)), hold on
datetick

% Detrend and demean -- upcast
[ODH,dfo] = despike(STERIC_HEIGHT_ANOMALY(UPDOWN_CAST==-1), thresh, runber);
figure;plot(TIME_STERIC_HEIGHT(UPDOWN_CAST==-1) + datenum(1950,1,1,0,0,0),detrend(-ODH-nanmean(-ODH),'omitnan'));hold on
plot(TIME_STERIC_HEIGHT1 + datenum(1950,1,1,0,0,0),STERIC_HEIGHT_ANOMALY1-nanmean(STERIC_HEIGHT_ANOMALY1)), hold on
datetick


% --------------- TEMP/SALT

figure
scatter(TIME(500000:800000),DEPTH(500000:800000),40,PSAL(500000:800000),'filled')

figure
scatter(TIME1(50000:70000),DEPTH1(50000:70000),40,PSAL1(50000:70000),'filled')

