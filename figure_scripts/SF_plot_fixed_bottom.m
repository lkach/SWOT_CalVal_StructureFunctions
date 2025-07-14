%% Figure and numerical results: interpolated to fixed depth and at densely spaced times
%% FIGURE # - Structure functions from moorings, gliders and SWOT

%% Load MOORINGS

MP_DIR = dir( '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/PROF/*.nc');
MF_DIR = dir('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/FIXED/*.nc');
% % % Separate files for GPS mooring data:
% SIO_GPS_DIR = dir('/Users/kachelein/Documents/JPL/work/07.2023_QC_moorings/PROFILERS/SWOT*.nc');
% PMEL_GPS_DIR = dir('/Users/kachelein/Documents/JPL/work/09.2024/SWO_hrly/SWO*_hrly.txt');

% % % Separate files for GPS mooring data:
SIO_GPS_DIR = dir('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/GPS_Moorings/SWOT*GPS.nc');
PMEL_GPS_DIR = dir('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/GPS_Gliders/SWO*_hrly.txt');
warning('Before publication, be sure to package the GPS data with the other data (or ideally SIO and PMEL will send us an updated version).')

MM_list = {'P1','P2','P3',    'P5','P6','P7',    'S1','S2','S3','S4'};

warning('Requires custom function ncreadall (like "load" but for NetCDF files)')
warning('"ncreadall" has the following license:')
warning(['For use to meet course requirements and perform academic research at ' ...
         'degree-granting institutions only. It is not available for government, ' ...
         'commercial, or other organizational use, unless granted permission by the original author.'])
disp(' ')
warning(['I, the aforementioned author, grant anyone the permission to use "ncreadall" to ' ...
         'verify the findings of this study for the purposes of investigating replicability ' ...
         'and checking our work. This permission applies to all other functions written by ' ...
         'Luke Kachelein that bear the aforementioned license and are used in this study.'])

% Eliminate P4, which is still present in M_GPS_DIR and MP_DIR:
UNUSABLE_MOORING_IND_P   = [];
UNUSABLE_MOORING_IND_GPS = [];
for mi = 1:length(MP_DIR)
    if contains(MP_DIR(mi).name,'-P4_')
        UNUSABLE_MOORING_IND_P = [UNUSABLE_MOORING_IND_P; mi];
    else
    end
    if contains(SIO_GPS_DIR(mi).name,'-P4_')
        UNUSABLE_MOORING_IND_GPS = [UNUSABLE_MOORING_IND_GPS; mi];
    else
    end
end
MP_DIR    =    MP_DIR(setdiff(1:length(MP_DIR)   , UNUSABLE_MOORING_IND_P  ));
SIO_GPS_DIR = SIO_GPS_DIR(setdiff(1:length(SIO_GPS_DIR), UNUSABLE_MOORING_IND_GPS));


EXCLUDED_VARS_MP = {'DEPTH','TIME', 'RHO',  'CNDC','PRES','PSAL','TEMP','PROF_NUM'};
EXCLUDED_VARS_MF = {'DEPTH',        'RHO',  'CNDC','CNDC_QC','PRES','PRES_QC','PSAL','PSAL_QC','TEMP','TEMP_QC'};

T0 = datenum('1950-01-01');

for mi = 1:length(MM_list)
    MM = MM_list{mi};
    % PROF
    eval([MM '.P = ncreadall(''' MP_DIR(mi).folder '/' MP_DIR(mi).name ''' , EXCLUDED_VARS_MP);']);
    % FIXED
    eval([MM '.F = ncreadall(''' MF_DIR(mi).folder '/' MF_DIR(mi).name ''' , EXCLUDED_VARS_MF);']);
    % GPS locaitons (separately handled):
    if strcmp(MM(1),'S')
        eval([MM '.GPS.TIME_GPS_SURFACE_BUOY = ncread('''      SIO_GPS_DIR(mi).folder '/' SIO_GPS_DIR(mi).name ''',''TIME_GPS_SURFACE_BUOY'');']);
        eval([MM '.GPS.LATITUDE_GPS_SURFACE_BUOY = ncread('''  SIO_GPS_DIR(mi).folder '/' SIO_GPS_DIR(mi).name ''',''LATITUDE_GPS_SURFACE_BUOY'');']);
        eval([MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY = ncread(''' SIO_GPS_DIR(mi).folder '/' SIO_GPS_DIR(mi).name ''',''LONGITUDE_GPS_SURFACE_BUOY'');']);
    elseif strcmp(MM(1),'P')
        disp([MM ' - ' [PMEL_GPS_DIR(mi).folder '/SWO' MM(2) '_hrly.txt']])
        P_GPS_temp = readtable([PMEL_GPS_DIR(mi).folder '/SWO' MM(2) '_hrly.txt']);
        eval([MM '.GPS.TIME_GPS_SURFACE_BUOY = datenum(P_GPS_temp.time) - T0;']);
        eval([MM '.GPS.LATITUDE_GPS_SURFACE_BUOY = P_GPS_temp.latitude;']);
        eval([MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY = P_GPS_temp.longitude;']);
        clear P_GPS_temp
    else
    end
end
% Workspace variables take up 1348470499 bytes of memory = 1.3485 GB of memory.

%% Correct the GPS locations by taking info from .GPS and putting it in .P

interp1nan = @(X,V,Xq) interp1(X(isfinite(X)&isfinite(V)),V(isfinite(X)&isfinite(V)),Xq);
for mi = 1:length(MM_list)
    MM = MM_list{mi};
        eval([replace('#.P.LATITUDE_GPSSB_STERIC_HEIGHT','#',MM) '=' ...
              replace('interp1nan(#.GPS.TIME_GPS_SURFACE_BUOY, #.GPS.LATITUDE_GPS_SURFACE_BUOY, #.P.TIME_STERIC_HEIGHT);','#',MM)]);
        eval([replace('#.P.LONGITUDE_GPSSB_STERIC_HEIGHT','#',MM) '=' ...
              replace('interp1nan(#.GPS.TIME_GPS_SURFACE_BUOY, #.GPS.LONGITUDE_GPS_SURFACE_BUOY, #.P.TIME_STERIC_HEIGHT);','#',MM)]);
end

figure;
plot(P1.P.TIME_STERIC_HEIGHT, P1.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-'); hold on
plot(P2.P.TIME_STERIC_HEIGHT, P2.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P3.P.TIME_STERIC_HEIGHT, P3.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P5.P.TIME_STERIC_HEIGHT, P5.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P6.P.TIME_STERIC_HEIGHT, P6.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P7.P.TIME_STERIC_HEIGHT, P7.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S1.P.TIME_STERIC_HEIGHT, S1.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S2.P.TIME_STERIC_HEIGHT, S2.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S3.P.TIME_STERIC_HEIGHT, S3.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S4.P.TIME_STERIC_HEIGHT, S4.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
figure;
plot(P1.P.TIME_STERIC_HEIGHT, P1.P.LONGITUDE_GPSSB_STERIC_HEIGHT,'.-'); hold on
plot(P2.P.TIME_STERIC_HEIGHT, P2.P.LONGITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P3.P.TIME_STERIC_HEIGHT, P3.P.LONGITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P5.P.TIME_STERIC_HEIGHT, P5.P.LONGITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P6.P.TIME_STERIC_HEIGHT, P6.P.LONGITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P7.P.TIME_STERIC_HEIGHT, P7.P.LONGITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S1.P.TIME_STERIC_HEIGHT, S1.P.LONGITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S2.P.TIME_STERIC_HEIGHT, S2.P.LONGITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S3.P.TIME_STERIC_HEIGHT, S3.P.LONGITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S4.P.TIME_STERIC_HEIGHT, S4.P.LONGITUDE_GPSSB_STERIC_HEIGHT,'.-')
figure;
plot(P1.P.LONGITUDE_GPSSB_STERIC_HEIGHT, P1.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-'); hold on
plot(P2.P.LONGITUDE_GPSSB_STERIC_HEIGHT, P2.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P3.P.LONGITUDE_GPSSB_STERIC_HEIGHT, P3.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P5.P.LONGITUDE_GPSSB_STERIC_HEIGHT, P5.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P6.P.LONGITUDE_GPSSB_STERIC_HEIGHT, P6.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(P7.P.LONGITUDE_GPSSB_STERIC_HEIGHT, P7.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S1.P.LONGITUDE_GPSSB_STERIC_HEIGHT, S1.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S2.P.LONGITUDE_GPSSB_STERIC_HEIGHT, S2.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S3.P.LONGITUDE_GPSSB_STERIC_HEIGHT, S3.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')
plot(S4.P.LONGITUDE_GPSSB_STERIC_HEIGHT, S4.P.LATITUDE_GPSSB_STERIC_HEIGHT,'.-')

%% Load GLIDERS

G_DIR = '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/GLIDERS/';
% % % 
ru32_600 = ncreadall([G_DIR 'ru32_L2_Processed_600m_20230412T10_to_20230713T15.nc']);
ru32_1000 = ncreadall([G_DIR 'ru32_L2_Processed_1000m_20230412T10_to_20230713T15.nc']);
% % % 
ru38_600 = ncreadall([G_DIR 'ru38_L2_Processed_600m_20230430T06_to_20230711T23.nc']);
ru38_1000 = ncreadall([G_DIR 'ru38_L2_Processed_1000m_20230430T06_to_20230711T23.nc']);

Glider_List = {'ru32_600','ru32_1000','ru38_600','ru38_1000'};
for gi = 1:length(Glider_List)
    GG = Glider_List{gi};
    % Correct the steric height discrepancy due to inconsistent reference
    % densities (RHO_0):
    eval([GG '.STERIC_HEIGHT_ANOMALY = ' GG '.STERIC_HEIGHT_ANOMALY*[' GG '.RHO_0/S1.P.RHO_0];']);
    eval([GG '.RHO_0 = ' GG '.RHO_0*[S1.P.RHO_0/' GG '.RHO_0];']);
end

%% Load SWOT

% Load SWOT data:
% SWOT = load('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/SWOTdata_CCS.mat',...
SWOT = load('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/SWOTdata_CCS_33_385.mat',...
... SWOT = load('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/SWOTdata_L3_CCS_33_385.mat',...
            'time','ssha_karin_2','ssha_karin_2_qual','name','lat','lon','internal_tide_hret','height_cor_xover','height_cor_xover_qual');
% ^ Omitted: distance_to_coast, mean_sea_surface_cnescls
% ~800 MB

error('Force stop for loading only')

%% Visualize Profiler depths to determine adequate cutoffs

close all

T0 = datenum('1950-01-01');

% S: Visual adequate profiler total depth cutoffs:
figure
plot(S1.P.STERIC_HEIGHT_MAX_PROF_DEPTH - S1.P.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-'); hold on
plot(S2.P.STERIC_HEIGHT_MAX_PROF_DEPTH - S2.P.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
plot(S3.P.STERIC_HEIGHT_MAX_PROF_DEPTH - S3.P.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
plot(S4.P.STERIC_HEIGHT_MAX_PROF_DEPTH - S4.P.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
% Verdict: >465 for the S profilers

% P: Visual adequate profiler total depth cutoffs:
figure
plot(P1.P.STERIC_HEIGHT_MAX_PROF_DEPTH - P1.P.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-'); hold on
plot(P2.P.STERIC_HEIGHT_MAX_PROF_DEPTH - P2.P.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
plot(P3.P.STERIC_HEIGHT_MAX_PROF_DEPTH - P3.P.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
% plot(P4.P.STERIC_HEIGHT_MAX_PROF_DEPTH - P4.P.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
plot(P5.P.STERIC_HEIGHT_MAX_PROF_DEPTH - P5.P.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
plot(P6.P.STERIC_HEIGHT_MAX_PROF_DEPTH - P6.P.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
plot(P7.P.STERIC_HEIGHT_MAX_PROF_DEPTH - P7.P.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
% Verdict: >410 for the P profilers

% Gliders:
figure
plot(T0+ru32_600.TIME_STERIC_HEIGHT, ru32_600.STERIC_HEIGHT_MAX_PROF_DEPTH - ru32_600.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-'); hold on
plot(T0+ru32_1000.TIME_STERIC_HEIGHT, ru32_1000.STERIC_HEIGHT_MAX_PROF_DEPTH - ru32_1000.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
plot(T0+ru38_600.TIME_STERIC_HEIGHT, ru38_600.STERIC_HEIGHT_MAX_PROF_DEPTH - ru38_600.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
plot(T0+ru38_1000.TIME_STERIC_HEIGHT, ru38_1000.STERIC_HEIGHT_MAX_PROF_DEPTH - ru38_1000.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
title('GLIDERS: MAX-MIN depths')
legend('ru32 600','ru32 1000','ru38 600','ru38 1000')



SP_RANGE_LIM = 477;
PP_RANGE_LIM = 414;
G_RANGE_LIM  = 570; % 570 for the most data to be used, 590 for more restriction

% For reference, here are the nominal depths of S moorings:
% S1.F.DEPTH_NOMINAL'
% 605         804        1205        1814        2523        3432        4427
% S2
% 604         803        1204        1706        2361        3268        4373
% S3
% 605         804        1205        1814        2523        3432        4428
% S4
% 608         798        1199        1813        2522        3432        4327

%% Visualize the depths of the fixed CTDs

close all 
COLOR_ORDER = colororder;
figure
plot(T0 + S1.F.TIME, S1.F.DEPTH_AUGMENTED, '.-', 'Color',COLOR_ORDER(1,:),'HandleVisibility','off'); hold on
plot(T0 + S2.F.TIME, S2.F.DEPTH_AUGMENTED, '.-', 'Color',COLOR_ORDER(2,:),'HandleVisibility','off')
plot(T0 + S3.F.TIME, S3.F.DEPTH_AUGMENTED, '.-', 'Color',COLOR_ORDER(3,:),'HandleVisibility','off')
plot(T0 + S4.F.TIME, S4.F.DEPTH_AUGMENTED, '.-', 'Color',COLOR_ORDER(4,:),'HandleVisibility','off')
plot(nan,nan, '.-', 'Color',COLOR_ORDER(1,:))
plot(nan,nan, '.-', 'Color',COLOR_ORDER(2,:))
plot(nan,nan, '.-', 'Color',COLOR_ORDER(3,:))
plot(nan,nan, '.-', 'Color',COLOR_ORDER(4,:))
legend('S1','S2','S3','S4')
datetick

figure
plot(T0 + P1.F.TIME, P1.F.DEPTH_AUGMENTED, '.-', 'Color',COLOR_ORDER(1,:),'HandleVisibility','off'); hold on
plot(T0 + P2.F.TIME, P2.F.DEPTH_AUGMENTED, '.-', 'Color',COLOR_ORDER(2,:),'HandleVisibility','off')
plot(T0 + P3.F.TIME, P3.F.DEPTH_AUGMENTED, '.-', 'Color',COLOR_ORDER(3,:),'HandleVisibility','off')
% plot(T0 + P4.F.TIME, P4.F.DEPTH_AUGMENTED, '.-', 'Color',COLOR_ORDER(4,:),'HandleVisibility','off')
plot(T0 + P5.F.TIME, P5.F.DEPTH_AUGMENTED, '.-', 'Color',COLOR_ORDER(5,:),'HandleVisibility','off')
plot(T0 + P6.F.TIME, P6.F.DEPTH_AUGMENTED, '.-', 'Color',COLOR_ORDER(6,:),'HandleVisibility','off')
plot(T0 + P7.F.TIME, P7.F.DEPTH_AUGMENTED, '.-', 'Color',COLOR_ORDER(7,:),'HandleVisibility','off')
plot(nan,nan, '.-', 'Color',COLOR_ORDER(1,:))
plot(nan,nan, '.-', 'Color',COLOR_ORDER(2,:))
plot(nan,nan, '.-', 'Color',COLOR_ORDER(3,:))
% plot(nan,nan, '.-', 'Color',COLOR_ORDER(4,:))
plot(nan,nan, '.-', 'Color',COLOR_ORDER(5,:))
plot(nan,nan, '.-', 'Color',COLOR_ORDER(6,:))
plot(nan,nan, '.-', 'Color',COLOR_ORDER(7,:))
legend('P1','P2','P3',     'P5','P6','P7')
datetick

%% Calculate steric height down to a fixed (not moving) depth:

for mi = 1:length(MM_list)
    MM = MM_list{mi};
    RHO_0 = eval([MM '.P.RHO_0;']); % Same for all moorings

    INTERP_FIXEDDEPTH_TO_PROFTIME = nan(size(eval([MM '.P.TIME_STERIC_HEIGHT']),1),...
                                        size(eval([MM '.F.DEPTH_AUGMENTED']),2));
    INTERP_FIXEDRHO_TO_PROFTIME = nan(size(eval([MM '.P.TIME_STERIC_HEIGHT']),1),...
                                      size(eval([MM '.F.DEPTH_AUGMENTED']),2));
    % INTERP_FIXEDRHO_TO_PROFTIME is [PROF_TIME by DEPTH_AUG] in size

    % Define time series as their own variables:
    warning(['NOTE: These solutions are integrated from the bottom-to-the-top, whereas the previous ' ...
             'steric height analysis was from the top-down, therefore to keep these congruent, we ' ...
             'subtract the corrections when we would have added them if integration were the other way:'])
    if strcmp(MM(1),'P')
        for ii = 1:size(INTERP_FIXEDRHO_TO_PROFTIME,2)
            INTERP_FIXEDDEPTH_TO_PROFTIME(:,ii) = interp1( ...
                eval([MM '.F.TIME']), ...
                eval([MM '.F.DEPTH_AUGMENTED(:,ii)']), ...
                eval([MM '.P.TIME_STERIC_HEIGHT']));
            INTERP_FIXEDRHO_TO_PROFTIME(:,ii) = interp1( ...
                eval([MM '.F.TIME(~sum(' MM '.F.QC_FLAG,2))']), ...
                eval([MM '.F.RHO_AUGMENTED(~sum(' MM '.F.QC_FLAG,2),ii)']), ...
                eval([MM '.P.TIME_STERIC_HEIGHT'])) - ... remove mean
                interp1( ...
                eval([MM '.F.RHO_TIMEMEAN_DEPTH']), ...
                eval([MM '.F.RHO_TIMEMEAN']), ...
                INTERP_FIXEDDEPTH_TO_PROFTIME(:,ii));
        end
        % Use the trend line to get the area of the trapezoid:
        eval([MM '_sh_600  = ' MM '.P.STERIC_HEIGHT_ANOMALY - ' ...
             '(0.5/RHO_0)*(600 - INTERP_FIXEDDEPTH_TO_PROFTIME(:,1)).*' ...
             '(2*INTERP_FIXEDRHO_TO_PROFTIME(:,1) + ' ...
             '[([INTERP_FIXEDRHO_TO_PROFTIME(:,2) - INTERP_FIXEDRHO_TO_PROFTIME(:,1)]./' ...
             '  [INTERP_FIXEDDEPTH_TO_PROFTIME(:,2) - INTERP_FIXEDDEPTH_TO_PROFTIME(:,1)]).*' ...
             '(600 - INTERP_FIXEDDEPTH_TO_PROFTIME(:,1))]);']);
        % eval([MM '_sh_600  = ' MM '_sh_600(' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM);']);
    elseif strcmp(MM(1),'S')
        QC_FLAG_AUG = eval(['[zeros(size(' MM '.F.TIME)) , ' MM '.F.QC_FLAG]']);
        for ii = 1:size(INTERP_FIXEDRHO_TO_PROFTIME,2)
            INTERP_FIXEDDEPTH_TO_PROFTIME(:,ii) = interp1( ...
                eval([MM '.F.TIME']), ...
                eval([MM '.F.DEPTH_AUGMENTED(:,ii)']), ...
                eval([MM '.P.TIME_STERIC_HEIGHT']));
            INTERP_FIXEDRHO_TO_PROFTIME(:,ii) = interp1( ...
                eval([MM '.F.TIME(~QC_FLAG_AUG(:,ii))']), ...
                eval([MM '.F.RHO_AUGMENTED(~QC_FLAG_AUG(:,ii),ii)']), ...
                eval([MM '.P.TIME_STERIC_HEIGHT'])) - ... remove mean
                interp1( ...
                eval([MM '.F.RHO_TIMEMEAN_DEPTH']), ...
                eval([MM '.F.RHO_TIMEMEAN']), ...
                INTERP_FIXEDDEPTH_TO_PROFTIME(:,ii));
        end
        eval([MM '_sh_600  = nan(size(' MM '.P.TIME_STERIC_HEIGHT));']);
        eval([MM '_sh_1800 = nan(size(' MM '.P.TIME_STERIC_HEIGHT));']);
        for ti = 1:length(eval([MM '.P.TIME_STERIC_HEIGHT']))
            % 600: (if-statement because sometimes the minimum observed
            % depth in the chain is below the desired level and we are
            % choosing not to use the line-fitting method here)
            if min(INTERP_FIXEDDEPTH_TO_PROFTIME(ti,isfinite(INTERP_FIXEDDEPTH_TO_PROFTIME(ti,:)))) < 600
                RHO_at_600 = interp1(INTERP_FIXEDDEPTH_TO_PROFTIME(ti,isfinite(INTERP_FIXEDDEPTH_TO_PROFTIME(ti,:))), ...
                                     INTERP_FIXEDRHO_TO_PROFTIME(  ti,isfinite(INTERP_FIXEDDEPTH_TO_PROFTIME(ti,:))),600);
                eval([MM '_sh_600(ti)  = ' MM '.P.STERIC_HEIGHT_ANOMALY(ti) + ' ...
                    '-trapz([INTERP_FIXEDDEPTH_TO_PROFTIME(ti,INTERP_FIXEDDEPTH_TO_PROFTIME(ti,:) < 600) 600]'',' ...
                    '       [INTERP_FIXEDRHO_TO_PROFTIME(  ti,INTERP_FIXEDDEPTH_TO_PROFTIME(ti,:) < 600) RHO_at_600]''/RHO_0);']);
            else
                eval([MM '_sh_600(ti) = NaN;']);
            end
            
            % 1800: (In our data, we always have data above and below 1800,
            % so do not do the same thing as for 600 here)
            RHO_at_1800 = interp1(INTERP_FIXEDDEPTH_TO_PROFTIME(ti,isfinite(INTERP_FIXEDDEPTH_TO_PROFTIME(ti,:))), ...
                                 INTERP_FIXEDRHO_TO_PROFTIME(  ti,isfinite(INTERP_FIXEDDEPTH_TO_PROFTIME(ti,:))),1800);
            eval([MM '_sh_1800(ti)  = ' MM '.P.STERIC_HEIGHT_ANOMALY(ti) + ' ...
                 '-trapz([INTERP_FIXEDDEPTH_TO_PROFTIME(ti,INTERP_FIXEDDEPTH_TO_PROFTIME(ti,:) < 1800) 1800]'',' ...
                 '       [INTERP_FIXEDRHO_TO_PROFTIME(  ti,INTERP_FIXEDDEPTH_TO_PROFTIME(ti,:) < 1800) RHO_at_1800]''/RHO_0);']);
        end
        % eval([MM '_sh_600  = ' MM '_sh_600( ' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM);']);
        % eval([MM '_sh_1800 = ' MM '_sh_1800(' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM);']);
    else
    end
end

%% Visualize the steric height of FIXED and PROF together:
close all

T0 = datenum('1950-01-01');

% Interval-averaging width (i.e. where SH is averaged in T(i) +- dT/2)
dT_interval_avg = [1/24]; % units of days

% Two options for interpolation: regular linear, or interval-average:
interp1nan = @(X,V,Xq) interp1(X(isfinite(X)&isfinite(V)),V(isfinite(X)&isfinite(V)),Xq);
% See auxiliary function(s) at the bottom, where MATLAB requires it to be.
interp1IA  = @(X,V,Xq) interp1IA_function(X,V,Xq,dT_interval_avg);
interp1_custom = interp1IA; % assign which one you want (interval average is better, the other is just for testing)

for mi = 1:length(MM_list)
    MM = MM_list{mi};
    figure('Color',[1 1 1])

    % % % UN-FLAGGED:
    % plot(eval([MM '.P.TIME_STERIC_HEIGHT']) , eval([       MM '.P.STERIC_HEIGHT_ANOMALY']) , '.-'); hold on
    % plot(eval([MM '.F.TIME'])               , eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED,2)']) , '.-')
    
    % % % FLAGGED BY DEPTH RANGE (for profilers) and QC_FLAG (for fixed):
    SH_PROF = eval([       MM '.P.STERIC_HEIGHT_ANOMALY(' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)']);
    plot(T0 + eval([       MM '.P.TIME_STERIC_HEIGHT('    MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)']) , ...
         SH_PROF, '.-'); hold on
    % plot(T0 + eval(replace(['#.P.TIME_STERIC_HEIGHT(   #.P.STERIC_HEIGHT_MAX_PROF_DEPTH - #.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)'],'#',MM)) , ...
    %           eval(replace(['#.P.STERIC_HEIGHT_ANOMALY(#.P.STERIC_HEIGHT_MAX_PROF_DEPTH - #.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)'],'#',MM)) , '.-'); hold on
    if strcmp(MM(1),'P')
        plot(T0 + eval([MM '.F.TIME']), eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , '.-')
        SH_FIXED_INTERP = interp1_custom( T0 + eval([MM '.F.TIME']), ...
                         eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
                         T0 + eval([MM '.P.TIME_STERIC_HEIGHT']));
        plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT']), SH_FIXED_INTERP, '.-')
        plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT']), eval([MM '_sh_600']) - eval([MM '.P.STERIC_HEIGHT_ANOMALY']), '.-k')
        
        
        legend([MM ' PROF'], [MM ' FIXED'], [MM ' FIXED interpolated to PROF times'], ...
               'Actually fixed 600')
    elseif strcmp(MM(1),'S')
        plot(T0 + eval([MM '.F.TIME']), eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , '.-')

        plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT']), ...
             interp1_custom( T0 + eval([MM '.F.TIME']), ...
                         eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
                         T0 + eval([MM '.P.TIME_STERIC_HEIGHT'])) , ...
             '.-')


        plot(T0 + eval([MM '.F.TIME']), eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2, ''omitnan'')']) , '.-')
        plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT']), ...
             interp1_custom( T0 + eval([MM '.F.TIME']), ...
                         eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2, ''omitnan'')']) , ...
                         T0 + eval([MM '.P.TIME_STERIC_HEIGHT'])) , ...
             '.-')
        plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT']), eval([MM '_sh_600']) -  eval([MM '.P.STERIC_HEIGHT_ANOMALY']), '.-k')
        plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT']), eval([MM '_sh_1800']) - eval([MM '.P.STERIC_HEIGHT_ANOMALY']), '.-','Color',[1 1 1]*0.5)

        legend([MM ' PROF'], [MM ' FIXED (uppermost layer, ~600-~500 m)'], [MM ' FIXED (uppermost layer, ~600-~500 m) interpolated to PROF times'] , [MM ' FIXED (full column)'], [MM ' FIXED interpolated to PROF times'], ...
               'Actually fixed 600','Actually fixed 1800')

    else
    end

    datetick
    title(MM); set(gca,'FontSize',16)
end

%% Calculate SFs using different parameters

% Recall:        SF = mean([Y(x) - Y(x + r)].^2)

% Each realization of [Y(x) - Y(x + r)].^2 should be assumed to
% approximately simultaneous. The time interval therefore matters.

% PROF and FIXED steric heights are on different time grids, so use the one
% with the lower time resultion (profilers) to interpolate the other
% (FIXED), thus saving computational resources.

%% Visualize the data that go into a structure function

close all

% % First, all data together to ensure that nothing is strange:
figure('Color',[1 1 1])
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT('       MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)']) , ...
              eval([MM '_sh_600('                     MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)']) , ...
              '.-'); hold on
end
datetick

figure('Color',[1 1 1])
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    % plot3(T0 + eval([MM '.P.TIME_STERIC_HEIGHT('       MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)']) , ...
    %            eval([MM '.P.LATITUDE*ones(sum('        MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM),1)']), ...
    %            eval([MM '.P.STERIC_HEIGHT_ANOMALY('    MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)']) , ...
    %            '.-'); hold on
    scatter3(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)']) , ...
                  eval([MM '.P.LATITUDE*ones(sum('  MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM),1)/100']), ...
                  eval([MM '_sh_600('               MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)']) , ...
                  5,eval([MM '_sh_600('             MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)']),'filled'); hold on
end
datetick
axis vis3d
grid on

%% SPIKES
% Obviously that spike in S2 is wrong and was somehow not filtered out.
% eliminate it and other possible outliers through a jump filter, wherein
% jumps of some minimum value indicate that the value should be set to NaN.

close all

figure('Color',[1 1 1])
% figure('Color',[1 1 1])

SPIKE_IND = {};
for mi = 1:length(MM_list)
    MM = MM_list{mi};

    figure(1)
    SPIKE_THRESHOLD = 0.04;
    % TT = eval([MM '.P.TIME_STERIC_HEIGHT(   ' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)']) + T0;
    % TS = eval([MM '.P.STERIC_HEIGHT_ANOMALY(' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM)']);
    TT = eval([MM '.P.TIME_STERIC_HEIGHT']) + T0;
    % TS = eval([MM '.P.STERIC_HEIGHT_ANOMALY']);
    TS = eval([MM '_sh_600']);
    SPIKE_IND{mi} = [abs([0;diff(TS)]) > SPIKE_THRESHOLD & sign([0;diff(TS)])./sign([0;0;diff(TS(1:[end-1]))]) == -1 ];
    TS(SPIKE_IND{mi}) = NaN;
    plot(TT,TS,'.-'); hold on

    figure(2)
    % plot(TT, diff([TS;NaN]), '.-'); hold on
    histogram(diff([TS;NaN]), [-0.08]:[0.001]:[0.08]); hold on
end
figure(1)
datetick

% % % For testing purposes, reset SPIKE_IND{...} to all false:
% for mi = 1:length(MM_list)
%     SPIKE_IND{mi} = 0*SPIKE_IND{mi};
% end


%% PLOT TIME SERIES - Visualize the data within a specified interval

close all

TEST_INTERVAL = [datenum('2023-02-15 00:00:00') datenum('2023-10-15 00:00:00')] - T0; % all times
% TEST_INTERVAL = [datenum('2023-06-01 00:00:00') datenum('2023-06-02 00:00:00')] - T0;
% TEST_INTERVAL = [datenum('2023-08-07 00:00:00') datenum('2023-08-08 00:00:00')] - T0;

MM_list_lats = [];
for mi=1:length(MM_list)
    MM = MM_list{mi};
    MM_list_lats(mi) = eval([MM '.P.LATITUDE;']);
end
[~,LatOrder] = sort(MM_list_lats);
LatOrder = flip(LatOrder);
COLOR = turbo(length(MM_list));

% First plot in time all the steric height:
figure('Color',[1 1 1])
II = 1;
LEG_TEXT = {};
for mi = LatOrder%1:length(MM_list)
    MM = MM_list{mi};
    % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT('    MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
    %           eval([MM '.P.STERIC_HEIGHT_ANOMALY(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
    %           '.:', 'MarkerSize',10); hold on
    if strcmp(MM(1),'P')
        % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
        %           interp1_custom( T0 + eval([MM '.F.TIME']), ...
        %                  eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
        %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})'])), ...
        %           '.:', 'MarkerSize',10);
        
        P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ' ... <-- Time interval
                   '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                      '~SPIKE_IND{mi}']); % <-- Spike removal
        plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(   P_IND)']) , ...
                  100*[eval([MM '.P.STERIC_HEIGHT_ANOMALY(P_IND)']) + ...
                  interp1_custom( T0 + eval([MM '.F.TIME']), ...
                         eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
                         T0 + eval([MM '.P.TIME_STERIC_HEIGHT(P_IND)'])) ], ...
                  '.-','Color',COLOR(II,:),'HandleVisibility','off'); hold on
        plot(nan,nan,'.','Color',COLOR(II,:),'MarkerSize',50)
    elseif strcmp(MM(1),'S')
        % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})']), ...
        %      interp1_custom( T0 + eval([MM '.F.TIME']), ...
        %                  eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
        %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})'])) , ...
        %      '.-')

        P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ' ... <-- Time interval
                    '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                    '~SPIKE_IND{mi}']); % <-- Spike removal
        plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(   P_IND)']), ...
                 100*[eval([MM '.P.STERIC_HEIGHT_ANOMALY(P_IND)']) + ...
                 interp1_custom( T0 + eval([MM '.F.TIME']), ...
                             eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
                             T0 + eval([MM '.P.TIME_STERIC_HEIGHT(P_IND)'])) ], ...
                 '.-','Color',COLOR(II,:),'HandleVisibility','off'); hold on
        plot(nan,nan,'.','Color',COLOR(II,:),'MarkerSize',50)
    else
    end
    LEG_TEXT = {LEG_TEXT{:}, MM};
    II = II + 1;
end
legend(LEG_TEXT,'Orientation','horizontal','Location','northwest')
datetick
% title(['Steric height data from ' datestr(T0 + TEST_INTERVAL(1)) ' to ' datestr(T0 + TEST_INTERVAL(2))])
ylabel(['Steric height anomaly' newline ' from 600 m (cm)'])
set(gca,'XLim',T0 + TEST_INTERVAL([1 end]))
set(gca,'FontSize',30)

set(gcf,'Position',[       -1685         324        1395         540])

disp(['24-Feb-2023 22:49:11 is the first time with 6 records at once'])
disp(['15-Sep-2023 14:51:09 is the last  time with 6 records at once'])

% %%
% exportgraphics(gcf,...
%     '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/STM_POSTER_TIMESERIES.pdf',...
%     'BackgroundColor','none','ContentType','vector')

%% Scatter them by separation and 2nd order SF:

close all

DATA_in_interval = [];
XY_in_interval = [];
figure('Color',[1 1 1])
TEST_INTERVAL = [datenum('2023-05-01 00:00:00') datenum('2023-05-01 02:00:00')] - T0;
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    % % % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT('    MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
    % % %           eval([MM '.P.STERIC_HEIGHT_ANOMALY(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
    % % %           '.:', 'MarkerSize',10); hold on
    if strcmp(MM(1),'P')
        % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
        %           interp1_custom( T0 + eval([MM '.F.TIME']), ...
        %                  eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
        %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})'])), ...
        %           '.:', 'MarkerSize',10);
        
        P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ' ... <-- Time interval
                   '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                   '~SPIKE_IND{mi}']); % <-- Spike removal
        DATA_in_interval = [DATA_in_interval ; ...
                  eval([MM '.P.STERIC_HEIGHT_ANOMALY(P_IND)']) + ...
                  interp1_custom( T0 + eval([MM '.F.TIME']), ...
                         eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
                         T0 + eval([MM '.P.TIME_STERIC_HEIGHT(P_IND)'])) ];
        % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
        %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
        XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)'])];
    elseif strcmp(MM(1),'S')
        % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})']), ...
        %      interp1_custom( T0 + eval([MM '.F.TIME']), ...
        %                  eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
        %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})'])) , ...
        %      '.-')

        P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ' ... <-- Time interval
                    '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                    '~SPIKE_IND{mi}']); % <-- Spike removal
        DATA_in_interval = [DATA_in_interval ; ...
                  eval([MM '.P.STERIC_HEIGHT_ANOMALY(P_IND)']) + ...
                 interp1_custom( T0 + eval([MM '.F.TIME']), ...
                             eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
                             T0 + eval([MM '.P.TIME_STERIC_HEIGHT(P_IND)'])) ];
        % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
        %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
        XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)'])];
    else
    end
end
SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;
% SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % naive cartesian way
DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
                              [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
% 
scatter(SEP_MAT_in_interval(:), SF_MAT_in_interval(:), 'filled'); hold on
    BIN_EDGES = [(-0.045):0.09:0.945]*DD;
plot([BIN_EDGES(1:[end-1]) + BIN_EDGES(2:end)]/2,...
     interval_avg(SEP_MAT_in_interval(:), SF_MAT_in_interval(:), BIN_EDGES), '.-','MarkerSize',10)
xlabel('|x_1 - x_2|')
ylabel('[\eta(x_1) - \eta(x_2)].^2')



% figure
% subplot(121)
% imagesc(SEP_MAT_in_interval); colorbar
% title('|x_1 - x_2|')
% subplot(122)
% imagesc(SF_MAT_in_interval); colorbar
% title('mean([\eta(x_1) - \eta(x_2)].^2)')



%%
%%
%%
%%
%%
%% SP.600m.1hr - Go through each time interval and calculate the squared differences, store together, and then fit

close all

SF_TIME_STEP = 1/24; % days
% TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];

DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_VEC = [];
SF_VEC = [];

tic
for ti = 1:[length(TIME_BIN_EDGES) - 1]
    DATA_in_interval = [];
    XY_in_interval = [];
    SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
    for mi = 1:length(MM_list)
        MM = MM_list{mi};
        if strcmp(MM(1),'P')
            % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
            %           interp1_custom( T0 + eval([MM '.F.TIME']), ...
            %                  eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
            %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])), ...
            %           '.:', 'MarkerSize',10);
            
            P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                       '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                       '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                 mean(eval([MM '_sh_600(P_IND)']) , 'omitnan') ];
            % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
            %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
             mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        elseif strcmp(MM(1),'S')
            % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']), ...
            %      interp1_custom( T0 + eval([MM '.F.TIME']), ...
            %                  eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
            %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])) , ...
            %      '.-')
    
            P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                        '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                        '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                mean(eval([MM '_sh_600(P_IND)']) , 'omitnan') ];
            % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
            %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        else
        end
    end
    % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
                              [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
    end
end
toc

% To reiterate, the results from this section are for S an P moorings and
% combining profilers with the upper-most fixed CTD (~600 m).
Separation_S_P_600m_1hr = SEP_VEC;
StructFunc_S_P_600m_1hr = SF_VEC;

%% SP.600m.24hr
% allps=interval_avg(Separation_S_P_600m_24hr,StructFunc_S_P_600m_24hr,BIN_EDGES);
% % figure;plot(BIN_CENTR,nop1,'.-'); hold on;
% % plot(BIN_CENTR,nop2,'.-');plot(BIN_CENTR,nop3,'.-');plot(BIN_CENTR,nop5,'.-')
% % plot(BIN_CENTR,nop6,'.-');plot(BIN_CENTR,nop7,'.-');plot(BIN_CENTR,allps,'k.--');

close all

SF_TIME_STEP = 24/24; % days
TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];

DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_VEC = [];
SF_VEC = [];

tic
for ti = 1:[length(TIME_BIN_EDGES) - 1]
    DATA_in_interval = [];
    XY_in_interval = [];
    SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
    for mi = 1:length(MM_list)
        MM = MM_list{mi};
        if strcmp(MM(1),'P')

            P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                       '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                       '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        elseif strcmp(MM(1),'S')
            P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                        '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                        '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        else
        end
    end
    % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
                              [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
    end
end
toc

Separation_S_P_600m_24hr = SEP_VEC;
StructFunc_S_P_600m_24hr = SF_VEC;

clear SF_VEC SEP_VEC

%% SP.600m.1hr.NoP1

close all

SF_TIME_STEP = 1/24; % days
% TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];

DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_VEC = [];
SF_VEC = [];

tic
for ti = 1:[length(TIME_BIN_EDGES) - 1]
    DATA_in_interval = [];
    XY_in_interval = [];
    SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
    for mi = 1:length(MM_list)
        MM = MM_list{mi};
        if strcmp(MM(1),'P') & ~[strcmp(MM,'P1')]
            % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
            %           interp1_custom( T0 + eval([MM '.F.TIME']), ...
            %                  eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
            %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])), ...
            %           '.:', 'MarkerSize',10);
            
            P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                       '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                       '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
            %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
             mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        elseif strcmp(MM(1),'S')
            % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']), ...
            %      interp1_custom( T0 + eval([MM '.F.TIME']), ...
            %                  eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
            %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])) , ...
            %      '.-')
    
            P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                        '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                        '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
            %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        else
        end
    end
    % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
                              [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
    end
end
toc

% To reiterate, the results from this section are for S an P moorings and
% combining profilers with the upper-most fixed CTD (~600 m).
Separation_S_P_600m_1hr_noP1 = SEP_VEC;
StructFunc_S_P_600m_1hr_noP1 = SF_VEC;

%% SP.600m.24hr.NoP1

close all

SF_TIME_STEP = 24/24; % days
% TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];

DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_VEC = [];
SF_VEC = [];

tic
for ti = 1:[length(TIME_BIN_EDGES) - 1]
    DATA_in_interval = [];
    XY_in_interval = [];
    SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
    for mi = 1:length(MM_list)
        MM = MM_list{mi};
        if strcmp(MM(1),'P') & ~[strcmp(MM,'P1')]
            % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
            %           interp1_custom( T0 + eval([MM '.F.TIME']), ...
            %                  eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
            %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])), ...
            %           '.:', 'MarkerSize',10);
            
            P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                       '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                       '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
            %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
             mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        elseif strcmp(MM(1),'S')
            % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']), ...
            %      interp1_custom( T0 + eval([MM '.F.TIME']), ...
            %                  eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
            %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])) , ...
            %      '.-')
    
            P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                        '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                        '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
            %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        else
        end
    end
    % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
                              [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
    end
end
toc

% To reiterate, the results from this section are for S an P moorings and
% combining profilers with the upper-most fixed CTD (~600 m).
Separation_S_P_600m_24hr_noP1 = SEP_VEC;
StructFunc_S_P_600m_24hr_noP1 = SF_VEC;

%% %%%% S.600m.1hr
% allps=interval_avg(Separation_S_P_600m_24hr,StructFunc_S_P_600m_24hr,BIN_EDGES);
% % figure;plot(BIN_CENTR,nop1,'.-'); hold on;
% % plot(BIN_CENTR,nop2,'.-');plot(BIN_CENTR,nop3,'.-');plot(BIN_CENTR,nop5,'.-')
% % plot(BIN_CENTR,nop6,'.-');plot(BIN_CENTR,nop7,'.-');plot(BIN_CENTR,allps,'k.--');

close all

SF_TIME_STEP = 1/24; % days
TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];

DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_VEC = [];
SF_VEC = [];

TIME_BIN_CENTERS = [TIME_BIN_EDGES(2:end) + TIME_BIN_EDGES(1:[end-1])]/2;

tic
for swot_ti = 1:length(SWOT.time)
    if isfinite(median(SWOT.time{swot_ti},'omitnan'))
        ti = dsearchn(TIME_BIN_CENTERS',median(SWOT.time{swot_ti},'omitnan'));
        DATA_in_interval = [];
        XY_in_interval = [];
        SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
        for mi = 1:length(MM_list)
            MM = MM_list{mi};
            if strcmp(MM(1),'P') & ~[strcmp(MM,'P1')]
            elseif strcmp(MM(1),'S')
                P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                    '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                    '~SPIKE_IND{mi}']); % <-- Spike removal
                DATA_in_interval = [DATA_in_interval ; ...
                    mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
                XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                    mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
            else
            end
        end
        % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
        SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
            [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
        SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

        SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
        SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

        % if ~mod(ti,10)
        % disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
        % end
        disp([num2str(ti) '   /   ' datestr(TIME_BIN_EDGES(ti)) '   /   ' datestr(median(SWOT.time{swot_ti},'omitnan'))])
    else
    end
end
toc

Separation_S_600m_1hr = SEP_VEC;
StructFunc_S_600m_1hr = SF_VEC;

clear SF_VEC SEP_VEC


%% SP.600m.24hr
% allps=interval_avg(Separation_S_P_600m_24hr,StructFunc_S_P_600m_24hr,BIN_EDGES);
% % figure;plot(BIN_CENTR,nop1,'.-'); hold on;
% % plot(BIN_CENTR,nop2,'.-');plot(BIN_CENTR,nop3,'.-');plot(BIN_CENTR,nop5,'.-')
% % plot(BIN_CENTR,nop6,'.-');plot(BIN_CENTR,nop7,'.-');plot(BIN_CENTR,allps,'k.--');

close all

SF_TIME_STEP = 24/24; % days
TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];

DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_VEC = [];
SF_VEC = [];

tic
for ti = 1:[length(TIME_BIN_EDGES) - 1]
    DATA_in_interval = [];
    XY_in_interval = [];
    SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
    for mi = 1:length(MM_list)
        MM = MM_list{mi};
        if strcmp(MM(1),'P')

            P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                       '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                       '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        elseif strcmp(MM(1),'S')
            P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                        '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                        '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        else
        end
    end
    % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
                              [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
    end
end
toc

Separation_S_P_600m_24hr = SEP_VEC;
StructFunc_S_P_600m_24hr = SF_VEC;

clear SF_VEC SEP_VEC



%% S.1800m.1hr

% S1.F.DEPTH_NOMINAL'
% S2.F.DEPTH_NOMINAL'
% S3.F.DEPTH_NOMINAL'
% S4.F.DEPTH_NOMINAL'
%          605         804        1205        1814        2523        3432        4427
%          604         803        1204        1706        2361        3268        4373
%          605         804        1205        1814        2523        3432        4428
%          608         798        1199        1813        2522        3432        4327
i_deepest_F = 4;

close all

SF_TIME_STEP = 1/24; % days
TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];

DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_VEC = [];
SF_VEC = [];

% S_MM_list = {'S1','S2','S3','S4'};

tic
for ti = 1:[length(TIME_BIN_EDGES) - 1]
    DATA_in_interval = [];
    XY_in_interval = [];
    SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
    for mi = 1:length(MM_list)
        MM = MM_list{mi};
        if strcmp(MM(1),'P')
            % % % P-MOORINGS OMITTED FROM THIS ANALYSIS
        elseif strcmp(MM(1),'S')
            P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                        '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                        '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                mean(eval([MM '_sh_1800(P_IND)']), 'omitnan') ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        else
        end
    end
    % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
                              [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
    end
end
toc

Separation_S_deep_1hr = SEP_VEC;
StructFunc_S_deep_1hr = SF_VEC;

%% S.1800m.24hr


close all

SF_TIME_STEP = 24/24; % days
TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];

DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_VEC = [];
SF_VEC = [];

% S_MM_list = {'S1','S2','S3','S4'};

tic
for ti = 1:[length(TIME_BIN_EDGES) - 1]
    DATA_in_interval = [];
    XY_in_interval = [];
    SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
    for mi = 1:length(MM_list)
        MM = MM_list{mi};
        if strcmp(MM(1),'P')
            % % % P-MOORINGS OMITTED FROM THIS ANALYSIS
        elseif strcmp(MM(1),'S')
            P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                        '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                        '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                mean(eval([MM '_sh_1800(P_IND)']), 'omitnan') ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        else
        end
    end
    % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
                              [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
    end
end
toc

Separation_S_deep_24hr = SEP_VEC;
StructFunc_S_deep_24hr = SF_VEC;


%% SWOT_at_S_P_G

% Do the same kind of analysis for SWOT data at S and P mooring locations

% SWOT:
    %                  time: {1×180 cell}
    %          ssha_karin_2: {1×180 cell}
    %     ssha_karin_2_qual: {1×180 cell}
    %                  name: {1×180 cell}
    %                   lat: {1×180 cell}
    %                   lon: {1×180 cell}
    %    internal_tide_hret: {1×180 cell}
    %      height_cor_xover: {1×180 cell}
    % height_cor_xover_qual: {1×180 cell}

% Set vector of nominal locations of moorings:
M_lonlat_nominal = nan(length(MM_list),2);
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    M_lonlat_nominal(mi,:) = [eval([MM '.P.LONGITUDE']) , eval([MM '.P.LATITUDE'])];
end

MM_GG_list = {'P1','P2','P3',    'P5','P6','P7',    'S1','S2','S3','S4',     'ru32_600','ru38_600'};
% MM_GG_list = {'P1','P2','P3',    'P5','P6','P7',    'S1','S2','S3','S4'};

SEP_VEC = [];
SF_VEC = [];
tic
for ti = 1:length(SWOT.time)

    SWOT_M_time = [];
    SWOT_M_ssha = [];
    SWOT_M_ssha_qc = [];
    SWOT_M_rollerror = [];
    SWOT_M_rollerror_qc = [];
    SWOT_M_hret = [];
    SWOT_M_ll = [];
    % for mi = 1:length(MM_list)
        % MM = MM_list{mi};
    for mi = 1:length(MM_GG_list)
        MM = MM_GG_list{mi};
        if strcmp(MM(1),'S') || strcmp(MM(1),'P') % at mooring location
            MM_GPS_LON = eval([[MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY']]);
            MM_GPS_LAT = eval([[MM '.GPS.LATITUDE_GPS_SURFACE_BUOY']]);
            Mooring_Loc_Mat = abs(SWOT.lon{ti} + 1i*SWOT.lat{ti} - [...
                MM_GPS_LON(dsearchn(eval([MM '.GPS.TIME_GPS_SURFACE_BUOY']) + T0, median(SWOT.time{ti}(:)))) + ...
                MM_GPS_LAT(dsearchn(eval([MM '.GPS.TIME_GPS_SURFACE_BUOY']) + T0, median(SWOT.time{ti}(:))))*1i] );
            % ^ the use of median(SWOT.time{ti}(:)) will mean that the true
            % time will be off by a maximum of 3-6 minutes, which should
            % be ok considering the time steps of profilers is 0.5-4 hr
        elseif strcmp(MM(1),'r') % at Rutgers glider location
            GG_GPS_LON = eval([[MM '.LONGITUDE_PROFILE']]);
            GG_GPS_LAT = eval([[MM '.LATITUDE_PROFILE']]);
            Mooring_Loc_Mat = abs(SWOT.lon{ti} + 1i*SWOT.lat{ti} - [...
                GG_GPS_LON(dsearchn(eval([MM '.TIME_STERIC_HEIGHT']) + T0, median(SWOT.time{ti}(:)))) + ...
                GG_GPS_LAT(dsearchn(eval([MM '.TIME_STERIC_HEIGHT']) + T0, median(SWOT.time{ti}(:))))*1i] );
        else
        end
        SWOT_M_time_ti = SWOT.time{ti}( ~~sum(Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)),1)' );
        SWOT_M_time = [SWOT_M_time ; SWOT_M_time_ti];

        SWOT_M_ssha_ti = SWOT.ssha_karin_2{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_ssha = [SWOT_M_ssha ; SWOT_M_ssha_ti];

        SWOT_M_ssha_qc_ti = SWOT.ssha_karin_2_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_ssha_qc = [SWOT_M_ssha_qc ; SWOT_M_ssha_qc_ti];

        SWOT_M_rollerror_ti = SWOT.height_cor_xover{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_rollerror = [SWOT_M_rollerror ; SWOT_M_rollerror_ti];

        SWOT_M_rollerror_qc_ti = SWOT.height_cor_xover_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_rollerror_qc = [SWOT_M_rollerror_qc ; SWOT_M_rollerror_qc_ti];

        SWOT_M_hret_ti = SWOT.internal_tide_hret{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_hret = [SWOT_M_hret ; SWOT_M_hret_ti];

        SWOT_M_ll_ti = SWOT.lon{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) ) + ...
                       SWOT.lat{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) )*1i;
        SWOT_M_ll = [SWOT_M_ll; SWOT_M_ll_ti];
    end

    DATA_in_interval = [SWOT_M_ssha + SWOT_M_hret + SWOT_M_rollerror].* ...
                      ~[SWOT_M_ssha_qc + SWOT_M_rollerror_qc]./~[SWOT_M_ssha_qc + SWOT_M_rollerror_qc];
    XY_in_interval = SWOT_M_ll;



    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
        [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str( length(SWOT.time) )])
    end
end
toc
warning('This calculation adds internal_tide_hret to ssha_karin_2 and applies flags.')
warning('SWOT xover passed are every 11 and 13 hours, but the variable is labeled as "12hr".')

Separation_SWOT_moor_12hr = SEP_VEC;
StructFunc_SWOT_moor_12hr = SF_VEC;

% % % OLD user manual:
% ssha_karin_2 = ssh_karin_2 – mean_sea_surface_sol1 – solid_earth_tide – ...
%                ocean_tide_sol1 – pole_tide - dac 

% % % PO.DAAC: <https://podaac.github.io/tutorials/notebooks/datasets/Localmachine_SWOT_Oceanography.html>
% Sea surface height anomaly from the KaRIn measurement = ssh_karin_2 -
% mean_sea_surface_cnescls - solid_earth_tide - ocean_tide_fes –
% internal_tide_hret - pole_tide - dac.

%% SWOT_at_S

% Do the same kind of analysis for SWOT data at S and P mooring locations

% SWOT:
    %                  time: {1×180 cell}
    %          ssha_karin_2: {1×180 cell}
    %     ssha_karin_2_qual: {1×180 cell}
    %                  name: {1×180 cell}
    %                   lat: {1×180 cell}
    %                   lon: {1×180 cell}
    %    internal_tide_hret: {1×180 cell}
    %      height_cor_xover: {1×180 cell}
    % height_cor_xover_qual: {1×180 cell}

% Set vector of nominal locations of moorings:
M_lonlat_nominal = nan(length(MM_list),2);
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    M_lonlat_nominal(mi,:) = [eval([MM '.P.LONGITUDE']) , eval([MM '.P.LATITUDE'])];
end

MM_GG_list = {'S1','S2','S3','S4'};
% MM_GG_list = {'P1','P2','P3',    'P5','P6','P7',    'S1','S2','S3','S4'};

SEP_VEC = [];
SF_VEC = [];
tic
for ti = 1:length(SWOT.time)

    SWOT_M_time = [];
    SWOT_M_ssha = [];
    SWOT_M_ssha_qc = [];
    SWOT_M_rollerror = [];
    SWOT_M_rollerror_qc = [];
    SWOT_M_hret = [];
    SWOT_M_ll = [];
    % for mi = 1:length(MM_list)
        % MM = MM_list{mi};
    for mi = 1:length(MM_GG_list)
        MM = MM_GG_list{mi};
        if strcmp(MM(1),'S') || strcmp(MM(1),'P') % at mooring location
            MM_GPS_LON = eval([[MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY']]);
            MM_GPS_LAT = eval([[MM '.GPS.LATITUDE_GPS_SURFACE_BUOY']]);
            Mooring_Loc_Mat = abs(SWOT.lon{ti} + 1i*SWOT.lat{ti} - [...
                MM_GPS_LON(dsearchn(eval([MM '.GPS.TIME_GPS_SURFACE_BUOY']) + T0, median(SWOT.time{ti}(:)))) + ...
                MM_GPS_LAT(dsearchn(eval([MM '.GPS.TIME_GPS_SURFACE_BUOY']) + T0, median(SWOT.time{ti}(:))))*1i] );
            % ^ the use of median(SWOT.time{ti}(:)) will mean that the true
            % time will be off by a maximum of 3-6 minutes, which should
            % be ok considering the time steps of profilers is 0.5-4 hr
        elseif strcmp(MM(1),'r') % at Rutgers glider location
            GG_GPS_LON = eval([[MM '.LONGITUDE_PROFILE']]);
            GG_GPS_LAT = eval([[MM '.LATITUDE_PROFILE']]);
            Mooring_Loc_Mat = abs(SWOT.lon{ti} + 1i*SWOT.lat{ti} - [...
                GG_GPS_LON(dsearchn(eval([MM '.TIME_STERIC_HEIGHT']) + T0, median(SWOT.time{ti}(:)))) + ...
                GG_GPS_LAT(dsearchn(eval([MM '.TIME_STERIC_HEIGHT']) + T0, median(SWOT.time{ti}(:))))*1i] );
        else
        end
        SWOT_M_time_ti = SWOT.time{ti}( ~~sum(Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)),1)' );
        SWOT_M_time = [SWOT_M_time ; SWOT_M_time_ti];

        SWOT_M_ssha_ti = SWOT.ssha_karin_2{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_ssha = [SWOT_M_ssha ; SWOT_M_ssha_ti];

        SWOT_M_ssha_qc_ti = SWOT.ssha_karin_2_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_ssha_qc = [SWOT_M_ssha_qc ; SWOT_M_ssha_qc_ti];

        SWOT_M_rollerror_ti = SWOT.height_cor_xover{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_rollerror = [SWOT_M_rollerror ; SWOT_M_rollerror_ti];

        SWOT_M_rollerror_qc_ti = SWOT.height_cor_xover_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_rollerror_qc = [SWOT_M_rollerror_qc ; SWOT_M_rollerror_qc_ti];

        SWOT_M_hret_ti = SWOT.internal_tide_hret{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_hret = [SWOT_M_hret ; SWOT_M_hret_ti];

        SWOT_M_ll_ti = SWOT.lon{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) ) + ...
                       SWOT.lat{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) )*1i;
        SWOT_M_ll = [SWOT_M_ll; SWOT_M_ll_ti];
    end

    DATA_in_interval = [SWOT_M_ssha + SWOT_M_hret + SWOT_M_rollerror].* ...
                      ~[SWOT_M_ssha_qc + SWOT_M_rollerror_qc]./~[SWOT_M_ssha_qc + SWOT_M_rollerror_qc];
    XY_in_interval = SWOT_M_ll;



    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
        [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str( length(SWOT.time) )])
    end
end
toc
warning('This calculation adds internal_tide_hret to ssha_karin_2 and applies flags.')
warning('SWOT xover passed are every 11 and 13 hours, but the variable is labeled as "12hr".')

Separation_SWOTatS_moor_12hr = SEP_VEC;
StructFunc_SWOTatS_moor_12hr = SF_VEC;

%% SWOT_at_S_P_G and other crossover diamonds

% Set vector of nominal locations of moorings:
MM_list_ordered = {'S1','P1','P2','S2','P3','S3','P5','P6','S4','P7'};
M_lonlat_nominal = nan(length(MM_list_ordered),2);
for mi = 1:length(MM_list_ordered)
    MM = MM_list_ordered{mi};
    M_lonlat_nominal(mi,:) = [eval([MM '.P.LONGITUDE']) , eval([MM '.P.LATITUDE'])];
end
dLON = median(abs(diff(M_lonlat_nominal(:,1))));
dLAT = median(abs(diff(M_lonlat_nominal(:,2))));
DIMAMOND_OFFSET = [[-4*dLON 0 4*dLON -4*dLON  4*dLON 0      0       8*dLON -8*dLON]' , ...
                   [-4*dLAT 0 4*dLAT  4*dLAT -4*dLAT 8*dLAT -8*dLAT 0      0      ]'];
M_lonlat_nominal = [M_lonlat_nominal ; ... E diamond (original)
    [[-125.4000 37.1651] + DIMAMOND_OFFSET] ; ... N diamond, hand-picked center
    [[-125.3915 34.1991] + DIMAMOND_OFFSET] ; ... S diamond, hand-picked center
    [[-125.7818 35.6967] + DIMAMOND_OFFSET] ; ... W diamond, hand-picked center
    ];
% % figure;plot(M_lonlat_nominal(:,1),M_lonlat_nominal(:,2),'-*');axis equal


% % Alternative: one straight line that only looks along one track (the
% % longer one):
% ASCENDING_DESCENDING = nan(length(SWOT),1);
% for jj = 1:length(SWOT)
%     ASCENDING_DESCENDING(jj) = sign(diff((SWOT.lat{jj}(50,1:2))));
% end

SEP_VEC = [];
SF_VEC = [];
tic
for ti = 1:length(SWOT.time)

    SWOT_M_time = [];
    SWOT_M_ssha = [];
    SWOT_M_ssha_qc = [];
    SWOT_M_rollerror = [];
    SWOT_M_rollerror_qc = [];
    SWOT_M_hret = [];
    SWOT_M_ll = [];
    % for mi = 1:length(MM_list)
        % MM = MM_list{mi};
    for mi = 1:size(M_lonlat_nominal,1)
            % At nominal locations, which is easier for comparing to fake
            % moorings in other xover diamonds

            MM_GPS_LON = M_lonlat_nominal(mi,1);
            MM_GPS_LAT = M_lonlat_nominal(mi,2);
            Mooring_Loc_Mat = abs(SWOT.lon{ti} + 1i*SWOT.lat{ti} - [MM_GPS_LON + MM_GPS_LAT*1i] );

        SWOT_M_time_ti = SWOT.time{ti}( ~~sum(Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)),1)' );
        SWOT_M_time = [SWOT_M_time ; SWOT_M_time_ti];

        SWOT_M_ssha_ti = SWOT.ssha_karin_2{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_ssha = [SWOT_M_ssha ; SWOT_M_ssha_ti];

        SWOT_M_ssha_qc_ti = SWOT.ssha_karin_2_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_ssha_qc = [SWOT_M_ssha_qc ; SWOT_M_ssha_qc_ti];

        SWOT_M_rollerror_ti = SWOT.height_cor_xover{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_rollerror = [SWOT_M_rollerror ; SWOT_M_rollerror_ti];

        SWOT_M_rollerror_qc_ti = SWOT.height_cor_xover_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_rollerror_qc = [SWOT_M_rollerror_qc ; SWOT_M_rollerror_qc_ti];

        SWOT_M_hret_ti = SWOT.internal_tide_hret{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_hret = [SWOT_M_hret ; SWOT_M_hret_ti];

        SWOT_M_ll_ti = SWOT.lon{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) ) + ...
                       SWOT.lat{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) )*1i;
        SWOT_M_ll = [SWOT_M_ll; SWOT_M_ll_ti];
    end

    DATA_in_interval = [SWOT_M_ssha + SWOT_M_hret + SWOT_M_rollerror].* ...
                      ~[SWOT_M_ssha_qc + SWOT_M_rollerror_qc]./~[SWOT_M_ssha_qc + SWOT_M_rollerror_qc];
    XY_in_interval = SWOT_M_ll;



    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
        [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str( length(SWOT.time) )])
    end
end
toc
warning('This calculation adds internal_tide_hret to ssha_karin_2 and applies flags.')
warning('SWOT xover passed are every 11 and 13 hours, but the variable is labeled as "12hr".')

Separation_SWOT_extended_12hr = SEP_VEC;
StructFunc_SWOT_extended_12hr = SF_VEC;


%% SPG.600m.1hr - Moorings AND Gliders

close all

% Glider-mooring offset (add this TO GLIDERS). This was examined and estimated in <MooringGlider_Offset.m>:
GM_offset = 0; 1.2083/100;

SF_TIME_STEP = 1/24; % days
% TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];

DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_VEC = [];
SF_VEC = [];

MM_GG_list = {'P1','P2','P3',    'P5','P6','P7',    'S1','S2','S3','S4',     'ru32_600','ru38_600'};
% MM_GG_list = {'ru32_600','ru38_600'}; % for testing purposes

tic
for ti = 1:[length(TIME_BIN_EDGES) - 1]
    DATA_in_interval = [];
    XY_in_interval = [];
    SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
    for mi = 1:length(MM_GG_list)
        MM = MM_GG_list{mi};
        if strcmp(MM(1),'P') %& ~[strcmp(MM,'P1')]
            % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
            %           interp1_custom( T0 + eval([MM '.F.TIME']), ...
            %                  eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
            %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])), ...
            %           '.:', 'MarkerSize',10);
            
            P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                       '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                       '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
            %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        elseif strcmp(MM(1),'S')
            % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']), ...
            %      interp1_custom( T0 + eval([MM '.F.TIME']), ...
            %                  eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
            %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])) , ...
            %      '.-')
    
            P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                        '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                        '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
            %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        elseif strcmp(MM(1),'r') % Rutgers gliders
            P_IND = eval([  MM '.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                        '[' MM '.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.STERIC_HEIGHT_MIN_PROF_DEPTH] > G_RANGE_LIM' ... <-- Glider depth requirement
                        ]); % <-- No spike removal for gliders
            DATA_in_interval = [DATA_in_interval ; ...
                      mean(eval([MM '.STERIC_HEIGHT_ANOMALY(P_IND) + GM_offset']), 'omitnan') ]; % no distinction between P and F for gliders
            % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
            %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.LONGITUDE_PROFILE(P_IND) + 1i*' MM '.LATITUDE_PROFILE(P_IND)']), 'omitnan') ];
        else
        end
    end
    % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
                              [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
    end
end
toc

Separation_S_P_G_600m_1hr = SEP_VEC;
StructFunc_S_P_G_600m_1hr = SF_VEC;

%% SPG.600m.24hr - Moorings AND Gliders

close all

% GM_offset = 0;

SF_TIME_STEP = 24/24; % days
TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];

DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_VEC = [];
SF_VEC = [];

MM_GG_list = {'P1','P2','P3',    'P5','P6','P7',    'S1','S2','S3','S4',     'ru32_600','ru38_600'};
% MM_GG_list = {'ru32_600','ru38_600'}; % for testing purposes

tic
for ti = 1:[length(TIME_BIN_EDGES) - 1]
    DATA_in_interval = [];
    XY_in_interval = [];
    SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
    for mi = 1:length(MM_GG_list)
        MM = MM_GG_list{mi};
        if strcmp(MM(1),'P') %& ~[strcmp(MM,'P1')]

            P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                       '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                       '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        elseif strcmp(MM(1),'S')

            P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                        '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                        '~SPIKE_IND{mi}']); % <-- Spike removal
            DATA_in_interval = [DATA_in_interval ; ...
                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
        elseif strcmp(MM(1),'r') % Rutgers gliders
            P_IND = eval([  MM '.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                        '[' MM '.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.STERIC_HEIGHT_MIN_PROF_DEPTH] > G_RANGE_LIM' ... <-- Glider depth requirement
                        ]); % <-- No spike removal for gliders
            DATA_in_interval = [DATA_in_interval ; ...
                      mean(eval([MM '.STERIC_HEIGHT_ANOMALY(P_IND) + GM_offset']), 'omitnan') ]; % no distinction between P and F for gliders
            XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                  mean(eval([MM '.LONGITUDE_PROFILE(P_IND) + 1i*' MM '.LATITUDE_PROFILE(P_IND)']), 'omitnan') ];
        else
        end
    end
    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
                              [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
    end
end
toc

Separation_S_P_G_600m_24hr = SEP_VEC;
StructFunc_S_P_G_600m_24hr = SF_VEC;

%% Save workspace so you don't have to redo all these calculations
% FOO(:,3) = interval_avg(Separation_S_P_600m_1hr,StructFunc_S_P_600m_1hr,BIN_EDGES);

error(['Do not run this section if you have already saved the large file "SF_plot_fixed_bottom_workspace.mat"'])

README = ...
    ['This is the result of running all the processing in SF_plot_fixed_bottom.m ' ...
     'and WavenumberSpectrum_SWOT.m, for the onvenience of not having to re-run ' ...
     'everytime a figure needs to be slightly adjusted'];
save(['/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/' ...
      'SF_plot_fixed_bottom_workspace.mat'])
disp(['Done saving workspace.'])

%% Load the workspace if needed

% error('Only load if this is needed.')

INPUT = input(['Do you want to load ~1.3 GB into memory? Enter any number for "yes" or push\n' ...
    '"Enter" with a blank input for "No".\n']);

if isempty(INPUT)
    disp(['No data loaded'])
else
    load(['/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/' ...
        'SF_plot_fixed_bottom_workspace.mat'])
    disp('Loaded data')
end

%% Plot the multiple SF's together:

close all

figure('Color',[1 1 1])
dBIN = 0.1;

% % % SP.600m.1hr
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_EDGES = interval_avg(Separation_S_P_600m_1hr,Separation_S_P_600m_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_EDGES(1) - [BIN_EDGES(2) - BIN_EDGES(1)]/2 ], ...
               [BIN_EDGES(2:end) + BIN_EDGES(1:[end-1])]/2, ...
               [BIN_EDGES(end) + [BIN_EDGES(end) - BIN_EDGES(end-1)]/2] ];
% plot(Separation_S_P_600m_1hr, StructFunc_S_P_600m_1hr, '.', 'MarkerSize',10); hold on
plot([BIN_EDGES(1:[end-1]) + BIN_EDGES(2:end)]/2,...
     interval_avg(Separation_S_P_600m_1hr,StructFunc_S_P_600m_1hr,BIN_EDGES), '.-', 'MarkerSize',10); hold on

% % % SP.600m.24hr
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_EDGES = interval_avg(Separation_S_P_600m_24hr,Separation_S_P_600m_24hr,BIN_EDGES);
BIN_EDGES = [[ BIN_EDGES(1) - [BIN_EDGES(2) - BIN_EDGES(1)]/2 ], ...
               [BIN_EDGES(2:end) + BIN_EDGES(1:[end-1])]/2, ...
               [BIN_EDGES(end) + [BIN_EDGES(end) - BIN_EDGES(end-1)]/2] ];
% plot(Separation_S_P_600m_24hr, StructFunc_S_P_600m_24hr, '.', 'MarkerSize',10); hold on
plot([BIN_EDGES(1:[end-1]) + BIN_EDGES(2:end)]/2,...
     interval_avg(Separation_S_P_600m_24hr,StructFunc_S_P_600m_24hr,BIN_EDGES), '.-', 'MarkerSize',10); hold on

% % % S.4400m.1hr
BIN_EDGES = [(-3*dBIN/2):[3*dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_EDGES = interval_avg(Separation_S_deep_1hr,Separation_S_deep_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_EDGES(1) - [BIN_EDGES(2) - BIN_EDGES(1)]/2 ], ...
               [BIN_EDGES(2:end) + BIN_EDGES(1:[end-1])]/2, ...
               [BIN_EDGES(end) + [BIN_EDGES(end) - BIN_EDGES(end-1)]/2] ];
% plot(Separation_S_deep_1hr, StructFunc_S_deep_1hr, '.', 'MarkerSize',10); hold on
plot([BIN_EDGES(1:[end-1]) + BIN_EDGES(2:end)]/2,...
     interval_avg(Separation_S_deep_1hr,StructFunc_S_deep_1hr,BIN_EDGES), '.-', 'MarkerSize',10); hold on

% % % S.4400m.24hr
BIN_EDGES = [(-3*dBIN/2):[3*dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_EDGES = interval_avg(Separation_S_deep_24hr,Separation_S_deep_24hr,BIN_EDGES);
BIN_EDGES = [[ BIN_EDGES(1) - [BIN_EDGES(2) - BIN_EDGES(1)]/2 ], ...
               [BIN_EDGES(2:end) + BIN_EDGES(1:[end-1])]/2, ...
               [BIN_EDGES(end) + [BIN_EDGES(end) - BIN_EDGES(end-1)]/2] ];
% plot(Separation_S_deep_24hr, StructFunc_S_deep_24hr, '.', 'MarkerSize',10); hold on
plot([BIN_EDGES(1:[end-1]) + BIN_EDGES(2:end)]/2,...
     interval_avg(Separation_S_deep_24hr,StructFunc_S_deep_24hr,BIN_EDGES), '.-', 'MarkerSize',10); hold on

% % % SWOT_at_S_P
BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_EDGES = interval_avg(Separation_SWOT_moor_12hr,Separation_SWOT_moor_12hr,BIN_EDGES);
BIN_EDGES = [[ BIN_EDGES(1) - [BIN_EDGES(2) - BIN_EDGES(1)]/2 ], ...
               [BIN_EDGES(2:end) + BIN_EDGES(1:[end-1])]/2, ...
               [BIN_EDGES(end) + [BIN_EDGES(end) - BIN_EDGES(end-1)]/2] ];
% plot(Separation_SWOT_moor_12hr, StructFunc_SWOT_moor_12hr, '.', 'MarkerSize',10); hold on
plot([BIN_EDGES(1:[end-1]) + BIN_EDGES(2:end)]/2,...
     interval_avg(Separation_SWOT_moor_12hr,StructFunc_SWOT_moor_12hr,BIN_EDGES), '.-', 'MarkerSize',10); hold on

% % % SPG.600m.1hr
BIN_EDGES = [(-dBIN/2/2):[dBIN/2]:[14*dBIN/2/2]]; BIN_EDGES =  [BIN_EDGES, [BIN_EDGES(end)+dBIN]:[dBIN]:1.1]*(0.92)*DD;% initial bin edges
BIN_EDGES = interval_avg(Separation_S_P_G_600m_1hr,Separation_S_P_G_600m_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_EDGES(1) - [BIN_EDGES(2) - BIN_EDGES(1)]/2 ], ...
               [BIN_EDGES(2:end) + BIN_EDGES(1:[end-1])]/2, ...
               [BIN_EDGES(end) + [BIN_EDGES(end) - BIN_EDGES(end-1)]/2] ];
% % % plot(Separation_S_P_G_600m_1hr, StructFunc_S_P_G_600m_1hr, '.', 'MarkerSize',10); hold on
plot([BIN_EDGES(1:[end-1]) + BIN_EDGES(2:end)]/2,...
     interval_avg(Separation_S_P_G_600m_1hr,StructFunc_S_P_G_600m_1hr,BIN_EDGES), '.-', 'MarkerSize',10); hold on


legend('1hr (S+P 600m)', '24hr (S+P 600m)', ...
       '1hr (S 4400m)', '24hr (S 4400m)', ...
       'SWOT at S+P (~12hr)',...
       '1hr (S+P + Gliders 600m)')

xlabel('Distance (km)')
ylabel('2nd order S.F. (m^2)')

error(['forced stop'])





%% Plot the most recently calculated structure function (testing purposes)
close all
figure('Color',[1 1 1])
dBIN = 0.1;
% BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges (for all 10)
% BIN_EDGES = [(-3*dBIN/2):[3*dBIN]:1.1]*(0.92)*DD; % initial bin edges (for only the four S)
BIN_EDGES = [(-dBIN/4):[dBIN/2]:[16*dBIN/4]]; BIN_EDGES =  [BIN_EDGES, [BIN_EDGES(end)+dBIN]:[dBIN]:1.1]*(0.92)*DD; % initial bin edges (for moorings and gliders together)

% BIN_EDGES = [(-dBIN/2):[dBIN]:1.2]*(0.92)*DD; 

BIN_EDGES = interval_avg(SEP_VEC,SEP_VEC,BIN_EDGES);
BIN_EDGES = [[ BIN_EDGES(1) - [BIN_EDGES(2) - BIN_EDGES(1)]/2 ], ...
               [BIN_EDGES(2:end) + BIN_EDGES(1:[end-1])]/2, ...
               [BIN_EDGES(end) + [BIN_EDGES(end) - BIN_EDGES(end-1)]/2] ];

plot(SEP_VEC, SF_VEC, '.', 'MarkerSize',10); hold on
plot([BIN_EDGES(1:[end-1]) + BIN_EDGES(2:end)]/2,...
     interval_median(SEP_VEC,SF_VEC,BIN_EDGES), '.-', 'MarkerSize',10)

interval_avg(SEP_VEC,SF_VEC,BIN_EDGES)

% % % EXPERIMENTAL: eliminate outliers
%     STD_FAC = 1.5;
% plot(SEP_VEC(SF_VEC < [mean(SF_VEC) + STD_FAC*std(SF_VEC)]), SF_VEC(SF_VEC < [mean(SF_VEC) + STD_FAC*std(SF_VEC)]), '.', 'MarkerSize',10); hold on
% plot([BIN_EDGES(1:[end-1]) + BIN_EDGES(2:end)]/2,...
%      interval_avg(SEP_VEC(SF_VEC < [mean(SF_VEC) + STD_FAC*std(SF_VEC)]),SF_VEC(SF_VEC < [mean(SF_VEC) + STD_FAC*std(SF_VEC)]),BIN_EDGES), '.-', 'MarkerSize',10)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Non-Publication Figure %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL TOGETHER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
figure('Color',[1 1 1])
dBIN = 0.1;
LEG_TEXT = {};


% % % SP.600m.1hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (600m, 1hr)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_1hr,Separation_S_P_600m_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_P_600m_1hr, StructFunc_S_P_600m_1hr, '.k', 'MarkerSize',10); hold on
plot(BIN_CENTR,... plot([BIN_EDGES(1:[end-1]) + BIN_EDGES(2:end)]/2,...
     interval_avg(Separation_S_P_600m_1hr,StructFunc_S_P_600m_1hr,BIN_EDGES), '.-', 'MarkerSize',25,'LineWidth',2); hold on


% % SP.600m.24hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (600m, 24hr)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_24hr,Separation_S_P_600m_24hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_P_600m_24hr, StructFunc_S_P_600m_24hr, 'k.', 'MarkerSize',10); hold on
plot(BIN_CENTR,...
     interval_avg(Separation_S_P_600m_24hr,StructFunc_S_P_600m_24hr,BIN_EDGES), '.-', 'MarkerSize',25,'LineWidth',2); hold on


% % % S.4400m.1hr
LEG_TEXT = {LEG_TEXT{:},'S Moorings (4400m, 1hr)'};
BIN_EDGES = [(-3*dBIN/2):[3*dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_deep_1hr,Separation_S_deep_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_deep_1hr, StructFunc_S_deep_1hr, '.', 'MarkerSize',10); hold on
plot(BIN_CENTR,...
     interval_avg(Separation_S_deep_1hr,StructFunc_S_deep_1hr,BIN_EDGES), '.-', 'MarkerSize',25,'LineWidth',2); hold on


% % % S.4400m.24hr
LEG_TEXT = {LEG_TEXT{:},'S Moorings (4400m, 24hr)'};
BIN_EDGES = [(-3*dBIN/2):[3*dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_deep_24hr,Separation_S_deep_24hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_deep_24hr, StructFunc_S_deep_24hr, '.', 'MarkerSize',10); hold on
plot(BIN_CENTR,...
     interval_avg(Separation_S_deep_24hr,StructFunc_S_deep_24hr,BIN_EDGES), '.-', 'MarkerSize',25,'LineWidth',2); hold on



% % % SPG.600m.1hr
LEG_TEXT = {LEG_TEXT{:},'Moorings+Gliders (600m, 1hr)'};
BIN_EDGES = [(-dBIN/2/2):[dBIN/2]:[14*dBIN/2/2]]; BIN_EDGES =  [BIN_EDGES, [BIN_EDGES(end)+dBIN]:[dBIN]:1.1]*(0.92)*DD;% initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_G_600m_1hr,Separation_S_P_G_600m_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_P_G_600m_1hr, StructFunc_S_P_G_600m_1hr, '.k', 'MarkerSize',10); hold on
plot(BIN_CENTR,...
     interval_avg(Separation_S_P_G_600m_1hr,StructFunc_S_P_G_600m_1hr,BIN_EDGES), '.-', 'MarkerSize',25,'LineWidth',2); hold on



% % % SPG.600m.24hr
LEG_TEXT = {LEG_TEXT{:},'Moorings+Gliders (600m, 24hr)'};
BIN_EDGES = [(-dBIN/2/2):[dBIN/2]:[14*dBIN/2/2]]; BIN_EDGES =  [BIN_EDGES, [BIN_EDGES(end)+dBIN]:[dBIN]:1.1]*(0.92)*DD;% initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_G_600m_24hr,Separation_S_P_G_600m_24hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% % % plot(Separation_S_P_G_600m_24hr, StructFunc_S_P_G_600m_24hr, '.', 'MarkerSize',10); hold on
plot(BIN_CENTR,...
     interval_avg(Separation_S_P_G_600m_24hr,StructFunc_S_P_G_600m_24hr,BIN_EDGES), '.-', 'MarkerSize',25,'LineWidth',2); hold on


% % % SWOT_at_S_P_G
LEG_TEXT = {LEG_TEXT{:},'SWOT cal/val at Moorings+Gliders'};
BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_SWOT_moor_12hr,Separation_SWOT_moor_12hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_SWOT_moor_12hr, StructFunc_SWOT_moor_12hr, '.', 'MarkerSize',10); hold on
plot(BIN_CENTR,...
     interval_avg(Separation_SWOT_moor_12hr,StructFunc_SWOT_moor_12hr,BIN_EDGES), '.-', 'MarkerSize',25,'LineWidth',2); hold on




% % % LEG = legend('Moorings (600m, 1hr)', 'Moorings (600m, 24hr)', ...
% % %              'S Moorings (4400m, 1hr)', 'S Moorings (4400m, 24hr)', ...
% % %              'Moorings+Gliders (600m, 1hr)',...
% % %              'SWOT cal/val at Moorings');
LEG = legend(LEG_TEXT);
LEG.FontSize = 14;
LEG.Position = [0.6547 0.2546 0.1526 0.1130];

xlabel('Binned Separation (km)','Interpreter','latex','FontSize',16)
ylabel('$$\overline{\big(\eta(x) - \eta(x+s)\big)^2}$$ \quad (m$^2$)','Interpreter','latex','FontSize',16)
set(gca,'FontSize',14)

grid on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Publication Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PANELED SF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('This is now a supplemental figure; see ./SF_plot_fixed_bottom_atSWOTtimes.m')

close all
figure('Color',[1 1 1])
dBIN = 0.1;
LEG_TEXT = {};
CO = colororder;
m2_to_cm2 = 100*100; % =10^4 if you are to plot cm^2; =1 if m^2
XSCALE = 'lin'; YSCALE = 'lin';

tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');

nexttile;%subplot(121) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % SWOT_at_S_P_G
% LEG_TEXT = {LEG_TEXT{:},['SWOT cal/val at' newline 'Moorings+Gliders']};
LEG_TEXT = {LEG_TEXT{:},['SWOT cal/val at Moorings+Gliders']};
BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
% BIN_EDGES = [(-dBIN/2/2):[dBIN/2]:[14*dBIN/2/2]]; BIN_EDGES =  [BIN_EDGES, [BIN_EDGES(end)+dBIN]:[dBIN]:1.1]*(0.92)*DD;% initial bin edges
BIN_CENTR = interval_avg(Separation_SWOT_moor_12hr,Separation_SWOT_moor_12hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_SWOT_moor_12hr, StructFunc_SWOT_moor_12hr, '.', 'MarkerSize',10); hold on
[Mean_SF_Separation_SWOT_moor_12hr] = ...
    errorbar_stdofmean(BIN_CENTR, Separation_SWOT_moor_12hr, m2_to_cm2*StructFunc_SWOT_moor_12hr, BIN_EDGES, [1 1 1]*0.7); hold on


% % % SP.600m.1hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (600m, 1hr)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_1hr,Separation_S_P_600m_1hr,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_P_600m_1hr, StructFunc_S_P_600m_1hr, '.k', 'MarkerSize',10); hold on
[Mean_SF_Separation_S_P_600m_1hr] = ...
    errorbar_stdofmean(BIN_CENTR, Separation_S_P_600m_1hr, m2_to_cm2*StructFunc_S_P_600m_1hr, BIN_EDGES, CO(1,:)); hold on
BIN_CENTR_S_P = BIN_CENTR;


% % SP.600m.24hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (600m, 24hr)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_24hr,Separation_S_P_600m_24hr,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_P_600m_24hr, StructFunc_S_P_600m_24hr, 'k.', 'MarkerSize',10); hold on
[Mean_SF_Separation_S_P_600m_24hr] = ...
    errorbar_stdofmean(BIN_CENTR, Separation_S_P_600m_24hr, m2_to_cm2*StructFunc_S_P_600m_24hr, BIN_EDGES, CO(2,:));
set(gca,'FontSize',18); grid on


% % % SPG.600m.1hr
LEG_TEXT = {LEG_TEXT{:},'Moorings+Gliders (600m, 1hr)'};
BIN_EDGES = [(-dBIN/2/2):[dBIN/2]:[14*dBIN/2/2]]; BIN_EDGES =  [BIN_EDGES, [BIN_EDGES(end)+dBIN]:[dBIN]:1.1]*(0.92)*DD;% initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_G_600m_1hr,Separation_S_P_G_600m_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_P_G_600m_1hr, StructFunc_S_P_G_600m_1hr, '.k', 'MarkerSize',10); hold on
[Mean_SF_Separation_S_P_G_600m_1hr] = ...
    errorbar_stdofmean(BIN_CENTR, Separation_S_P_G_600m_1hr, m2_to_cm2*StructFunc_S_P_G_600m_1hr, BIN_EDGES, CO(3,:)); hold on


% % % SP.600m.1hr.noP1
LEG_TEXT = {LEG_TEXT{:},'Moorings, no P1 (600m, 1hr)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_1hr_noP1,Separation_S_P_600m_1hr_noP1,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_P_600m_1hr_noP1, StructFunc_S_P_600m_1hr_noP1, '.k', 'MarkerSize',10); hold on
[Mean_SF_Separation_S_P_600m_1hr_noP1] = ...
    errorbar_stdofmean(BIN_CENTR, Separation_S_P_600m_1hr_noP1, m2_to_cm2*StructFunc_S_P_600m_1hr_noP1, BIN_EDGES, CO(1,:), '--'); hold on


% % SP.600m.24hr.noP1
LEG_TEXT = {LEG_TEXT{:},'Moorings, no P1 (600m, 24hr)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_24hr_noP1,Separation_S_P_600m_24hr_noP1,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_P_600m_24hr_noP1, StructFunc_S_P_600m_24hr_noP1, 'k.', 'MarkerSize',10); hold on
[Mean_SF_Separation_S_P_600m_24hr_noP1] = ...
    errorbar_stdofmean(BIN_CENTR, Separation_S_P_600m_24hr_noP1, m2_to_cm2*StructFunc_S_P_600m_24hr_noP1, BIN_EDGES, CO(2,:), '--');
set(gca,'FontSize',18); grid on
set(gca,'yscale',YSCALE,'xscale',XSCALE)

xlabel('Separation (km)')
ylabel('$$\overline{\big(\eta''(x) - \eta''(x+s)\big)^2}$$ \quad (cm$^2$)','Interpreter','latex','FontSize',16)

set(gca,'YLim',[0 50])
set(gca,'FontSize',20); grid on
LEG = legend(LEG_TEXT,'Location','southeast');
% LEG.Position = LEG.Position + [-0.1 0 0 0];
LEG_TEXT = {};


nexttile;%subplot(122) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % SWOT_at_S
LEG_TEXT = {LEG_TEXT{:},'SWOT cal/val at S-Moorings'};
BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_SWOTatS_moor_12hr,Separation_SWOTatS_moor_12hr,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_SWOTatS_moor_12hr, StructFunc_SWOTatS_moor_12hr, '.', 'MarkerSize',10); hold on
[Mean_SF_Separation_SWOTatS_moor_12hr] = ...
    errorbar_stdofmean(BIN_CENTR, Separation_SWOTatS_moor_12hr, m2_to_cm2*StructFunc_SWOTatS_moor_12hr, BIN_EDGES, [1 1 1]*0.7); hold on


% % % S.4400m.1hr
LEG_TEXT = {LEG_TEXT{:},'S Moorings (1800m, 1hr)'};
BIN_EDGES = [(-3*dBIN/2):[3*dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_deep_1hr,Separation_S_deep_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_deep_1hr, StructFunc_S_deep_1hr, '.', 'MarkerSize',10); hold on
[Mean_SF_Separation_S_deep_1hr] = ...
    errorbar_stdofmean(BIN_CENTR, Separation_S_deep_1hr, m2_to_cm2*StructFunc_S_deep_1hr, BIN_EDGES, CO(1,:)); hold on
BIN_CENTR_S = BIN_CENTR;


% % % S.4400m.24hr
LEG_TEXT = {LEG_TEXT{:},'S Moorings (1800m, 24hr)'};
BIN_EDGES = [(-3*dBIN/2):[3*dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_deep_24hr,Separation_S_deep_24hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_deep_24hr, StructFunc_S_deep_24hr, '.', 'MarkerSize',10); hold on
[Mean_SF_Separation_S_deep_24hr] = ...
    errorbar_stdofmean(BIN_CENTR, Separation_S_deep_24hr, m2_to_cm2*StructFunc_S_deep_24hr, BIN_EDGES, CO(2,:)); hold on
% xlabel('Binned Separation (km)','Interpreter','latex','FontSize',16)
% ylabel('$$\overline{\big(\eta(x) - \eta(x+s)\big)^2}$$ \quad (m$^2$)','Interpreter','latex','FontSize',16)


% % % S.600m_1hr
LEG_TEXT = {LEG_TEXT{:},'S Moorings (600m, 1hr)'};
BIN_EDGES = [(-3*dBIN/2):[3*dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_600m_1hr,Separation_S_600m_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_S_600m_1hr, StructFunc_S_600m_1hr, '.', 'MarkerSize',10); hold on
[Mean_SF_Separation_S_600m_1hr] = ...
    errorbar_stdofmean(BIN_CENTR, Separation_S_600m_1hr, m2_to_cm2*StructFunc_S_600m_1hr, BIN_EDGES, CO(1,:), '.--'); hold on

set(gca,'YLim',[0 50])
set(gca,'yscale',YSCALE,'xscale',XSCALE)
set(gca,'FontSize',20); grid on
LEG = legend(LEG_TEXT,'Location','southeast');
% LEG.Position = LEG.Position + [-0.1 0 0 0];
LEG_TEXT = {}; %nexttile;%subplot(224)

set(gcf,'Position', [54         195        1291         505])

INPUT = input(['Do you want to save this figure? Enter any number for "yes" or\n push ' ...
    '"Enter" with a blank input for "No".\n']);

if isempty(INPUT)
    disp(['Image not saved'])
else
    figure(1)
    exportgraphics(gcf,...
        '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/extra_FigD.pdf',...
        'BackgroundColor','none','ContentType','vector')
    disp(['Image saved'])
end

%% Figure showing the differences in structure function values as a function
% of separation:
close all
figure('Color',[1 1 1])

plot(BIN_CENTR_S_P,Mean_SF_Separation_S_P_600m_1hr - Mean_SF_Separation_S_P_600m_24hr, ...
    '.-', 'MarkerSize',25,'LineWidth',2); hold on
plot(BIN_CENTR_S,Mean_SF_Separation_S_deep_1hr - Mean_SF_Separation_S_deep_24hr, ...
    '.-', 'MarkerSize',25,'LineWidth',2)

title('') % title('Mean\_SF\_Separation\_S\_P\_600m\_1hr-Mean\_SF\_Separation\_S\_P\_600m\_24hr')
xlabel('Separation (km)')
ylabel('$D_\eta(s)_{1\mathrm{hr}}$ - $D_\eta(s)_{24\mathrm{hr}}$ \quad (cm$^2$)','Interpreter','latex','FontSize',16)
legend('600 m (S & P)','1800 m (S)')

set(gca,'YLim',[0 6])
set(gca,'yscale',YSCALE,'xscale',XSCALE)
set(gca,'FontSize',20); grid on

INPUT = input(['Do you want to save this figure? Enter any number for "yes" or\n push ' ...
    '"Enter" with a blank input for "No".\n']);

if isempty(INPUT)
    disp(['Image not saved'])
else
    figure(1)
    % exportgraphics(gcf,...
    %     '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/F_SF_diff.pdf',...
    %     'BackgroundColor','none','ContentType','vector')
    exportgraphics(gcf,...
    '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/extra_FigC.pdf',...
    'BackgroundColor','none','ContentType','vector')
    disp(['Image saved'])
end

% exportgraphics(gcf,...
%     '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/extra_FigC.pdf',...
%     'BackgroundColor','none','ContentType','vector')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Log plot of S.F. of SWOT %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure('Color',[1 1 1])
dBIN = 0.1;
LEG_TEXT = {};
CO = colororder;
m2_to_cm2 = 100*100; % =10^4 if you are to plot cm^2; =1 if m^2
XSCALE = 'lin'; YSCALE = 'lin';

tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');

nexttile;%subplot(121) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % SWOT_ALL
% % % This structure function must be created by running the last three
% % % sections of <./WavenumberSpectrum_SWOT.m>
warning('This structure function must be created by running the last three sections of <./WavenumberSpectrum_SWOT.m>')
% LEG_TEXT = {LEG_TEXT{:},'SWOT cal/val all data (29.9N-50.2N)'};
LEG_TEXT = {LEG_TEXT{:},['SWOT cal/val all data (' num2str(round(min([SWOT.lat{end}(:) ; SWOT.lat{end-1}(:)]),1)) 'N-' ...
                                                   num2str(round(max([SWOT.lat{end}(:) ; SWOT.lat{end-1}(:)]),1)) 'N)']};
fill([flip(DIST_VEC); DIST_VEC], ...
     [flip(mean(SF_MAT_swot,2,'omitnan') + std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))); ...
           mean(SF_MAT_swot,2,'omitnan') - std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))], ...
     [1 1 1]*0.7, 'EdgeAlpha',0,'HandleVisibility','off'); hold on
plot(DIST_VEC, mean(SF_MAT_swot,2,'omitnan')*m2_to_cm2, 'k.-'); hold on
Mean_SF_MAT_swot = mean(SF_MAT_swot,2,'omitnan')*m2_to_cm2;


% % % SWOT_at_S_P_G
LEG_TEXT = {LEG_TEXT{:},'SWOT cal/val at Moorings+Gliders'};
BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_SWOT_moor_12hr,Separation_SWOT_moor_12hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_SWOT_moor_12hr, StructFunc_SWOT_moor_12hr, '.', 'MarkerSize',10); hold on
errorbar_stdofmean(BIN_CENTR, Separation_SWOT_moor_12hr, m2_to_cm2*StructFunc_SWOT_moor_12hr, BIN_EDGES, [1 1 1]*0.7); hold on
Mean_SWOT_at_S_P_G = interval_avg(Separation_SWOT_moor_12hr, m2_to_cm2*StructFunc_SWOT_moor_12hr, BIN_EDGES);
% plot(BIN_CENTR,Mean_SWOT_at_S_P_G,'r') % verify that this is the quantity of interest

% % % % SWOT_at_S_P_G
% LEG_TEXT = {LEG_TEXT{:},['SWOT cal/val at Moorings+Gliders' newline 'and three other crossover ' newline 'diamonds']};
% BIN_EDGES = [(-dBIN/2):[dBIN]:1.1 , (1.1+2*dBIN):(2*dBIN):(2.0) , (2.0+5*dBIN):(5*dBIN):(5.5)]*(0.92)*DD; % initial bin edges
% BIN_CENTR = interval_avg(Separation_SWOT_extended_12hr,Separation_SWOT_extended_12hr,BIN_EDGES);
% BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
%                [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
%                [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% % plot(Separation_SWOT_extended_12hr, StructFunc_SWOT_extended_12hr, '.', 'MarkerSize',10); hold on
% errorbar_stdofmean(BIN_CENTR, Separation_SWOT_extended_12hr, m2_to_cm2*StructFunc_SWOT_extended_12hr, BIN_EDGES, [1 1 1]*0.3); hold on

xlabel('Separation (km)')
ylabel('$\overline{\big(\mathrm{SSHA}(x) - \mathrm{SSHA}(x+s)\big)^2}$ \quad (cm$^2$)','Interpreter','latex','FontSize',16)

set(gca,'XLim',DIST_VEC([1 end]))
set(gca,'YLim',[0.1 100]); set(gca,'XLim',[0 1000])
set(gca,'yscale',YSCALE,'xscale',XSCALE)
set(gca,'yscale','log','xscale','log')
set(gca,'FontSize',20); grid on
LEG = legend(LEG_TEXT,'Location','southeast');
% LEG.Position = LEG.Position + [-0.1 0 0 0];
LEG_TEXT = {};


% Fitting parameters
FitBasisFunction = {'a*t.^b','a','b'};
FitParameters0 = [0.1,1.5];
FitBasisDerivatives = {'t.^b','a*log(t).*t.^b'};
ToleranceLevel = 0.0001; % NLLSF parameter, tells change in parameter magnitude, if small then the result was stable
N_nllsf_iterations = 20; % number of iterations before we give up


nexttile;%subplot(122) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SWOT_fitted_SpecSlopes','var')
    % i.e. only run this if we haven't yet calculated the slope as a
    % function of maximum fitting interval.
    Nx_max = 100; % Maximum number of lags to fit out to (approx. = 2 km each integer)
    IND_fitting_lags = [4:2:10 12:4:20 25:5:50 60:10:100 150]';
    SWOT_fitted_SpecSlopes            = nan(size(IND_fitting_lags));
    SWOT_fitted_SpecSlopes_covariance = nan(size(IND_fitting_lags));
    for xi = 1:length(IND_fitting_lags)
        XI = IND_fitting_lags(xi);
        [NLLSF_COEF,JACOBIAN,SF_FIT,ConvergenceRecord] = ...
        nonlinear_lsqf(mean(SF_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2, DIST_VEC(1:XI), ...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations);
        NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ [mean(SF_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2] - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
        SWOT_fitted_SpecSlopes(xi) = -NLLSF_COEF(2) - 1;
        SWOT_fitted_SpecSlopes_covariance(xi) = NLLSF_COEF_cov(2,2);
        disp(ConvergenceRecord(end))
    end
else
end; disp(' ')

plot(DIST_VEC(IND_fitting_lags), SWOT_fitted_SpecSlopes, 'k.-', 'MarkerSize',20,'LineWidth',1.5); hold on
errorbar(DIST_VEC(IND_fitting_lags), SWOT_fitted_SpecSlopes, SWOT_fitted_SpecSlopes_covariance, 'MarkerSize',20,'LineWidth',1.5,'Color','k')

% Do the same (get slope as a function of max. evaluated separation) for
% SWOT at instrument locations:
SWOT_CalVal_fitted_SpecSlopes            = nan(size(BIN_CENTR));
SWOT_CalVal_fitted_SpecSlopes_maxlag     = nan(size(BIN_CENTR));
SWOT_CalVal_fitted_SpecSlopes_covariance = nan(size(BIN_CENTR));
for xi = 3:length(BIN_CENTR)
    XI = IND_fitting_lags(xi);
    [NLLSF_COEF,JACOBIAN,SF_FIT,ConvergenceRecord] = ...
        nonlinear_lsqf(Mean_SWOT_at_S_P_G(1:xi), BIN_CENTR(1:xi), ...
        FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations);
    NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ Mean_SWOT_at_S_P_G(1:xi)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
    SWOT_CalVal_fitted_SpecSlopes(xi) = -NLLSF_COEF(2) - 1;
    SWOT_CalVal_fitted_SpecSlopes_maxlag(xi) = BIN_CENTR(xi);
    % The covariance of gamma (SF order) is the same as that of lambda
    % (spectral slope) because both are just the other times -1 and shifted -1:
    SWOT_CalVal_fitted_SpecSlopes_covariance(xi) = NLLSF_COEF_cov(2,2);
    disp(ConvergenceRecord(end))
end
plot(SWOT_CalVal_fitted_SpecSlopes_maxlag, SWOT_CalVal_fitted_SpecSlopes, ...
    '.-', 'Color', [1 1 1]*0.7, 'MarkerSize',20,'LineWidth',1.5)
errorbar(SWOT_CalVal_fitted_SpecSlopes_maxlag, SWOT_CalVal_fitted_SpecSlopes, SWOT_CalVal_fitted_SpecSlopes_covariance, ...
    'Color', [1 1 1]*0.7, 'MarkerSize',20,'LineWidth',1.5)

xlabel({'Max. separation (km) of';'fitted structure function'})
ylabel('Inferred Spectral Slope \lambda')
set(gca,'FontSize',20); grid on

% set(gcf,'Position',[-1759 392 1466 442])
set(gcf,'Position',[1 355 1041 442])


INPUT = input(['Do you want to save this figure? Enter any number for "yes" or\n push ' ...
    '"Enter" with a blank input for "No".\n']);

if isempty(INPUT)
    disp(['Image not saved'])
else
    exportgraphics(gcf,...
        '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/F4.pdf',...
        'BackgroundColor','none','ContentType','vector')
    disp(['Image saved'])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Spectral slope of average SFs from above figure: %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% /Users/kachelein/Documents/MATLAB/Luke/GaussNewtonParameterUncertaintyTest.m
% https://www8.cs.umu.se/kurser/5DA001/HT07/lectures/lsq-handouts.pdf
% https://www.incertitudes.fr/book.pdf

clc ; close all

FitBasisFunction = {'a*t.^b','a','b'};
FitParameters0 = [0.1,1.5];
FitBasisDerivatives = {'t.^b','a*log(t).*t.^b'};
ToleranceLevel = 0.00001; % NLLSF parameter, tells change in parameter magnitude, if small then the result was stable
N_nllsf_iterations = 10; % number of iterations before we give up




% % % SP.600m.24hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (600m, 24hr)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_24hr,Separation_S_P_600m_24hr,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_P_600m_24hr, m2_to_cm2*StructFunc_S_P_600m_24hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:[end-1]),BIN_CENTR(2:[end-1]),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:[end-1])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/24hr/M to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
% [NLLSF_COEF,~,~,~] = ...
%             nonlinear_lsqf(SF_avg(5:[end-1]),BIN_CENTR(5:[end-1]),...
%             FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
% disp(['600m24hr/M 40 to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1)])
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:5),BIN_CENTR(2:5),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:5)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/24hr/M to ' num2str(BIN_CENTR(5)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')


% % % SP.600m.1hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (600m, 1hr)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_1hr,Separation_S_P_600m_1hr,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_P_600m_1hr, m2_to_cm2*StructFunc_S_P_600m_1hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:[end-1]),BIN_CENTR(2:[end-1]),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:[end-1])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/1hr/M to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
% [NLLSF_COEF,~,~,~] = ...
%             nonlinear_lsqf(SF_avg(5:[end-1]),BIN_CENTR(5:[end-1]),...
%             FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
% disp(['600m/1hr/M 40 to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1)])
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:5),BIN_CENTR(2:5),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:5)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/1hr/M to ' num2str(BIN_CENTR(5)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')


% % % SPG.600m.1hr
LEG_TEXT = {LEG_TEXT{:},'Moorings+Gliders (600m, 1hr)'};
BIN_EDGES = [(-dBIN/2/2):[dBIN/2]:[14*dBIN/2/2]]; BIN_EDGES =  [BIN_EDGES, [BIN_EDGES(end)+dBIN]:[dBIN]:1.1]*(0.92)*DD;% initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_G_600m_1hr,Separation_S_P_G_600m_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_P_G_600m_1hr, m2_to_cm2*StructFunc_S_P_G_600m_1hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:[end-1]),BIN_CENTR(2:[end-1]),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:[end-1])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/1hr/MG to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
% [NLLSF_COEF,~,~,~] = ...
%             nonlinear_lsqf(SF_avg(8:[end-1]),BIN_CENTR(8:[end-1]),...
%             FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
% disp(['600m24hr/MG 40 to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1)])
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:8),BIN_CENTR(2:8),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:8)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/1hr/MG to ' num2str(BIN_CENTR(8)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')




% % % SP.1800m.24hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (1800m, 24hr)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_deep_24hr, Separation_S_deep_24hr,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_deep_24hr, m2_to_cm2*StructFunc_S_deep_24hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:end),BIN_CENTR(2:end),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:end)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['1800m/24hr/M to ' num2str(BIN_CENTR(end)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
% 
BIN_CENTR = interval_avg(Separation_S_deep_1hr, Separation_S_deep_1hr,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_deep_1hr, m2_to_cm2*StructFunc_S_deep_1hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:end),BIN_CENTR(2:end),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:end)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['1800m/1hr/M to ' num2str(BIN_CENTR(end)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')


% % % SWOT_at_S_P_G
LEG_TEXT = {LEG_TEXT{:},'SWOT cal/val at Moorings+Gliders'};
BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_SWOT_moor_12hr,Separation_SWOT_moor_12hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_SWOT_moor_12hr, m2_to_cm2*StructFunc_SWOT_moor_12hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:[end-1]),BIN_CENTR(2:[end-1]),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:[end-1])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
% The covariance of gamma (SF order) is the same as that of lambda
% (spectral slope) because both are just the other times -1 and shifted -1.
disp(['SWOT to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
% [NLLSF_COEF,~,~,~] = ...
%             nonlinear_lsqf(SF_avg(5:[end-1]),BIN_CENTR(5:[end-1]),...
%             FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
% disp(['SWOT 40 to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1)])
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:5),BIN_CENTR(2:5),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:5)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['SWOT to ' num2str(BIN_CENTR(5)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')

[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(Mean_SF_MAT_swot(2:dsearchn(DIST_VEC,90)), DIST_VEC(2:dsearchn(DIST_VEC,90)),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ Mean_SF_MAT_swot(2:dsearchn(DIST_VEC,90)) - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['SWOT (all) to ' num2str(DIST_VEC(dsearchn(DIST_VEC,90))) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(Mean_SF_MAT_swot(2:dsearchn(DIST_VEC,40)), DIST_VEC(2:dsearchn(DIST_VEC,40)),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ Mean_SF_MAT_swot(2:dsearchn(DIST_VEC,40)) - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['SWOT (all) to ' num2str(DIST_VEC(dsearchn(DIST_VEC,40))) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')




%% Time series of evolving spectral slopes from structure functions

% close all
% 
% error('Do not rerun unless you are sure (takes ~1 min 20 sec)')
% 
% % SP.600m.Xhr
% SF_TIME_STEP_LIST = [1 6 12 24]./24; % days
% SF_TIME_SERIES = struct;
% ToleranceLevel = 0.001; % NLLSF parameter, tells change in parameter magnitude, if small then the result was stable
% N_nllsf_iterations = 10; % number of iterations before we give up
% 
% for ii=1:length(SF_TIME_STEP_LIST)
%     SF_TIME_STEP = SF_TIME_STEP_LIST(ii); % days
%     TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
% 
%     DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
%     SEP_VEC = [];
%     SF_VEC = [];
%     Slope_from_SF = [];
% 
%     tic
%     warning('off','all')
%     for ti = 1:[length(TIME_BIN_EDGES) - 1]
%         DATA_in_interval = [];
%         XY_in_interval = [];
%         SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
%         for mi = 1:length(MM_list)
%             MM = MM_list{mi};
%             if strcmp(MM(1),'P') %& ~[strcmp(MM,'P1')]
%                 P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
%                     '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
%                     '~SPIKE_IND{mi}']); % <-- Spike removal
%                 DATA_in_interval = [DATA_in_interval ; ...
%                     mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
%                 XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
%                     mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
%             elseif strcmp(MM(1),'S')
%                 P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
%                     '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
%                     '~SPIKE_IND{mi}']); % <-- Spike removal
%                 DATA_in_interval = [DATA_in_interval ; ...
%                     mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
%                 XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
%                     mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
%             else
%             end
%         end
%         SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
%             [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
%         SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;
% 
%         SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
%         SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];
% 
%         % % % Average this time-interval's structure function and do a NLLSF
%         BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
%         BIN_CENTR = interval_avg(SEP_VEC,SEP_VEC,BIN_EDGES);
%         BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
%             [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
%             [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
%         % MEAN_SF_ti = interval_avg(SEP_VEC,SF_VEC,BIN_EDGES); % figure;plot(BIN_CENTR,MEAN_SF_ti,'.-')
%         MEAN_SF_ti = interval_avg(SEP_MAT_in_interval(:),SF_MAT_in_interval(:),BIN_EDGES); % figure;plot(BIN_CENTR,MEAN_SF_ti,'.-')
% 
%         % NNLSF:
%         if sum(~isfinite(MEAN_SF_ti)) < 7 % check for enough data
%             [NLLSF_COEF,~,~,ConvergenceRecord] = ...
%                 nonlinear_lsqf(MEAN_SF_ti(2:[end-1]),BIN_CENTR(2:[end-1]),...
%                 {'a*t.^b','a','b'},...
%                 [0.00001,1.5],{'t.^b','a*log(t).*t.^b'},ToleranceLevel,   N_nllsf_iterations);
%             if ConvergenceRecord(end) < ToleranceLevel % good result (converged)
%                 Slope_from_SF = [Slope_from_SF ; (-1-NLLSF_COEF(2))];
%             else % bad result (didn't converge)
%                 Slope_from_SF = [Slope_from_SF ; NaN];
%             end
%         else
%             Slope_from_SF = [Slope_from_SF ; NaN];
%         end
% 
%         if ~mod(ti,10)
%             disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
%         end
%     end
% 
%     eval(['SF_TIME_SERIES.Step_' num2str(SF_TIME_STEP*24) 'hr.TimeStep = SF_TIME_STEP;'])
%     eval(['SF_TIME_SERIES.Step_' num2str(SF_TIME_STEP*24) 'hr.Time = [TIME_BIN_EDGES(1:[end-1]) + TIME_BIN_EDGES(2:end)]/2;'])
%     eval(['SF_TIME_SERIES.Step_' num2str(SF_TIME_STEP*24) 'hr.Slope_from_SF = Slope_from_SF;'])
% 
%     warning('on','all')
%     toc
% end
% figure
% BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
% % BIN_CENTR = interval_avg(SEP_VEC,SEP_VEC,BIN_EDGES);
% % BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
% %                [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
% %                [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(SEP_VEC, m2_to_cm2*SF_VEC, '.', 'MarkerSize',10); hold on
% errorbar_stdofmean(BIN_CENTR, SEP_VEC, m2_to_cm2*SF_VEC, BIN_EDGES, [1 1 1]*0.7); hold on

%% Time series of evolving spectral slopes, centered at SWOT times

close all

% SP.600m.Xhr
SF_TIME_STEP_LIST = [1 6 12 24]./24; % days
if exist('SF_TIME_SERIES','var')
else
    SF_TIME_SERIES = struct;
end
ToleranceLevel = 0.001; % NLLSF parameter, tells change in parameter magnitude, if small then the result was stable
N_nllsf_iterations = 10; % number of iterations before we give up

for ti = 1:length(SWOT.time)
    median_SWOT_time(ti) = median(SWOT.time{ti},'omitnan');
end

for ii=1:length(SF_TIME_STEP_LIST)
    SF_TIME_STEP = SF_TIME_STEP_LIST(ii); % days

    DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
    SEP_VEC = [];
    SF_VEC = [];
    Slope_from_SF = [];

    tic
    warning('off','all')
    for ti = 1:length(SWOT.time)
        DATA_in_interval = [];
        XY_in_interval = [];
        SF_TIME_INTERVAL = median(SWOT.time{ti},'omitnan') + [-1 1]*SF_TIME_STEP/2 - T0;
        for mi = 1:length(MM_list)
            MM = MM_list{mi};
            if strcmp(MM(1),'P') %& ~[strcmp(MM,'P1')]
                P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                    '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                    '~SPIKE_IND{mi}']); % <-- Spike removal
                DATA_in_interval = [DATA_in_interval ; ...
                    mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
                XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                    mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
            elseif strcmp(MM(1),'S')
                P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
                    '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                    '~SPIKE_IND{mi}']); % <-- Spike removal
                DATA_in_interval = [DATA_in_interval ; ...
                    mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
                XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
                    mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
            else
            end
        end
        SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
            [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
        SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

        SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
        SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

        % % % Average this time-interval's structure function and do a NLLSF
        BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
        BIN_CENTR = interval_avg(SEP_VEC,SEP_VEC,BIN_EDGES);
        BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
            [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
            [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
        % MEAN_SF_ti = interval_avg(SEP_VEC,SF_VEC,BIN_EDGES); % figure;plot(BIN_CENTR,MEAN_SF_ti,'.-')
        MEAN_SF_ti = interval_avg(SEP_MAT_in_interval(:),SF_MAT_in_interval(:),BIN_EDGES); % figure;plot(BIN_CENTR,MEAN_SF_ti,'.-')

        % NNLSF:
        if sum(~isfinite(MEAN_SF_ti)) < 7 % check for enough data
            [NLLSF_COEF,~,~,ConvergenceRecord] = ...
                nonlinear_lsqf(MEAN_SF_ti(2:[end-1]),BIN_CENTR(2:[end-1]),...
                {'a*t.^b','a','b'},...
                [0.00001,1.5],{'t.^b','a*log(t).*t.^b'},ToleranceLevel,   N_nllsf_iterations);
            if ConvergenceRecord(end) < ToleranceLevel % good result (converged)
                Slope_from_SF = [Slope_from_SF ; (-1-NLLSF_COEF(2))];
            else % bad result (didn't converge)
                Slope_from_SF = [Slope_from_SF ; NaN];
            end
        else
            Slope_from_SF = [Slope_from_SF ; NaN];
        end

        if ~mod(ti,10)
            disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
        end
    end

    eval(['SF_TIME_SERIES.Step_' num2str(SF_TIME_STEP*24) 'hr.TimeStep = SF_TIME_STEP;'])
    eval(['SF_TIME_SERIES.Step_' num2str(SF_TIME_STEP*24) 'hr.Time = median_SWOT_time;'])
    eval(['SF_TIME_SERIES.Step_' num2str(SF_TIME_STEP*24) 'hr.Slope_from_SF = Slope_from_SF;'])

    warning('on','all')
    toc
end
figure
BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
% BIN_CENTR = interval_avg(SEP_VEC,SEP_VEC,BIN_EDGES);
% BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
%                [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
%                [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
plot(SEP_VEC, m2_to_cm2*SF_VEC, '.', 'MarkerSize',10); hold on
errorbar_stdofmean(BIN_CENTR, SEP_VEC, m2_to_cm2*SF_VEC, BIN_EDGES, [1 1 1]*0.7); hold on

%% Time series of evolving spectral slopes from SWOT SF's at mooring location:
% Set vector of nominal locations of moorings:
M_lonlat_nominal = nan(length(MM_list),2);
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    M_lonlat_nominal(mi,:) = [eval([MM '.P.LONGITUDE']) , eval([MM '.P.LATITUDE'])];
end

MM_GG_list = {'P1','P2','P3',    'P5','P6','P7',    'S1','S2','S3','S4',     'ru32_600','ru38_600'};

SEP_VEC = [];
SF_VEC = [];
Slope_from_SF = [];
Time_for_SF = [];
tic
for ti = 1:length(SWOT.time)

    SWOT_M_time = [];
    SWOT_M_ssha = [];
    SWOT_M_ssha_qc = [];
    SWOT_M_rollerror = [];
    SWOT_M_rollerror_qc = [];
    SWOT_M_hret = [];
    SWOT_M_ll = [];
    % for mi = 1:length(MM_list)
        % MM = MM_list{mi};
    for mi = 1:length(MM_GG_list)
        MM = MM_GG_list{mi};
        if strcmp(MM(1),'S') || strcmp(MM(1),'P') % at mooring location
            MM_GPS_LON = eval([[MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY']]);
            MM_GPS_LAT = eval([[MM '.GPS.LATITUDE_GPS_SURFACE_BUOY']]);
            Mooring_Loc_Mat = abs(SWOT.lon{ti} + 1i*SWOT.lat{ti} - [...
                MM_GPS_LON(dsearchn(eval([MM '.GPS.TIME_GPS_SURFACE_BUOY']) + T0, median(SWOT.time{ti}(:)))) + ...
                MM_GPS_LAT(dsearchn(eval([MM '.GPS.TIME_GPS_SURFACE_BUOY']) + T0, median(SWOT.time{ti}(:))))*1i] );
            % ^ the use of median(SWOT.time{ti}(:)) will mean that the true
            % time will be off by a maximum of 3-6 minutes, which should
            % be ok considering the time steps of profilers is 0.5-4 hr
        elseif strcmp(MM(1),'r') % at Rutgers glider location
            GG_GPS_LON = eval([[MM '.LONGITUDE_PROFILE']]);
            GG_GPS_LAT = eval([[MM '.LATITUDE_PROFILE']]);
            Mooring_Loc_Mat = abs(SWOT.lon{ti} + 1i*SWOT.lat{ti} - [...
                GG_GPS_LON(dsearchn(eval([MM '.TIME_STERIC_HEIGHT']) + T0, median(SWOT.time{ti}(:)))) + ...
                GG_GPS_LAT(dsearchn(eval([MM '.TIME_STERIC_HEIGHT']) + T0, median(SWOT.time{ti}(:))))*1i] );
        else
        end
        SWOT_M_time_ti = SWOT.time{ti}( ~~sum(Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)),1)' );
        SWOT_M_time = [SWOT_M_time ; SWOT_M_time_ti];

        SWOT_M_ssha_ti = SWOT.ssha_karin_2{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_ssha = [SWOT_M_ssha ; SWOT_M_ssha_ti];

        SWOT_M_ssha_qc_ti = SWOT.ssha_karin_2_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_ssha_qc = [SWOT_M_ssha_qc ; SWOT_M_ssha_qc_ti];

        SWOT_M_rollerror_ti = SWOT.height_cor_xover{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_rollerror = [SWOT_M_rollerror ; SWOT_M_rollerror_ti];

        SWOT_M_rollerror_qc_ti = SWOT.height_cor_xover_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_rollerror_qc = [SWOT_M_rollerror_qc ; SWOT_M_rollerror_qc_ti];

        SWOT_M_hret_ti = SWOT.internal_tide_hret{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
        SWOT_M_hret = [SWOT_M_hret ; SWOT_M_hret_ti];

        SWOT_M_ll_ti = SWOT.lon{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) ) + ...
                       SWOT.lat{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) )*1i;
        SWOT_M_ll = [SWOT_M_ll; SWOT_M_ll_ti];
    end

    DATA_in_interval = [SWOT_M_ssha + SWOT_M_hret + SWOT_M_rollerror].* ...
                      ~[SWOT_M_ssha_qc + SWOT_M_rollerror_qc]./~[SWOT_M_ssha_qc + SWOT_M_rollerror_qc];
    XY_in_interval = SWOT_M_ll;

    SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
        [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
    SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;

    SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
    SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];

    % % % Average this time-interval's structure function and do a NLLSF
    BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
    BIN_CENTR = interval_avg(SEP_VEC,SEP_VEC,BIN_EDGES);
    BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
        [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
        [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
    % MEAN_SF_ti = interval_avg(SEP_VEC,SF_VEC,BIN_EDGES); % figure;plot(BIN_CENTR,MEAN_SF_ti,'.-')
    MEAN_SF_ti = interval_avg(SEP_MAT_in_interval(:),SF_MAT_in_interval(:),BIN_EDGES); % figure;plot(BIN_CENTR,MEAN_SF_ti,'.-')

    % NNLSF:
    if sum(~isfinite(MEAN_SF_ti)) < 7 % check for enough data
        [NLLSF_COEF,~,~,ConvergenceRecord] = ...
            nonlinear_lsqf(MEAN_SF_ti(2:[end-1]),BIN_CENTR(2:[end-1]),...
            {'a*t.^b','a','b'},...
            [0.00001,1.5],{'t.^b','a*log(t).*t.^b'},ToleranceLevel,   N_nllsf_iterations);
        if ConvergenceRecord(end) < ToleranceLevel % good result (converged)
            Slope_from_SF = [Slope_from_SF ; (-1-NLLSF_COEF(2))];
        else % bad result (didn't converge)
            Slope_from_SF = [Slope_from_SF ; NaN];
        end
    else
        Slope_from_SF = [Slope_from_SF ; NaN];
    end

    Time_for_SF = [Time_for_SF ; mean(SWOT_M_time,'omitnan')];

    if ~mod(ti,10)
        disp([num2str(ti) '/' num2str( length(SWOT.time) )])
    end
end

SF_TIME_SERIES.SWOT.Slope_from_SF = Slope_from_SF;
SF_TIME_SERIES.SWOT.Time = Time_for_SF;

%% PLOT TIME SERIES AND EVOLVING STERIC HEIGHTS

close all

TEST_INTERVAL = [datenum('2023-02-15 00:00:00') datenum('2023-10-15 00:00:00')] - T0; % all times
TEST_INTERVAL = [datenum('2023-04-01 00:00:00') datenum('2023-07-10 00:00:00')] - T0;
% TEST_INTERVAL = [datenum('2023-08-07 00:00:00') datenum('2023-08-08 00:00:00')] - T0;
% TEST_INTERVAL = [datenum('2023-03-01 00:00:00') datenum('2023-12-31 00:00:00')] - T0;

figure('Color',[1 1 1])
% tiledlayout(2,1, 'Padding', 'none', 'TileSpacing', 'compact');
tiledlayout(2,1, 'TileSpacing', 'compact');
nexttile

MM_list_lats = [];
for mi=1:length(MM_list)
    MM = MM_list{mi};
    MM_list_lats(mi) = eval([MM '.P.LATITUDE;']);
end
[~,LatOrder] = sort(MM_list_lats);
LatOrder = flip(LatOrder);
COLOR = turbo(length(MM_list));

% First plot in time all the steric height:
II = 1;
LEG_TEXT = {};
for mi = LatOrder%1:length(MM_list)
    MM = MM_list{mi};
    if strcmp(MM(1),'P')
        % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
        %           interp1_custom( T0 + eval([MM '.F.TIME']), ...
        %                  eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
        %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})'])), ...
        %           '.:', 'MarkerSize',10);
        
        P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ' ... <-- Time interval
                   '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                      '~SPIKE_IND{mi}']); % <-- Spike removal
        plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(   P_IND)']) , ...
                  100*[eval([MM '_sh_600(P_IND)']) ], ...
                  '.-','Color',COLOR(II,:),'HandleVisibility','off'); hold on
        plot(nan,nan,'.','Color',COLOR(II,:),'MarkerSize',50)

        disp(MM)
        disp(var(100*[eval([MM '_sh_600(P_IND)']) ],'omitnan'))
    elseif strcmp(MM(1),'S')
        % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})']), ...
        %      interp1_custom( T0 + eval([MM '.F.TIME']), ...
        %                  eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
        %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ~SPIKE_IND{mi})'])) , ...
        %      '.-')

        P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > TEST_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < TEST_INTERVAL(2) & ' ... <-- Time interval
                    '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                    '~SPIKE_IND{mi}']); % <-- Spike removal
        plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(   P_IND)']), ...
                 100*[eval([MM '_sh_600(P_IND)']) ], ...
                 '.-','Color',COLOR(II,:),'HandleVisibility','off'); hold on
        plot(nan,nan,'.','Color',COLOR(II,:),'MarkerSize',50)

        disp(MM)
        disp(var(100*[eval([MM '_sh_600(P_IND)']) ],'omitnan'))
    else
    end
    LEG_TEXT = {LEG_TEXT{:}, MM};
    II = II + 1;
end
legend(LEG_TEXT,'Orientation','horizontal','Location','northwest')
datetick;set(gca,'XTickLabel','')
% title(['Steric height data from ' datestr(T0 + TEST_INTERVAL(1)) ' to ' datestr(T0 + TEST_INTERVAL(2))])
ylabel(['Steric height anomaly' newline ' from 600 m (cm)'])
set(gca,'XLim',T0 + TEST_INTERVAL([1 end]))
set(gca,'FontSize',30)



disp(['24-Feb-2023 22:49:11 is the first time with 6 records at once'])
disp(['15-Sep-2023 14:51:09 is the last  time with 6 records at once'])

%% Plot evolving slopes

nexttile

plot(SF_TIME_SERIES.Step_24hr.Time, SF_TIME_SERIES.Step_24hr.Slope_from_SF, '.-','MarkerSize',40,'Color',[1 1 1]*0.7); hold on
plot(SF_TIME_SERIES.Step_12hr.Time, SF_TIME_SERIES.Step_12hr.Slope_from_SF, '.-','MarkerSize',30,'Color',[51 153 255]/255)
plot(SF_TIME_SERIES.Step_6hr.Time , SF_TIME_SERIES.Step_6hr.Slope_from_SF,  '.-','MarkerSize',20,'Color',[255 128 0]/255)
plot(SF_TIME_SERIES.Step_1hr.Time , SF_TIME_SERIES.Step_1hr.Slope_from_SF,  '.-','MarkerSize',10,'Color',[76 153 0]/255)
plot(SF_TIME_SERIES.SWOT.Time , SF_TIME_SERIES.SWOT.Slope_from_SF,  'k*','MarkerSize',10,'LineWidth',2)
%$ Calculate the time-evolving spectral slope here, if deemed necessary
legend('24 hr','12 hr','6 hr','1 hr','SWOT (S.F.)', 'SWOT (Welch)')
datetick
ylabel(['Equivalent Spectral Slope'])
set(gca,'XLim',T0 + TEST_INTERVAL([1 end]))
set(gca,'FontSize',30)
% set(gcf,'Position',[-1661         120        1395         835])
set(gcf,'Position',[-1744         120        1682         835])

% %%
% exportgraphics(gcf,...
%     '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/extra_FigA.pdf',...
%     'BackgroundColor','none','ContentType','vector')

%%
%%
%% Save variables as desired

error('Not yet adjusted for truly fixed 600 or 1800 m levels (if that''s even needed here)')

Mooring_list_withP4 =    {'P1','P2','P3',  'P4',   'P5','P6','P7',    'S1','S2','S3','S4'};
Mooring_list_withoutP4 = {'P1','P2','P3',          'P5','P6','P7',    'S1','S2','S3','S4'};
% Open P4 manually:
P4_file = ['/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/PROF/' ...
    'SWOTPOSTLAUNCH_L2_JPLQC_MOORING-P4_CTD-PROFILER-BOTH_START20230219_END20231002_DM_VER01_5s4s_alltimes.nc'];
P4.P.STERIC_HEIGHT_MAX_PROF_DEPTH = ncread(P4_file,'STERIC_HEIGHT_MAX_PROF_DEPTH');
P4.P.STERIC_HEIGHT_MIN_PROF_DEPTH = ncread(P4_file,'STERIC_HEIGHT_MIN_PROF_DEPTH');
P4.P.STERIC_HEIGHT_ANOMALY = ncread(P4_file,'STERIC_HEIGHT_ANOMALY');
P4.P.TIME_STERIC_HEIGHT = ncread(P4_file,'TIME_STERIC_HEIGHT');

P4_gpsfile = ['/Users/kachelein/Documents/JPL/work/07.2023_QC_moorings/PROFILERS/' ...
    'SWOTPOSTLAUNCH_L2_JPLQC_MOORING-P4_CTD-PROFILER_START20230220_END20230728_RT_VER002_5s4s_alltimes.nc'];
P4.GPS.TIME_GPS_SURFACE_BUOY      = ncread(P4_gpsfile,'TIME_GPS_SURFACE_BUOY');
P4.GPS.LATITUDE_GPS_SURFACE_BUOY  = ncread(P4_gpsfile,'LATITUDE_GPS_SURFACE_BUOY');
P4.GPS.LONGITUDE_GPS_SURFACE_BUOY = ncread(P4_gpsfile,'LONGITUDE_GPS_SURFACE_BUOY');


for mi = 1:length(Mooring_list_withP4)
    MM = Mooring_list_withP4{mi};
    eval(['SAVED_VARS.' MM '_PROF_STERIC_HEIGHT_ANOM = '  MM '.P.STERIC_HEIGHT_ANOMALY([' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' ...
                                            MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM);']);
    eval(['SAVED_VARS.' MM '_PROF_TIME               = '  MM '.P.TIME_STERIC_HEIGHT(['    MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' ...
                                            MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM);']);
    eval(['SAVED_VARS.' MM '_PROF_required_total_profile_depth = ' MM(1) 'P_RANGE_LIM;'])

    eval(['SAVED_VARS.' MM '_TIME_GPS_SURFACE_BUOY = ' MM '.GPS.TIME_GPS_SURFACE_BUOY;'])
    eval(['SAVED_VARS.' MM '_LONGITUDE_GPS_SURFACE_BUOY = ' MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY;'])
    eval(['SAVED_VARS.' MM '_LATITUDE_GPS_SURFACE_BUOY = ' MM '.GPS.LATITUDE_GPS_SURFACE_BUOY;'])
end
for mi = 1:length(Mooring_list_withoutP4)
    MM = Mooring_list_withoutP4{mi};
    eval(['SAVED_VARS.' MM '_FIXED_STERIC_HEIGHT_ANOM = ' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*' ...
                                       '[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG];']);
    eval(['SAVED_VARS.' MM '_FIXED_TIME              = '  MM '.F.TIME;'])
    eval(['SAVED_VARS.' MM '_FIXED_NOMINAL_DEPTH     = '  MM '.F.DEPTH_NOMINAL_AUGMENTED;'])
end
SAVED_VARS.PROF_RHO_MEAN = S1.P.RHO_TIMEMEAN;
SAVED_VARS.PROF_RHO_MEAN_DEPTH = S1.P.RHO_TIMEMEAN_DEPTH;
SAVED_VARS.FIXED_RHO_MEAN = S1.F.RHO_TIMEMEAN;
SAVED_VARS.FIXED_RHO_MEAN_DEPTH = S1.F.RHO_TIMEMEAN_DEPTH;

save('/Users/kachelein/Documents/JPL/work/06.2024/JPL_QC_StericHeight.mat','-struct','SAVED_VARS')

%% Gliders:
Glider_list = {'ru32_600','ru32_1000','ru38_600','ru38_1000'};


% for gi = 1:length(Glider_list)
%     P_IND = eval(['[' GG '.STERIC_HEIGHT_MAX_PROF_DEPTH - ' GG '.STERIC_HEIGHT_MIN_PROF_DEPTH] > G_RANGE_LIM' ... <-- Glider depth requirement
%                  ]); % <-- No spike removal for gliders
%     DATA_in_interval = [DATA_in_interval ; ...
%         mean(eval([GG '.STERIC_HEIGHT_ANOMALY(P_IND) + GM_offset']), 'omitnan') ]; % no distinction between P and F for gliders
%     XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
%         mean(eval([GG '.LONGITUDE_PROFILE(P_IND) + 1i*' GG '.LATITUDE_PROFILE(P_IND)']), 'omitnan') ];
% end

% save('/Users/kachelein/Documents/JPL/work/06.2024/JPL_QC_RUTGERS_GLIDERS_StericHeight.mat',...
%      'ru32_600','ru32_1000','ru38_600','ru38_1000')

save('/Users/kachelein/Documents/JPL/work/06.2024/JPL_QC_RU32_600m_StericHeight.mat', '-struct','ru32_600')
save('/Users/kachelein/Documents/JPL/work/06.2024/JPL_QC_RU32_1000m_StericHeight.mat', '-struct','ru32_1000')
save('/Users/kachelein/Documents/JPL/work/06.2024/JPL_QC_RU38_600m_StericHeight.mat', '-struct','ru38_600')
save('/Users/kachelein/Documents/JPL/work/06.2024/JPL_QC_RU38_1000m_StericHeight.mat', '-struct','ru38_1000')


% SP_RANGE_LIM = 477;
% PP_RANGE_LIM = 414;

%%
%%
%%
%% Plot a 4-panel comparison of the hour S moorings at 600 and 1800 m, and SWOT

Mooring_list_S_only = {'S1','S2','S3','S4'};
for ti = 1:length(SWOT.time)
    median_SWOT_time(ti) = median(SWOT.time{ti},'omitnan');
end
for mi = 1:length(Mooring_list_S_only)
    MM = Mooring_list_S_only{mi};
    SWOT_M_ssha = []; SWOT_M_hret = []; SWOT_M_rollerror = []; SWOT_M_ssha_qc = []; SWOT_M_rollerror_qc = [];
    for ti = 1:length(SWOT.time)
        MM_GPS_LON = eval([[MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY']]);
        MM_GPS_LAT = eval([[MM '.GPS.LATITUDE_GPS_SURFACE_BUOY']]);
        Mooring_Loc_Mat = abs(SWOT.lon{ti} + 1i*SWOT.lat{ti} - [...
            MM_GPS_LON(dsearchn(eval([MM '.GPS.TIME_GPS_SURFACE_BUOY']) + T0, median(SWOT.time{ti}(:)))) + ...
            MM_GPS_LAT(dsearchn(eval([MM '.GPS.TIME_GPS_SURFACE_BUOY']) + T0, median(SWOT.time{ti}(:))))*1i] );
        % ^ the use of median(SWOT.time{ti}(:)) will mean that the true
        % time will be off by a maximum of 3-6 minutes, which should
        % be ok considering the time steps of profilers is 0.5-4 hr

        SWOT_M_ssha         = [SWOT_M_ssha         ; SWOT.ssha_karin_2{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) )];
        SWOT_M_hret         = [SWOT_M_hret         ; SWOT.internal_tide_hret{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) )];
        SWOT_M_rollerror    = [SWOT_M_rollerror    ; SWOT.height_cor_xover{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) )];
        SWOT_M_ssha_qc      = [SWOT_M_ssha_qc      ; SWOT.ssha_karin_2_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) )];
        SWOT_M_rollerror_qc = [SWOT_M_rollerror_qc ; SWOT.height_cor_xover_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) )];

    end
    eval([MM '_locational_SWOT_ssha =  [SWOT_M_ssha + SWOT_M_hret + SWOT_M_rollerror].*' ...
        '~[SWOT_M_ssha_qc + SWOT_M_rollerror_qc]./~[SWOT_M_ssha_qc + SWOT_M_rollerror_qc];'])
end
%%
close all
figure('Color','w')
tiledlayout(2,2)%, 'TileSpacing', 'compact');
XLIM_4panel = datenum(['2023-03-01';'2023-10-01']);
for mi = 1:length(Mooring_list_S_only)
    MM = Mooring_list_S_only{mi};

    IND = eval([ '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
                   '~SPIKE_IND{mi+6}']); % <-- Spike removal

    % subplot(2,2,mi)
    nexttile
    plot(eval([MM '.P.TIME_MIDPROFILE(IND)']) + T0,...
         eval([MM '_sh_600(IND)']),'.-'); hold on
    plot(eval([MM '.P.TIME_MIDPROFILE(IND)']) + T0,...
         eval([MM '_sh_1800(IND)']),'.-');

    plot(median_SWOT_time, eval([MM '_locational_SWOT_ssha']), '.-','Color',[1 1 1]*0)
    datetick
    set(gca,'FontSize',16,'YLim',[-0.12 0.2])
    set(gca,'XLim',XLIM_4panel)
    ylabel('Steric Height/SSHA')
    title(MM)
end

%%
exportgraphics(gcf,...
    '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/extra_FigB.pdf',...
    'BackgroundColor','none','ContentType','vector')

%%
%%
%% Regular, no-frills variance

T1 = datenum('2023-04-01 00:00:00');
T2 = datenum('2023-07-10 23:00:00');

disp('variance at 600 for P1, P2, P3, P5, P6, P7,   S1, S2, S3, S4')
... convert to cm * Steric Height @ data that integrate enough depth & despiked & within the cal/val date range
VAR_VEC_600 = nan(length(MM_list),1);
LEN_VEC_600 = nan(length(MM_list),1);
for ii = 1:length(MM_list)
    MM = MM_list{ii};
    VAR_VEC_600(ii) = eval([ ' var( 100*' MM '_sh_600([' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii} & ' ...
                                   '[T0+' MM '.P.TIME_STERIC_HEIGHT] > T1 & [T0+' MM '.P.TIME_STERIC_HEIGHT] < T2),''omitnan'' ) '   ]);
    LEN_VEC_600(ii) = eval([ ' length( 100*' MM '_sh_600([' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii} & ' ...
                                     '[T0+' MM '.P.TIME_STERIC_HEIGHT] > T1 & [T0+' MM '.P.TIME_STERIC_HEIGHT] < T2)) '   ]);
end
[reshape([MM_list{:}],[2,numel([MM_list{:}])/2])' , repmat('    ',[10,1]), num2str(VAR_VEC_600)]


disp('variance at 1800 for S1, S2, S3, S4')


VAR_VEC_1800 = nan(length(MM_list),1);
LEN_VEC_1800 = nan(length(MM_list),1);
for ii = 7:length(MM_list)
    MM = MM_list{ii};
    VAR_VEC_1800(ii) = eval([ ' var( 100*' MM '_sh_1800([' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii} & ' ...
                                   '[T0+' MM '.P.TIME_STERIC_HEIGHT] > T1 & [T0+' MM '.P.TIME_STERIC_HEIGHT] < T2),''omitnan'' ) '   ]);
    LEN_VEC_1800(ii) = eval([ ' length( 100*' MM '_sh_1800([' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii} & ' ...
                                      '[T0+' MM '.P.TIME_STERIC_HEIGHT] > T1 & [T0+' MM '.P.TIME_STERIC_HEIGHT] < T2)) '   ]);
end
[reshape([MM_list{:}],[2,numel([MM_list{:}])/2])' , repmat('    ',[10,1]), num2str(VAR_VEC_1800)]

close all
figure
subplot(211)
for ii = 1:length(MM_list)
    MM = MM_list{ii};
    T_temp = eval([       MM '.P.TIME_STERIC_HEIGHT([' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii})']);
    H_temp = eval(['100*' MM '_sh_600(['               MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii})']);
    plot( T0 + T_temp , H_temp , '.-'); hold on
end
datetick
xlim([datenum('2023-02-01 00:00:00') datenum('2023-11-01 00:00:00')])
subplot(212)
for ii = 1:length(MM_list)
    MM = MM_list{ii};
    T_temp = eval([       MM '.P.TIME_STERIC_HEIGHT([' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii} & ' ...
                                     '[T0+' MM '.P.TIME_STERIC_HEIGHT] > T1 & [T0+' MM '.P.TIME_STERIC_HEIGHT] < T2)']);
    H_temp = eval(['100*' MM '_sh_600(['               MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii} & ' ...
                                     '[T0+' MM '.P.TIME_STERIC_HEIGHT] > T1 & [T0+' MM '.P.TIME_STERIC_HEIGHT] < T2)']);
    plot( T0 + T_temp , H_temp , '.-'); hold on
end
datetick
xlim([datenum('2023-02-01 00:00:00') datenum('2023-11-01 00:00:00')])

figure
subplot(211)
for ii = 7:length(MM_list)
    MM = MM_list{ii};
    T_temp = eval([       MM '.P.TIME_STERIC_HEIGHT([' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii})']);
    H_temp = eval(['100*' MM '_sh_1800(['               MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii})']);
    plot( T0 + T_temp , H_temp , '.-'); hold on
end
datetick
xlim([datenum('2023-02-01 00:00:00') datenum('2023-11-01 00:00:00')])
subplot(212)
for ii = 7:length(MM_list)
    MM = MM_list{ii};
    T_temp = eval([       MM '.P.TIME_STERIC_HEIGHT([' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii} & ' ...
                                     '[T0+' MM '.P.TIME_STERIC_HEIGHT] > T1 & [T0+' MM '.P.TIME_STERIC_HEIGHT] < T2)']);
    H_temp = eval(['100*' MM '_sh_1800(['               MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ~SPIKE_IND{ii} & ' ...
                                     '[T0+' MM '.P.TIME_STERIC_HEIGHT] > T1 & [T0+' MM '.P.TIME_STERIC_HEIGHT] < T2)']);
    plot( T0 + T_temp , H_temp , '.-'); hold on
end
datetick
xlim([datenum('2023-02-01 00:00:00') datenum('2023-11-01 00:00:00')])


%% Auxiliary function

function Vq = interp1IA_function(X,V,Xq,dX)
    Vq = nan(size(Xq));
    for ii = 1:length(Xq)
        Vq(ii) = mean(V(X > [Xq(ii) - dX/2] & X <= [Xq(ii) + dX/2]),'omitnan');
    end
end

function V_avg = interval_avg(X,V,X_int)
    % X_int = edges of bins
    % X = location of data
    % V = data
    % V_avg = interval-averaged data
    V_avg = nan(size(X_int(2:end)));
    for ii = 1:[length(X_int) - 1]
        V_avg(ii) = mean(V(X > X_int(ii) & X < X_int(ii+1)),'omitnan');
    end
end

function V_stdofmean = interval_stdofmean(X,V,X_int)
    % X_int = edges of bins
    % X = location of data
    % V = data
    % V_stdofmean = std of mean = std/sqrt(N)
    V_stdofmean = nan(size(X_int(2:end)));
    for ii = 1:[length(X_int) - 1]
        V_stdofmean(ii) = std(V(X > X_int(ii) & X < X_int(ii+1)),'omitnan')/...
            sqrt(sum(isfinite(V(X > X_int(ii) & X < X_int(ii+1)))));
    end
end

function mean_SF_out = errorbar_stdofmean(bin_center, sep_vec, sf_vec, bin_edges, varargin)
    % Custom plotting function so that the plotting window doesn't get co
    % cluttered.
    if nargin == 4
        errorbar(bin_center, ...
            interval_avg(sep_vec,sf_vec,bin_edges), ...
            interval_stdofmean(sep_vec,sf_vec,bin_edges), ...
            interval_stdofmean(sep_vec,sf_vec,bin_edges), ...
            '.-', 'MarkerSize',25,'LineWidth',2)
    elseif nargin == 5
        errorbar(bin_center, ...
            interval_avg(sep_vec,sf_vec,bin_edges), ...
            interval_stdofmean(sep_vec,sf_vec,bin_edges), ...
            interval_stdofmean(sep_vec,sf_vec,bin_edges), ...
            '.-', 'MarkerSize',25,'LineWidth',2,'Color',varargin{1})
    elseif nargin == 6
        errorbar(bin_center, ...
            interval_avg(sep_vec,sf_vec,bin_edges), ...
            interval_stdofmean(sep_vec,sf_vec,bin_edges), ...
            interval_stdofmean(sep_vec,sf_vec,bin_edges), ...
            varargin{2}, 'MarkerSize',25,'LineWidth',1,'Color',varargin{1})
    else
        error('4, 5, or 6 inputs only')
    end
    mean_SF_out = interval_avg(sep_vec,sf_vec,bin_edges);
end

function V_avg = interval_median(X,V,X_int)
    % X_int = edges of bins
    % X = location of data
    % V = data
    % V_avg = interval-averaged data
    V_avg = nan(size(X_int(2:end)));
    for ii = 1:[length(X_int) - 1]
        V_avg(ii) = median(V(X > X_int(ii) & X < X_int(ii+1)),'omitnan');
    end
end






%% SWOT METADATA

% % % ncdisp('/mnt/flow/swot/KaRIn/SWOT_L2_LR_SSH_2.0_AllToDate/SWOT_L2_LR_SSH_Basic_498_013_20230422T043022_20230422T052129_PGC0_02.nc')
    % time                                 
    %        Size:       9866x1
    %        Dimensions: num_lines
    %        Datatype:   double
    %        Attributes:
    %                    _FillValue         = 9.969209968386869e+36
    %                    long_name          = 'time in UTC'
    %                    standard_name      = 'time'
    %                    calendar           = 'gregorian'
    %                    tai_utc_difference = 37
    %                    leap_second        = '0000-00-00T00:00:00Z'
    %                    units              = 'seconds since 2000-01-01 00:00:00.0'
    %                    comment            = 'Time of measurement in seconds in the UTC time scale since 1 Jan 2000 00:00:00 UTC. [tai_utc_difference] is the difference between TAI and UTC reference time (seconds) for the first measurement of the data set. If a leap second occurs within the data set, the attribute leap_second is set to the UTC time at which the leap second occurs.'
    %   latitude
    %        Size:       69x9866
    %        Dimensions: num_pixels,num_lines
    %        Datatype:   int32
    %        Attributes:
    %                    _FillValue    = 2147483647
    %                    long_name     = 'latitude (positive N, negative S)'
    %                    standard_name = 'latitude'
    %                    units         = 'degrees_north'
    %                    scale_factor  = 1e-06
    %                    valid_min     = -80000000
    %                    valid_max     = 80000000
    %                    comment       = 'Latitude of measurement [-80,80]. Positive latitude is North latitude, negative latitude is South latitude.'
    % longitude !!!!! We modified this to be +-180
    %        Size:       69x9866
    %        Dimensions: num_pixels,num_lines
    %        Datatype:   int32
    %        Attributes:
    %                    _FillValue    = 2147483647
    %                    long_name     = 'longitude (degrees East)'
    %                    standard_name = 'longitude'
    %                    units         = 'degrees_east'
    %                    scale_factor  = 1e-06
    %                    valid_min     = 0
    %                    valid_max     = 359999999
    %                    comment       = 'Longitude of measurement. East longitude relative to Greenwich meridian.'
    % ssha_karin_2                         
    %        Size:       69x9866
    %        Dimensions: num_pixels,num_lines
    %        Datatype:   int32
    %        Attributes:
    %                    _FillValue   = 2147483647
    %                    long_name    = 'sea surface height anomaly'
    %                    units        = 'm'
    %                    scale_factor = 0.0001
    %                    quality_flag = 'ssha_karin_2_qual'
    %                    valid_min    = -1000000
    %                    valid_max    = 1000000
    %                    coordinates  = 'longitude latitude'
    %                    comment      = UNSUPPORTED DATATYPE
    % ssha_karin_2_qual                    
    %        Size:       69x9866
    %        Dimensions: num_pixels,num_lines
    %        Datatype:   uint32
    %        Attributes:
    %                    _FillValue    = 4294967295
    %                    long_name     = 'sea surface height anomaly quality flag'
    %                    standard_name = 'status_flag'
    %                    flag_meanings = 'suspect_large_ssh_delta suspect_large_ssh_std suspect_large_ssh_window_std suspect_beam_used suspect_less_than_nine_beams suspect_ssb_out_of_range suspect_pixel_used suspect_num_pt_avg suspect_karin_telem suspect_orbit_control suspect_sc_event_flag suspect_tvp_qual suspect_volumetric_corr degraded_ssb_not_computable degraded_media_delays_missing degraded_beam_used degraded_large_attitude degraded_karin_ifft_overflow bad_karin_telem bad_very_large_attitude bad_tide_corrections_missing bad_outside_of_range degraded bad_not_usable'
    %                    flag_masks    = [1           2           4           8          16          64         128         256         512        1024        2048        4096        8192       32768       65536      131072      262144      524288    16777216    33554432    67108864   536870912  1073741824  2147483648]
    %                    valid_min     = 0
    %                    valid_max     = 3876569055
    %                    coordinates   = 'longitude latitude'
    %                    comment       = 'Quality flag for the SSHA from KaRIn in the ssha_karin_2 variable'
    % internal_tide_hret                   
    %        Size:       69x9866
    %        Dimensions: num_pixels,num_lines
    %        Datatype:   int16
    %        Attributes:
    %                    _FillValue   = 32767
    %                    long_name    = 'coherent internal tide (HRET)'
    %                    source       = 'Zaron (2019)'
    %                    units        = 'm'
    %                    scale_factor = 0.0001
    %                    valid_min    = -2000
    %                    valid_max    = 2000
    %                    coordinates  = 'longitude latitude'
    %                    comment      = 'Coherent internal ocean tide. This value is subtracted from the ssh_karin and ssh_karin_2 to compute ssha_karin and ssha_karin_2, respectively.'
    % height_cor_xover                     
    %        Size:       69x9866
    %        Dimensions: num_pixels,num_lines
    %        Datatype:   int32
    %        Attributes:
    %                    _FillValue   = 2147483647
    %                    long_name    = 'height correction from crossover calibration'
    %                    units        = 'm'
    %                    scale_factor = 0.0001
    %                    quality_flag = 'height_cor_xover_qual'
    %                    valid_min    = -100000
    %                    valid_max    = 100000
    %                    coordinates  = 'longitude latitude'
    %                    comment      = 'Height correction from crossover calibration. To apply this correction the value of height_cor_xover should be added to the value of ssh_karin, ssh_karin_2, ssha_karin, and ssha_karin_2.'
    % height_cor_xover_qual                
    %        Size:       69x9866
    %        Dimensions: num_pixels,num_lines
    %        Datatype:   uint8
    %        Attributes:
    %                    _FillValue    = 255
    %                    long_name     = 'quality flag for height correction from crossover calibration'
    %                    standard_name = 'status_flag'
    %                    flag_meanings = 'good suspect bad'
    %                    flag_values   = [0  1  2]
    %                    valid_min     = 0
    %                    valid_max     = 2
    %                    coordinates   = 'longitude latitude'
    %                    comment       = 'Flag indicating the quality of the height correction from crossover calibration.  Values of 0, 1, and 2 indicate that the correction is good, suspect, and bad, respectively.'

%% Discarded material

% SWOT_M_time = [];
% SWOT_M_ssha = [];
% SWOT_M_ssha_qc = [];
% SWOT_M_rollerror = [];
% SWOT_M_rollerror_qc = [];
% SWOT_M_hret = [];
% SWOT_M_ll = [];
% tic
% for ti = 1:length(SWOT.time)
%     for mi = 1:length(MM_list)
%         MM = MM_list{mi};
% 
%         MM_GPS_LON = eval([[MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY']]);
%         MM_GPS_LAT = eval([[MM '.GPS.LATITUDE_GPS_SURFACE_BUOY']]);
%         Mooring_Loc_Mat = abs(SWOT.lon{ti} + 1i*SWOT.lat{ti} - [...
%             MM_GPS_LON(dsearchn(eval([MM '.GPS.TIME_GPS_SURFACE_BUOY']) + T0, median(SWOT.time{ti}(:)))) + ...
%             MM_GPS_LAT(dsearchn(eval([MM '.GPS.TIME_GPS_SURFACE_BUOY']) + T0, median(SWOT.time{ti}(:))))*1i] );
%             % ^ the use of median(SWOT.time{ti}(:)) will mean that the true
%             % time will be off by a maximum of 3-6 minutes, which should
%             % be ok considering the time steps of profilers is 0.5-4 hr
%         SWOT_M_time_ti = SWOT.time{ti}( ~~sum(Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)),1)' );
%         SWOT_M_time = [SWOT_M_time ; SWOT_M_time_ti];
% 
%         SWOT_M_ssha_ti = SWOT.ssha_karin_2{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
%         SWOT_M_ssha = [SWOT_M_ssha ; SWOT_M_ssha_ti];
% 
%         SWOT_M_ssha_qc_ti = SWOT.ssha_karin_2_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
%         SWOT_M_ssha_qc = [SWOT_M_ssha_qc ; SWOT_M_ssha_qc_ti];
% 
%         SWOT_M_rollerror_ti = SWOT.height_cor_xover{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
%         SWOT_M_rollerror = [SWOT_M_rollerror ; SWOT_M_rollerror_ti];
% 
%         SWOT_M_rollerror_qc_ti = SWOT.height_cor_xover_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
%         SWOT_M_rollerror_qc = [SWOT_M_rollerror_qc ; SWOT_M_rollerror_qc_ti];
% 
%         SWOT_M_hret_ti = SWOT.internal_tide_hret{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
%         SWOT_M_hret = [SWOT_M_hret ; SWOT_M_hret_ti];
% 
%         SWOT_M_ll_ti = SWOT.lon{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) ) + ...
%                        SWOT.lat{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) )*1i;
%         SWOT_M_ll = [SWOT_M_ll; SWOT_M_ll_ti];
% 
%     end
% end
% toc





% function V_25 = interval_25prctile(X,V,X_int)
%     V_25 = nan(size(X_int(2:end)));
%     for ii = 1:[length(X_int) - 1]
%         V_25(ii) = prctile(V(X > X_int(ii) & X < X_int(ii+1)),25);
%         % V_25(ii) = mean(V(X > X_int(ii) & X < X_int(ii+1)),'omitnan') - ...
%         %     std(V(X > X_int(ii) & X < X_int(ii+1)),'omitnan'); % for -1std
%     end
% end
% function V_75 = interval_75prctile(X,V,X_int)
%     V_75 = nan(size(X_int(2:end)));
%     for ii = 1:[length(X_int) - 1]
%         V_75(ii) = prctile(V(X > X_int(ii) & X < X_int(ii+1)),75);
%         % V_75(ii) = mean(V(X > X_int(ii) & X < X_int(ii+1)),'omitnan') + ...
%         %     std(V(X > X_int(ii) & X < X_int(ii+1)),'omitnan'); % for +1std
%     end
% end
% function errorbar_25mean75(bin_center, sep_vec, sf_vec, bin_edges)
%     % Custom plotting function so that the plotting window doesn't get co
%     % cluttered.
%     errorbar(bin_center, ...
%              interval_avg(sep_vec,sf_vec,bin_edges), ...
%              interval_avg(sep_vec,sf_vec,bin_edges) - interval_25prctile(sep_vec,sf_vec,bin_edges), ...
%              interval_75prctile(sep_vec,sf_vec,bin_edges) - interval_avg(sep_vec,sf_vec,bin_edges), ...
%              '.-', 'MarkerSize',25,'LineWidth',2)
% end