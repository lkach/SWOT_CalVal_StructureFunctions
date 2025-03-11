%% Figure for Manuscript - spectral analysis of mooring data
%% Time-series analysis of mooring data to ascertain the contribution of internal waves to steric height
%% Load MOORINGS

MP_DIR = dir( '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/PROF/*.nc');
MF_DIR = dir('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/FIXED/*.nc');
% % % Separate files for GPS mooring data:
M_GPS_DIR = dir('/Users/kachelein/Documents/JPL/work/07.2023_QC_moorings/PROFILERS/SWOT*.nc');
warning('Before publication, be sure to package the GPS data with the other data (or ideally SIO will send us an updated version).')

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
    if contains(M_GPS_DIR(mi).name,'-P4_')
        UNUSABLE_MOORING_IND_GPS = [UNUSABLE_MOORING_IND_GPS; mi];
    else
    end
end
MP_DIR    =    MP_DIR(setdiff(1:length(MP_DIR)   , UNUSABLE_MOORING_IND_P  ));
M_GPS_DIR = M_GPS_DIR(setdiff(1:length(M_GPS_DIR), UNUSABLE_MOORING_IND_GPS));


EXCLUDED_VARS_MP = {'DEPTH','TIME','RHO',   'CNDC','PRES','PSAL','TEMP','PROF_NUM'};
EXCLUDED_VARS_MF = {'DEPTH',       'RHO',   'CNDC','CNDC_QC','PRES','PRES_QC','PSAL','PSAL_QC','TEMP','TEMP_QC'};

for mi = 1:length(MM_list)
    MM = MM_list{mi};
    % PROF
    eval([MM '.P = ncreadall(''' MP_DIR(mi).folder '/' MP_DIR(mi).name ''' , EXCLUDED_VARS_MP);']);
    % FIXED
    eval([MM '.F = ncreadall(''' MF_DIR(mi).folder '/' MF_DIR(mi).name ''' , EXCLUDED_VARS_MF);']);
    % GPS locaitons (separately hadled):
    eval([MM '.GPS.TIME_GPS_SURFACE_BUOY = ncread(''' M_GPS_DIR(mi).folder '/' M_GPS_DIR(mi).name ''',''TIME_GPS_SURFACE_BUOY'');']);
    eval([MM '.GPS.LATITUDE_GPS_SURFACE_BUOY = ncread(''' M_GPS_DIR(mi).folder '/' M_GPS_DIR(mi).name ''',''LATITUDE_GPS_SURFACE_BUOY'');']);
    eval([MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY = ncread(''' M_GPS_DIR(mi).folder '/' M_GPS_DIR(mi).name ''',''LONGITUDE_GPS_SURFACE_BUOY'');']);
end
% Workspace variables take up 283100496 bytes of memory = 283.1005 MB of memory.

%% Correct the GPS locations by taking info from .GPS and putting it in .P

interp1nan = @(X,V,Xq) interp1(X(isfinite(X)&isfinite(V)),V(isfinite(X)&isfinite(V)),Xq);
for mi = 1:length(MM_list)
    MM = MM_list{mi};
        eval([replace('#.P.LATITUDE_GPSSB_STERIC_HEIGHT','#',MM) '=' ...
              replace('interp1nan(#.GPS.TIME_GPS_SURFACE_BUOY, #.GPS.LATITUDE_GPS_SURFACE_BUOY, #.P.TIME_STERIC_HEIGHT);','#',MM)]);
        eval([replace('#.P.LONGITUDE_GPSSB_STERIC_HEIGHT','#',MM) '=' ...
              replace('interp1nan(#.GPS.TIME_GPS_SURFACE_BUOY, #.GPS.LONGITUDE_GPS_SURFACE_BUOY, #.P.TIME_STERIC_HEIGHT);','#',MM)]);
end

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



SP_RANGE_LIM = 477;
PP_RANGE_LIM = 414;

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
    plot(TT, diff([TS;NaN]), '.-'); hold on
end
figure(1)
datetick

% % % For testing purposes, reset SPIKE_IND{...} to all false:
% for mi = 1:length(MM_list)
%     SPIKE_IND{mi} = 0*SPIKE_IND{mi};
% end

%% Visualize the data that will be analyzed in the frequency domain:

close all

SF_TIME_STEP = 1; % days
% TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];
SF_TIME_INTERVAL = [TIME_BIN_EDGES(1), TIME_BIN_EDGES(end)]';

% First, all data (profilers only) together to ensure that nothing is strange:
figure('Color',[1 1 1])
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    IND = [MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM & ' ...
                    '~SPIKE_IND{mi} & [' MM '.P.TIME_STERIC_HEIGHT + T0] > SF_TIME_INTERVAL(1) & [' MM '.P.TIME_STERIC_HEIGHT + T0] < SF_TIME_INTERVAL(2)'];
    plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']) , ...
              eval([MM '.P.STERIC_HEIGHT_ANOMALY(' IND ')']) , ...
              '.-'); hold on
end
datetick
title('Profilers only')

% Now plot the (truly) fixed bottom steric heights
figure('Color',[1 1 1])
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    IND = [MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM & ' ...
                    '~SPIKE_IND{mi} & [' MM '.P.TIME_STERIC_HEIGHT + T0] > SF_TIME_INTERVAL(1) & [' MM '.P.TIME_STERIC_HEIGHT + T0] < SF_TIME_INTERVAL(2)'];
    plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']) , ...
              eval([MM '_sh_600(' IND ')']) , ...
              '.-'); hold on
end
datetick
title('Down to 600 m (truly fixed)')

% Now plot the (truly) fixed bottom steric heights
figure('Color',[1 1 1])
for mi = (length(MM_list)-3):length(MM_list)
    MM = MM_list{mi};
    IND = [MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM & ' ...
                    '~SPIKE_IND{mi} & [' MM '.P.TIME_STERIC_HEIGHT + T0] > SF_TIME_INTERVAL(1) & [' MM '.P.TIME_STERIC_HEIGHT + T0] < SF_TIME_INTERVAL(2)'];
    plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']) , ...
              eval([MM '_sh_1800(' IND ')']) , ...
              '.-'); hold on
end
datetick
title('Down to 1800 m (truly fixed)')

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

%%
%%
%% Plot the segments for validation of the method:

close all

% All spectra together
VAR_VEC = nan(length(MM_list),1);
SPEC_SUM_VEC = nan(length(MM_list),1);
SPEC_SUM_VEC_SUPERDIURNAL = nan(length(MM_list),1);
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    IND = [MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM & ' ...
                    '~SPIKE_IND{mi} & [' MM '.P.TIME_STERIC_HEIGHT + T0] > SF_TIME_INTERVAL(1) & [' MM '.P.TIME_STERIC_HEIGHT + T0] < SF_TIME_INTERVAL(2)'];
    % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']) , ...
    %           eval([MM '.P.STERIC_HEIGHT_ANOMALY(' IND ')']) , ...
    %           '.-');
    if [median(diff(eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')'])))\1] < 36
        OBS_PER_DAY = 24;
    elseif [median(diff(eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')'])))\1] > 36
        OBS_PER_DAY = 48;
    end
    [SPEC,FF,ERR] = ...
    nunanspectrum(eval([MM '_sh_600(' IND ')']), ...
                  eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']),'day',...
                  'Segments',3,...
                  'Window','hanning',     'Freq',[(1/30):(1/30):(OBS_PER_DAY/2)],...
                  'Plot',false,'PlotSegments',true); hold on
    VAR_VEC(mi) = var(eval([MM '_sh_600(' IND ')']),'omitnan');
    SPEC_SUM_VEC(mi) = sum(SPEC*FF(1));
    SPEC_SUM_VEC_SUPERDIURNAL(mi) = sum(SPEC(FF>=[24/24])*FF(1));
end
%%
%%
%% 600 m - TIME SERIES ANALYSIS, NUFFT AND LOMB-SCARGLE VIA MATLAB'S plomb

close all

TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];
SF_TIME_INTERVAL = [TIME_BIN_EDGES(1), TIME_BIN_EDGES(end)]';
SEGS = 3;

MM_list_char = [];
for ii = 1:length(MM_list)
    MM_list_char(ii,:) = [MM_list{ii} '  '];
end
MM_list_char = char(MM_list_char);

MM_list_lats = [];
for mi=1:length(MM_list)
    MM = MM_list{mi};
    MM_list_lats(mi) = eval([MM '.P.LATITUDE;']);
end
[~,LatOrder] = sort(MM_list_lats);
LatOrder = flip(LatOrder)';
COLOR = turbo(length(MM_list));

% All spectra together
figure('Color',[1 1 1])
figure('Color',[1 1 1])
VAR_VEC = nan(length(MM_list),1);

% For nufft method:
FF = [(1/30):(1/30):(12)]';
SPEC_SUM_VEC = nan(length(MM_list),1);
SPEC_SUM_VEC_SUPERDIURNAL = nan(length(MM_list),1);
SPEC_SUM_VEC_M2 = nan(length(MM_list),1);
SPEC_MEAN = nan(length(FF),length(MM_list));

% For plomb method:
FF_LS = FF;%[(1/90):(1/90):(12)]';
SPEC_SUM_VEC_LS = nan(length(MM_list),1);
SPEC_SUM_VEC_LS_SUPERDIURNAL = nan(length(MM_list),1);
SPEC_SUM_VEC_LS_M2 = nan(length(MM_list),1);
SPEC_LS_MEAN = nan(length(FF_LS),length(MM_list));

% For estimating uncertainty:
N_unc_dist = 10000;
Unc_dist_nufft_600m = nan(N_unc_dist,length(MM_list));
Unc_dist_plomb_600m = nan(N_unc_dist,length(MM_list));

for mi = 1:length(MM_list)
    MM = MM_list{mi};
    IND = [MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM & ' ...
                    '~SPIKE_IND{mi} & [' MM '.P.TIME_STERIC_HEIGHT + T0] > SF_TIME_INTERVAL(1) & [' MM '.P.TIME_STERIC_HEIGHT + T0] < SF_TIME_INTERVAL(2)'];

    % [SPEC_LS,~] = plomb(eval([MM '.P.STERIC_HEIGHT_ANOMALY(' IND ')']),...
    %                     eval([MM '.P.TIME_STERIC_HEIGHT('    IND ')']), FF_LS);
    % [SPEC_LS,~] = plomb(detrend(eval([MM '.P.STERIC_HEIGHT_ANOMALY(' IND ')'])),...
    %                             eval([MM '.P.TIME_STERIC_HEIGHT('    IND ')']), FF_LS);
    [SPEC_LS,~,ERR] = ...
    nunanspectrum(eval([MM '_sh_600(' IND ')']), ...
                  eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']),'day',...
                  'Segments',SEGS,...
                  'Window','hanning',     'Freq',FF_LS,...
                  'Plot',false,'PlotSegments',false,'Method','plomb'); hold on


    [SPEC,~,ERR] = ...
    nunanspectrum(eval([MM '_sh_600(' IND ')']), ...
                  eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']),'day',...
                  'Segments',SEGS,...
                  'Window','hanning',     'Freq',FF,...
                  'Plot',false,'PlotSegments',false,'Method','nufft'); hold on
    % VAR_VEC(mi) = var(detrend( eval([MM '.P.STERIC_HEIGHT_ANOMALY(' IND ')']) ) ,'omitnan');
    VAR_VEC(mi) = var( ( eval([MM '_sh_600(' IND ')']) ) ,'omitnan');
    figure(1)
        loglog(FF_LS,SPEC_LS,'-','Color',COLOR(dsearchn(LatOrder,mi),:),'LineWidth',1);hold on
    figure(2)
        loglog(FF,SPEC,'-','Color',COLOR(dsearchn(LatOrder,mi),:),'LineWidth',1);hold on

    SPEC_SUM_VEC(mi) = sum(SPEC*FF(1));
    SPEC_SUM_VEC_SUPERDIURNAL(mi) = sum(SPEC(FF>=[24/24])*FF(1));
        Unc_dist_nufft_600m(:,mi) = spec_sum_distribution(SPEC(FF>=[24/24])*FF(1), 2*SEGS, N_unc_dist);
    SPEC_SUM_VEC_M2(mi) = sum(SPEC(dsearchn(FF,[24/12.42]))*FF(1));
    SPEC_MEAN(:,mi) = SPEC;

    SPEC_SUM_VEC_LS(mi) = sum(SPEC_LS*FF_LS(1));
    SPEC_SUM_VEC_LS_SUPERDIURNAL(mi) = sum(SPEC_LS(FF_LS>=[24/24])*FF_LS(1));
        Unc_dist_plomb_600m(:,mi) = spec_sum_distribution(SPEC_LS(FF_LS>=[24/24])*FF_LS(1), 2*SEGS, N_unc_dist);
    SPEC_SUM_VEC_LS_M2(mi) = sum(SPEC_LS(dsearchn(FF_LS,[24/12.42]))*FF_LS(1));
    SPEC_LS_MEAN(:,mi) = SPEC_LS;
end
figure(1)
ylabel('Power spectral density (plomb)')
xlabel('Frequency (cpd)')
set(gca,'XScale','Log','YScale','Log','FontSize',18)
legend(MM_list)

figure(2)
ylabel('Power spectral density (nunanspectrum)')
xlabel('Frequency (cpd)')
set(gca,'XScale','Log','YScale','Log','FontSize',18)
legend(MM_list)

figure('Color',[1 1 1])
loglog(FF,mean(SPEC_MEAN,2,'omitnan'),'.-'); hold on%,'Color',[.2 .5 .8]
loglog(FF_LS,mean(SPEC_LS_MEAN,2,'omitnan'),'.-')%,'Color',[.8 .5 .2]
ylabel('Power spectral density (nunanspectrum)')
xlabel('Frequency (cpd)')
legend('<S(f)> = Welch with nufft','<S(f)> = Lomb-Scargle')

SPEC_MEAN_600 = SPEC_MEAN;
SPEC_LS_MEAN_600 = SPEC_LS_MEAN;


disp('STD in cm (NUNANSPECTRUM NUFFT):')
CM_str = repmat(' cm',length(SPEC_SUM_VEC_SUPERDIURNAL),1);
[MM_list_char num2str(sqrt([SPEC_SUM_VEC_SUPERDIURNAL])*100,'%.2f') CM_str]

disp('%%%%%%%%%%% SUPERDIURNAL %%%%%%%%%%%')
disp('STD in cm (PLOMB):')
CM_str = repmat(' cm',length(SPEC_SUM_VEC_LS_SUPERDIURNAL),1);
[MM_list_char num2str(sqrt([SPEC_SUM_VEC_LS_SUPERDIURNAL])*100,'%.2f') CM_str]
%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Variance in cm^2 (NUNANSPECTRUM NUFFT):')
CM_str = repmat(' cm^2, STD = ',length(SPEC_SUM_VEC_SUPERDIURNAL),1);
[MM_list_char num2str(SPEC_SUM_VEC_SUPERDIURNAL*10000,'%.2f') CM_str num2str(std(Unc_dist_nufft_600m,1)'*10000,'%.2f')]

disp('%%%%%%%%%%% SUPERDIURNAL %%%%%%%%%%%')
disp('Variance in cm^2 (PLOMB):')
CM_str = repmat(' cm^2, STD = ',length(SPEC_SUM_VEC_SUPERDIURNAL),1);
[MM_list_char num2str(SPEC_SUM_VEC_LS_SUPERDIURNAL*10000,'%.2f') CM_str num2str(std(Unc_dist_plomb_600m,1)'*10000,'%.2f')]


%% Variability in different frequency bands:

disp('%%%%%%%%%%%%%%%%%%%%%%% STDEV: UNITS OF CENTIMETERS %%%%%%%%%%%%%%%%%%%%%%%')
disp('diurnal and faster:')
sqrt(sum(mean(SPEC_LS_MEAN(dsearchn(FF_LS,1/1):end,:),2,'omitnan')*FF_LS(1)))*100
sqrt(sum(mean(SPEC_MEAN(   dsearchn(FF   ,1/1):end,:),2,'omitnan')*FF(   1)))*100

disp('M2 +- a few frequencies:')
% sum(mean(SPEC_LS_MEAN(dsearchn(FF_LS,24/12.42) + [-9:9],:),2,'omitnan')*FF_LS(1))
% sum(mean(SPEC_MEAN(   dsearchn(FF   ,24/12.42) + [-3:3],:),2,'omitnan')*FF(   1))
sqrt(sum(mean(SPEC_LS_MEAN(dsearchn(FF_LS,24/12.42 - 3/30):dsearchn(FF_LS,24/12.42 + 3/30),:),2,'omitnan')*FF_LS(1)))*100
sqrt(sum(mean(SPEC_MEAN(   dsearchn(FF   ,24/12.42 - 3/30):dsearchn(FF   ,24/12.42 + 3/30),:),2,'omitnan')*FF(   1)))*100



%%
%%
%% 1800 m - TIME SERIES ANALYSIS, LOMB-SCARGLE VIA MATLAB'S plomb

close all

MM_list_lats = [];
for mi=1:length(MM_list)
    MM = MM_list{mi};
    MM_list_lats(mi) = eval([MM '.P.LATITUDE;']);
end
[~,LatOrder] = sort(MM_list_lats);
LatOrder = flip(LatOrder)';
COLOR = turbo(length(MM_list));
SEGS = 3;

MM_list_char = [];
for ii = 1:length(MM_list)
    MM_list_char(ii,:) = [MM_list{ii} '  '];
end
MM_list_char = char(MM_list_char);

% All spectra together
figure('Color',[1 1 1])
figure('Color',[1 1 1])
VAR_VEC = nan(length(MM_list),1);

% For nufft method:
FF = [(1/30):(1/30):(12)]';
SPEC_SUM_VEC = nan(length(MM_list),1);
SPEC_SUM_VEC_SUPERDIURNAL_1800 = nan(length(MM_list),1);
SPEC_SUM_VEC_M2 = nan(length(MM_list),1);
SPEC_MEAN = nan(length(FF),length(MM_list));

% For plomb method:
FF_LS = FF;%[(1/90):(1/90):(12)]';
SPEC_SUM_VEC_LS = nan(length(MM_list),1);
SPEC_SUM_VEC_LS_SUPERDIURNAL_1800 = nan(length(MM_list),1);
SPEC_SUM_VEC_LS_M2 = nan(length(MM_list),1);
SPEC_LS_MEAN = nan(length(FF_LS),length(MM_list));

% For estimating uncertainty:
N_unc_dist = 10000;
Unc_dist_nufft_1800m = nan(N_unc_dist,length(MM_list));
Unc_dist_plomb_1800m = nan(N_unc_dist,length(MM_list));

for mi = (length(MM_list)-3):length(MM_list)
    MM = MM_list{mi};
    IND = [MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM & ' ...
                    '~SPIKE_IND{mi} & [' MM '.P.TIME_STERIC_HEIGHT + T0] > SF_TIME_INTERVAL(1) & [' MM '.P.TIME_STERIC_HEIGHT + T0] < SF_TIME_INTERVAL(2)'];

    [SPEC_LS,~,ERR] = ...
    nunanspectrum(eval([MM '_sh_1800(' IND ')']), ...
                  eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']),'day',...
                  'Segments',SEGS,...
                  'Window','hanning',     'Freq',FF_LS,...
                  'Plot',false,'PlotSegments',false,'Method','plomb'); hold on


    [SPEC,~,ERR] = ...
    nunanspectrum(eval([MM '_sh_1800(' IND ')']), ...
                  eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']),'day',...
                  'Segments',SEGS,...
                  'Window','hanning',     'Freq',FF,...
                  'Plot',false,'PlotSegments',false,'Method','nufft'); hold on
    VAR_VEC(mi) = var( ( eval([MM '_sh_1800(' IND ')']) ) ,'omitnan');
    figure(1)
        loglog(FF_LS,SPEC_LS,'-','Color',COLOR(dsearchn(LatOrder,mi),:),'LineWidth',1);hold on % f1
    figure(2)
        loglog(FF,SPEC,'-','Color',COLOR(dsearchn(LatOrder,mi),:),'LineWidth',1);hold on % f2

    SPEC_SUM_VEC(mi) = sum(SPEC*FF(1));
    SPEC_SUM_VEC_SUPERDIURNAL_1800(mi) = sum(SPEC(FF>=[24/24])*FF(1));
        Unc_dist_nufft_1800m(:,mi) = spec_sum_distribution(SPEC(FF>=[24/24])*FF(1), 2*SEGS, N_unc_dist);
    SPEC_SUM_VEC_M2(mi) = sum(SPEC(dsearchn(FF,[24/12.42]))*FF(1));
    SPEC_MEAN(:,mi) = SPEC;

    SPEC_SUM_VEC_LS(mi) = sum(SPEC_LS*FF_LS(1));
    SPEC_SUM_VEC_LS_SUPERDIURNAL_1800(mi) = sum(SPEC_LS(FF_LS>=[24/24])*FF_LS(1));
        Unc_dist_plomb_1800m(:,mi) = spec_sum_distribution(SPEC_LS(FF_LS>=[24/24])*FF_LS(1), 2*SEGS, N_unc_dist);
    SPEC_SUM_VEC_LS_M2(mi) = sum(SPEC_LS(dsearchn(FF_LS,[24/12.42]))*FF_LS(1));
    SPEC_LS_MEAN(:,mi) = SPEC_LS;
end
figure(1)
ylabel('Power spectral density (plomb)')
xlabel('Frequency (cpd)')
set(gca,'XScale','Log','YScale','Log','FontSize',18)
legend(MM_list{(length(MM_list)-3):length(MM_list)})

figure(2)
ylabel('Power spectral density (nunanspectrum)')
xlabel('Frequency (cpd)')
set(gca,'XScale','Log','YScale','Log','FontSize',18)
legend(MM_list{(length(MM_list)-3):length(MM_list)})

figure('Color',[1 1 1])
loglog(FF,mean(SPEC_MEAN,2,'omitnan'),'.-'); hold on%,'Color',[.2 .5 .8]
loglog(FF_LS,mean(SPEC_LS_MEAN,2,'omitnan'),'.-')%,'Color',[.8 .5 .2]
ylabel('Power spectral density (nunanspectrum)')
xlabel('Frequency (cpd)')
legend('<S(f)> = Welch with nufft','<S(f)> = Lomb-Scargle')

SPEC_MEAN_1800 = SPEC_MEAN;
SPEC_LS_MEAN_1800 = SPEC_LS_MEAN;


disp('STD in cm (NUNANSPECTRUM NUFFT):')
CM_str = repmat(' cm',length(SPEC_SUM_VEC_SUPERDIURNAL_1800),1);
[MM_list_char num2str(sqrt([SPEC_SUM_VEC_SUPERDIURNAL_1800])*100,'%.2f') CM_str]

disp('%%%%%%%%%%% SUPERDIURNAL %%%%%%%%%%%')
disp('STD in cm (PLOMB):')
CM_str = repmat(' cm',length(SPEC_SUM_VEC_LS_SUPERDIURNAL_1800),1);
[MM_list_char num2str(sqrt([SPEC_SUM_VEC_LS_SUPERDIURNAL_1800])*100,'%.2f') CM_str]
%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Variance in cm^2 (NUNANSPECTRUM NUFFT):')
CM_str = repmat(' cm^2, STD = ',length(SPEC_SUM_VEC_SUPERDIURNAL_1800),1);
[MM_list_char num2str(SPEC_SUM_VEC_SUPERDIURNAL_1800*10000,'%.2f') CM_str num2str(std(Unc_dist_nufft_1800m,1)'*10000,'%.2f')]

disp('%%%%%%%%%%% SUPERDIURNAL %%%%%%%%%%%')
disp('Variance in cm^2 (PLOMB):')
CM_str = repmat(' cm^2, STD = ',length(SPEC_SUM_VEC_LS_SUPERDIURNAL_1800),1);
[MM_list_char num2str(SPEC_SUM_VEC_LS_SUPERDIURNAL_1800*10000,'%.2f') CM_str num2str(std(Unc_dist_plomb_1800m,1)'*10000,'%.2f')]

%% Combine 600m and 1800 spectra averages onto the same plot

close all

LineWidth = 1;
YLIM = 10.^[-3 2];

EffectiveSegments = [2*3 - 1]*10;
err_high = 2*EffectiveSegments/chi2inv(.05/2,2*EffectiveSegments);
err_low = 2*EffectiveSegments/chi2inv(1-.05/2,2*EffectiveSegments);
Err = [err_low err_high]; Err = Err/Err(1);




figure('Color',[1 1 1])
loglog(FF,      10000*mean(SPEC_MEAN_600,2,'omitnan'),'.-','LineWidth',LineWidth); hold on
loglog(FF_LS,   10000*mean(SPEC_LS_MEAN_600,2,'omitnan'),'.-','LineWidth',LineWidth)

loglog(FF,      10000*mean(SPEC_MEAN_1800,2,'omitnan'),'.-','LineWidth',LineWidth)
loglog(FF_LS,   10000*mean(SPEC_LS_MEAN_1800,2,'omitnan'),'.-','LineWidth',LineWidth)

loglog([1.1 1.1]*FF(end), Err*10^-2, 'k', 'LineWidth',3)

loglog([1 1]./[12.4206012/24],YLIM,'k--')
text([24/12.4206012],3*YLIM(1),'M_2','FontSize',16)
loglog([1 1]./[6.210300601/24],YLIM,'k--')
text([24/6.210300601],3*YLIM(1),'M_4','FontSize',16)
OM = 1/(23.9344699);
loglog([1 1].*[2*OM*sind(mean(MM_list_lats))*24],YLIM,'k--')
text(24*[2*OM*sind(mean(MM_list_lats))],3*YLIM(1),'$f$','FontSize',20,'Interpreter','LaTeX')

set(gca,'FontSize',16,'XLim',[0.03 16],'YLim',YLIM)

% ylabel('\langle S(f) \rangle (cm^2)','Interpreter','TeX','FontSize',20)
ylabel('Power spectral density (cm^2 cpd^{-1})')
xlabel('Frequency (cpd)')
legend({'NUFFT, 600 m','Lomb-Scargle, 600 m', 'NUFFT, 1800 m','Lomb-Scargle, 1800 m'},...
        'Interpreter','TeX','FontSize',20,'NumColumns',2)

INPUT = input(['Do you want to save this figure? Enter any number for "yes" or\n push ' ...
    '"Enter" with a blank input for "No".\n']);

if isempty(INPUT)
    disp(['Image not saved'])
else
    exportgraphics(gcf,...
        '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/F5.pdf',...
        'BackgroundColor','none','ContentType','vector')
    disp(['Image saved'])
end

%%
%%
%% Error estimation of the variance

% Use the custom function "spec_sum_distribution" to estimate the
% distribution of the sum based on degrees of freedom of S(f).



%% Simple but inadequate approach: subtract a 24-hr running mean, and then
% % find the range of variance for 2-week intervals.
% 
% % close all;figure%$
% for mi = 1:length(MM_list)
%     MM = MM_list{mi};
%     FORNIGHTLY_SEGMENTS = datenum('2023-04-01'):14:datenum('2023-07-22');
% 
%     LOW_PASS = nan(size(eval([MM '_sh_600'])));
%     for ii=1:length(LOW_PASS)
%         PlusMinus24hr = dsearchn(eval([MM '.P.TIME_STERIC_HEIGHT']), eval([MM '.P.TIME_STERIC_HEIGHT(ii)']) - 12/24) : ...
%                         dsearchn(eval([MM '.P.TIME_STERIC_HEIGHT']), eval([MM '.P.TIME_STERIC_HEIGHT(ii)']) + 12/24);
%         LOW_PASS(ii) = mean(eval([MM '_sh_600([' num2str(PlusMinus24hr) '])']) , 'omitnan');
%     end
% 
%     VAR = [];
%     for si = 1:[length(FORNIGHTLY_SEGMENTS) - 1]
%         IND = [MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM & ' ...
%             '~SPIKE_IND{mi} & [' MM '.P.TIME_STERIC_HEIGHT + T0] > FORNIGHTLY_SEGMENTS(si) & [' MM '.P.TIME_STERIC_HEIGHT + T0] < FORNIGHTLY_SEGMENTS(si+1)'];
% 
%         % Linear fit approach % % % % % % % % % % % %
%         % HH = [ones(size(eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']))) , ...
%         %       eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']) - mean(eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']),'omitnan')];
%         % LINEAR_FIT = HH*[HH\eval([MM '_sh_600(' IND ')'])];
%         % 
%         % plot( eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']) , eval([MM '_sh_600(' IND ')']) - LINEAR_FIT, '.-'); hold on%$
%         % % % % % % % % % % % % % % % % % % % % % % %
% 
% 
%         % % Boxcar filter high-pass approach  % % % % %
%         % plot( eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']) , ...
%         %       eval([MM '_sh_600(' IND ')']) - eval(['LOW_PASS(' IND ')']), '.-'); hold on%$
%         % % % % % % % % % % % % % % % % % % % % % % % %
% 
%         VAR(si) = var(eval([MM '_sh_600(' IND ')']) - eval(['LOW_PASS(' IND ')']) , 'omitnan');
%     end
%     eval([MM '_sh_600_var = VAR'';']);
% end
% 
% for mi = 7:length(MM_list)
%     MM = MM_list{mi};
%     FORNIGHTLY_SEGMENTS = datenum('2023-04-01'):14:datenum('2023-07-22');
% 
%     LOW_PASS = nan(size(eval([MM '_sh_1800'])));
%     for ii=1:length(LOW_PASS)
%         PlusMinus24hr = dsearchn(eval([MM '.P.TIME_STERIC_HEIGHT']), eval([MM '.P.TIME_STERIC_HEIGHT(ii)']) - 12/24) : ...
%                         dsearchn(eval([MM '.P.TIME_STERIC_HEIGHT']), eval([MM '.P.TIME_STERIC_HEIGHT(ii)']) + 12/24);
%         LOW_PASS(ii) = mean(eval([MM '_sh_1800([' num2str(PlusMinus24hr) '])']) , 'omitnan');
%     end
% 
%     VAR = [];
%     for si = 1:[length(FORNIGHTLY_SEGMENTS) - 1]
%         IND = [MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM & ' ...
%             '~SPIKE_IND{mi} & [' MM '.P.TIME_STERIC_HEIGHT + T0] > FORNIGHTLY_SEGMENTS(si) & [' MM '.P.TIME_STERIC_HEIGHT + T0] < FORNIGHTLY_SEGMENTS(si+1)'];
% 
%         VAR(si) = var(eval([MM '_sh_1800(' IND ')']) - eval(['LOW_PASS(' IND ')']) , 'omitnan');
%     end
%     eval([MM '_sh_1800_var = VAR'';']);
% end
% %%
% close all
% 
% figure
% for mi = 1:length(MM_list)
%     MM = MM_list{mi};
%     subplot(10,1,mi)
%     plot(eval([MM '_sh_600_var'])*100*100,'.-'); ylim([0 2]); hold on
%     title([MM ': \mu = ' num2str(mean(eval([MM '_sh_600_var'])*100*100,'omitnan')) ...
%            ', \sigma = ' num2str(std(eval([MM '_sh_600_var'])*100*100,'omitnan'))])
% end
% % 600 m variance over 14 day periods
% 
% figure
% for mi = 7:length(MM_list)
%     MM = MM_list{mi};
%     subplot(4,1,mi-7+1)
%     plot(eval([MM '_sh_1800_var'])*100*100,'.-'); ylim([0 2]); hold on
%     title([MM ': \mu = ' num2str(mean(eval([MM '_sh_1800_var'])*100*100,'omitnan')) ...
%            ', \sigma = ' num2str(std(eval([MM '_sh_1800_var'])*100*100,'omitnan'))])
% end
% % 1800 m variance over 14 day periods

%%
%%
%%
%% 600 m but no segmenting or windowing, in order to assess those steps' effects

close all

warning('The results below are for unsegmented and unwindowed time series and are meant for verification of methodology.')

MM_list_char = [];
for ii = 1:length(MM_list)
    MM_list_char(ii,:) = [MM_list{ii} '  '];
end
MM_list_char = char(MM_list_char);

MM_list_lats = [];
for mi=1:length(MM_list)
    MM = MM_list{mi};
    MM_list_lats(mi) = eval([MM '.P.LATITUDE;']);
end
[~,LatOrder] = sort(MM_list_lats);
LatOrder = flip(LatOrder)';
COLOR = turbo(length(MM_list));

% All spectra together
figure('Color',[1 1 1])
figure('Color',[1 1 1])
VAR_VEC = nan(length(MM_list),1);

% For nufft method:
FF = [(1/90):(1/90):(12)]';
SPEC_SUM_VEC = nan(length(MM_list),1);
SPEC_SUM_VEC_SUPERDIURNAL = nan(length(MM_list),1);
SPEC_SUM_VEC_M2 = nan(length(MM_list),1);
SPEC_MEAN = nan(length(FF),length(MM_list));

% For plomb method:
FF_LS = FF;%[(1/90):(1/90):(12)]';
SPEC_SUM_VEC_LS = nan(length(MM_list),1);
SPEC_SUM_VEC_LS_SUPERDIURNAL = nan(length(MM_list),1);
SPEC_SUM_VEC_LS_M2 = nan(length(MM_list),1);
SPEC_LS_MEAN = nan(length(FF_LS),length(MM_list));

for mi = 1:length(MM_list)
    MM = MM_list{mi};
    IND = [MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH > ' MM(1) 'P_RANGE_LIM & ' ...
                    '~SPIKE_IND{mi} & [' MM '.P.TIME_STERIC_HEIGHT + T0] > SF_TIME_INTERVAL(1) & [' MM '.P.TIME_STERIC_HEIGHT + T0] < SF_TIME_INTERVAL(2)'];

    % [SPEC_LS,~] = plomb(eval([MM '.P.STERIC_HEIGHT_ANOMALY(' IND ')']),...
    %                     eval([MM '.P.TIME_STERIC_HEIGHT('    IND ')']), FF_LS);
    % [SPEC_LS,~] = plomb(detrend(eval([MM '.P.STERIC_HEIGHT_ANOMALY(' IND ')'])),...
    %                             eval([MM '.P.TIME_STERIC_HEIGHT('    IND ')']), FF_LS);
    [SPEC_LS,~,ERR] = ...
    nunanspectrum(eval([MM '_sh_600(' IND ')']), ...
                  eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']),'day',...
                  'Segments',1,...
                  'Window','rectwin',     'Freq',FF_LS,...
                  'Plot',false,'PlotSegments',false,'Method','plomb'); hold on


    [SPEC,~,ERR] = ...
    nunanspectrum(eval([MM '_sh_600(' IND ')']), ...
                  eval([MM '.P.TIME_STERIC_HEIGHT(' IND ')']),'day',...
                  'Segments',1,...
                  'Window','rectwin',     'Freq',FF,...
                  'Plot',false,'PlotSegments',false,'Method','nufft'); hold on
    % VAR_VEC(mi) = var(detrend( eval([MM '.P.STERIC_HEIGHT_ANOMALY(' IND ')']) ) ,'omitnan');
    VAR_VEC(mi) = var( ( eval([MM '_sh_600(' IND ')']) ) ,'omitnan');
    figure(1)
        loglog(FF_LS,SPEC_LS,'-','Color',COLOR(dsearchn(LatOrder,mi),:),'LineWidth',1);hold on
    figure(2)
        loglog(FF,SPEC,'-','Color',COLOR(dsearchn(LatOrder,mi),:),'LineWidth',1);hold on

    SPEC_SUM_VEC(mi) = sum(SPEC*FF(1));
    SPEC_SUM_VEC_SUPERDIURNAL(mi) = sum(SPEC(FF>=[24/24])*FF(1));
    SPEC_SUM_VEC_M2(mi) = sum(SPEC(dsearchn(FF,[24/12.42]))*FF(1));
    SPEC_MEAN(:,mi) = SPEC;

    SPEC_SUM_VEC_LS(mi) = sum(SPEC_LS*FF_LS(1));
    SPEC_SUM_VEC_LS_SUPERDIURNAL(mi) = sum(SPEC_LS(FF_LS>=[24/24])*FF_LS(1));
    SPEC_SUM_VEC_LS_M2(mi) = sum(SPEC_LS(dsearchn(FF_LS,[24/12.42]))*FF_LS(1));
    SPEC_LS_MEAN(:,mi) = SPEC_LS;
end
figure(1)
ylabel('Power spectral density (plomb)')
xlabel('Frequency (cpd)')
set(gca,'XScale','Log','YScale','Log','FontSize',18)
legend(MM_list)

figure(2)
ylabel('Power spectral density (nunanspectrum)')
xlabel('Frequency (cpd)')
set(gca,'XScale','Log','YScale','Log','FontSize',18)
legend(MM_list)

figure('Color',[1 1 1])
loglog(FF,mean(SPEC_MEAN,2,'omitnan'),'.-'); hold on%,'Color',[.2 .5 .8]
loglog(FF_LS,mean(SPEC_LS_MEAN,2,'omitnan'),'.-')%,'Color',[.8 .5 .2]
ylabel('Power spectral density (nunanspectrum)')
xlabel('Frequency (cpd)')
legend('<S(f)> = Welch with nufft','<S(f)> = Lomb-Scargle')

SPEC_MEAN_600 = SPEC_MEAN;
SPEC_LS_MEAN_600 = SPEC_LS_MEAN;

disp('STD in cm (NUNANSPECTRUM NUFFT):')
CM_str = repmat(' cm',length(SPEC_SUM_VEC_SUPERDIURNAL),1);
[MM_list_char num2str(sqrt([SPEC_SUM_VEC_SUPERDIURNAL])*100,'%.2f') CM_str]

disp('%%%%%%%%%%% SUPERDIURNAL %%%%%%%%%%%')
disp('STD in cm (PLOMB):')
CM_str = repmat(' cm',length(SPEC_SUM_VEC_LS_SUPERDIURNAL),1);
[MM_list_char num2str(sqrt([SPEC_SUM_VEC_LS_SUPERDIURNAL])*100,'%.2f') CM_str]

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Variance in cm^2 (NUNANSPECTRUM NUFFT):')
CM_str = repmat(' cm^2',length(SPEC_SUM_VEC_SUPERDIURNAL),1);
[MM_list_char num2str(SPEC_SUM_VEC_SUPERDIURNAL*10000,'%.2f') CM_str]

disp('%%%%%%%%%%% SUPERDIURNAL %%%%%%%%%%%')
disp('Variance in cm^2 (PLOMB):')
CM_str = repmat(' cm^2',length(SPEC_SUM_VEC_LS_SUPERDIURNAL),1);
[MM_list_char num2str(SPEC_SUM_VEC_LS_SUPERDIURNAL*10000,'%.2f') CM_str]

warning('The results above are for unsegmented and unwindowed time series.')

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

function errorbar_stdofmean(bin_center, sep_vec, sf_vec, bin_edges, varargin)
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

% Given a Spectrum, get the distribution of its sum based on the assumption
% that each of its elements has a chi-squared distribution but of different
% magnitudes.
% 
% IN:
% SS = Power spectrum, vector. Multiply by its fundamental frequency in
%      order to turn its sum into variance
% kk = number of degrees of freedom, which should be the same for each and
%      every element of SS. With overlapping segments multiplied by the
%      appropriate window, kk = 2*num_of_segments (2 for the 
% NN = Desired length of the output distribution.
% 
% OUT:
% DD = NN-element-long vector that is approximately distributed how we
%      would expect sum(SS) to be distributed
function DD = spec_sum_distribution(SS,kk,NN)
% We might find ourselves in the situation where length(SS)*NN is too
% large, so loop through to enable large NN:
DD = zeros(NN,1);
DIST_OUT = [0:[1/NN]:[1 - 1/NN]]';
for si = 1:length(SS)
    CHI2INV_si = chi2inv(DIST_OUT,kk)*SS(si)/kk;
    DD = DD + CHI2INV_si(randperm(length(CHI2INV_si)));
end
end
