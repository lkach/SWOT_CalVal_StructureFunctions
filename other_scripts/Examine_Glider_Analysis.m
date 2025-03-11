%% Supplemental Material
%% Visualize glider data
G_DIR = '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/GLIDERS';
% -rw-r--r--  1 kachelein  staff   316M Jun  7 16:47 ru32_L2_Processed_20230412T10_to_20230713T15.nc
% -rw-r--r--  1 kachelein  staff   193K Jun  7 16:47 ru32_L2_Processed_600m_20230412T10_to_20230713T15.nc
% -rw-r--r--  1 kachelein  staff   193K Jun  7 16:47 ru32_L2_Processed_1000m_20230412T10_to_20230713T15.nc
% -rw-r--r--  1 kachelein  staff   254M Jun  7 16:47 ru38_L2_Processed_20230430T06_to_20230711T23.nc
% -rw-r--r--  1 kachelein  staff   159K Jun  7 16:47 ru38_L2_Processed_600m_20230430T06_to_20230711T23.nc
% -rw-r--r--  1 kachelein  staff   159K Jun  7 16:47 ru38_L2_Processed_1000m_20230430T06_to_20230711T23.nc

ru32_600 = ncreadall('ru32_L2_Processed_600m_20230412T10_to_20230713T15.nc');
% ru32_1000 = ncreadall('ru32_L2_Processed_1000m_20230412T10_to_20230713T15.nc');

ru38_600 = ncreadall('ru38_L2_Processed_600m_20230430T06_to_20230711T23.nc');
% ru38_1000 = ncreadall('ru38_L2_Processed_1000m_20230430T06_to_20230711T23.nc');


%% Plot SH

close all

T0 = datenum('1950-01-01');

figure
plot(T0+ru32_600.TIME_STERIC_HEIGHT, ru32_600.STERIC_HEIGHT_ANOMALY, '.-'); hold on
plot(T0+ru32_1000.TIME_STERIC_HEIGHT, ru32_1000.STERIC_HEIGHT_ANOMALY, '.-')
datetick

figure
plot(T0+ru38_600.TIME_STERIC_HEIGHT, ru38_600.STERIC_HEIGHT_ANOMALY, '.-'); hold on
plot(T0+ru38_1000.TIME_STERIC_HEIGHT, ru38_1000.STERIC_HEIGHT_ANOMALY, '.-')
datetick

%% Plot SH

close all

T0 = datenum('1950-01-01');

figure
plot(T0+ru32_600.TIME_STERIC_HEIGHT, ru32_600.STERIC_HEIGHT_ANOMALY, '.-'); hold on
plot(T0+ru32_1000.TIME_STERIC_HEIGHT, ru32_1000.STERIC_HEIGHT_ANOMALY, '.-')
datetick

figure
plot(T0+ru38_600.TIME_STERIC_HEIGHT, ru38_600.STERIC_HEIGHT_ANOMALY, '.-'); hold on
plot(T0+ru38_1000.TIME_STERIC_HEIGHT, ru38_1000.STERIC_HEIGHT_ANOMALY, '.-')
datetick



%% Plot profile depths

close all

figure
subplot(211)
plot(T0+ru32_600.TIME_STERIC_HEIGHT, ru32_600.STERIC_HEIGHT_MAX_PROF_DEPTH, '.-'); hold on
plot(T0+ru32_1000.TIME_STERIC_HEIGHT, ru32_1000.STERIC_HEIGHT_MAX_PROF_DEPTH, '.-')
title('RU32: MAX (600 & 1000)')
subplot(212)
plot(T0+ru32_600.TIME_STERIC_HEIGHT, ru32_600.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-'); hold on
plot(T0+ru32_1000.TIME_STERIC_HEIGHT, ru32_1000.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')

figure
subplot(211)
plot(T0+ru38_600.TIME_STERIC_HEIGHT, ru38_600.STERIC_HEIGHT_MAX_PROF_DEPTH, '.-'); hold on
plot(T0+ru38_1000.TIME_STERIC_HEIGHT, ru38_1000.STERIC_HEIGHT_MAX_PROF_DEPTH, '.-')
title('RU32: MAX (600 & 1000)')
subplot(212)
plot(T0+ru38_600.TIME_STERIC_HEIGHT, ru38_600.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-'); hold on
plot(T0+ru38_1000.TIME_STERIC_HEIGHT, ru38_1000.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')

%% Plot full profile span

close all

figure
plot(T0+ru32_600.TIME_STERIC_HEIGHT, ru32_600.STERIC_HEIGHT_MAX_PROF_DEPTH - ru32_600.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-'); hold on
plot(T0+ru32_1000.TIME_STERIC_HEIGHT, ru32_1000.STERIC_HEIGHT_MAX_PROF_DEPTH - ru32_1000.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
title('RU32: MAX-MIN (600 & 1000)')

figure
plot(T0+ru38_600.TIME_STERIC_HEIGHT, ru38_600.STERIC_HEIGHT_MAX_PROF_DEPTH - ru38_600.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-'); hold on
plot(T0+ru38_1000.TIME_STERIC_HEIGHT, ru38_1000.STERIC_HEIGHT_MAX_PROF_DEPTH - ru38_1000.STERIC_HEIGHT_MIN_PROF_DEPTH, '.-')
title('RU38: MAX-MIN (600 & 1000)')


%%
%%
%% Additionally, verify basic information about the gliders:

% clear memory if necessary

G_DIR = '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/GLIDERS';
ru32 = ncreadall('ru32_L2_Processed_20230412T10_to_20230713T15.nc');
ru38 = ncreadall('ru38_L2_Processed_20230430T06_to_20230711T23.nc');
wsm % ~600 MB

%%
    %                         TIME: [3760581×1 double]
    %                    LONGITUDE: [3760581×1 double]
    %                     LATITUDE: [3760581×1 double]
    %                         PRES: [3760581×1 double]
    %                        DEPTH: [3760581×1 double]
    %                         TEMP: [3760581×1 double]
    %                         CNDC: [3760581×1 double]
    %                         PSAL: [3760581×1 double]
    %                          RHO: [3760581×1 double]
    %                        RHO_0: 1029
    %                      QC_FLAG: [3760581×1 double]
    %                  PROFILE_NUM: [3760581×1 double]
    %              TIME_MIDPROFILE: [1693×1 double]
    %             DURATION_PROFILE: [1693×1 double]
    %            LONGITUDE_PROFILE: [1693×1 double]
    %             LATITUDE_PROFILE: [1693×1 double]
    %          PRESSURE_MIDPROFILE: [1693×1 double]
    %                  UPDOWN_CAST: [1693×1 double]
    %           TIME_STERIC_HEIGHT: [1693×1 double]
    %        STERIC_HEIGHT_ANOMALY: [1693×1 double]
    % STERIC_HEIGHT_MAX_PROF_DEPTH: [1693×1 double]
    % STERIC_HEIGHT_MIN_PROF_DEPTH: [1693×1 double]
    %       PROF_GOOD_REALIZATIONS: [1693×1 double]
    %   STERIC_HEIGHT_CUTOFF_DEPTH: 500
% TIME:
% units     = 'days since 1950-01-01 00:00:00.0'
%%

close all
T0 = datenum('1950-01-01 00:00:00');

figure
plot(ru32.TIME + T0, ru32.DEPTH .*[~ru32.QC_FLAG]./[~ru32.QC_FLAG],'.-'); hold on
plot(ru38.TIME + T0, ru38.DEPTH .*[~ru38.QC_FLAG]./[~ru38.QC_FLAG],'.-')
datetick

disp('RU32 date range:')
disp(datestr(ru32.TIME([1 end]) + T0))
disp('RU38 date range:')
disp(datestr(ru38.TIME([1 end]) + T0))
disp(' ')
disp('RU32 flags (first and last):')
ru32.QC_FLAG([1 end])
disp('RU38 flags (first and last):')
ru38.QC_FLAG([1 end])
disp(' ')
disp('RU32 max. num. profiles:')
max(ru32.PROFILE_NUM)
disp('RU38 max. num. profiles:')
max(ru38.PROFILE_NUM)
