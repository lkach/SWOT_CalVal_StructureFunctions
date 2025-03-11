% Grid steric height contribution from both deep and profiler observations


%% Load P-moorings (which only have profilers)

% Add necessary file paths before the rest of the script
addpath ~; lkaddpath

% % % LOAD PROFILERS
% FOLDER_P = '/Users/kachelein/Documents/JPL/work/07.2023_QC_moorings/PROFILERS/';
% FOLDER_P = '/mnt/flow/swot/calval/Data/insitu_mooring/JPL_QC/PROFILERS/';
FOLDER_P = '/mnt/flow/swot/Analysis_Luke/Moorings2/DATA/PROF/';
FOLDER_P = '/Users/kachelein/Documents/JPL/work/02.2024/OSM/DATA/PROF/'; % local
DIR_P = dir([FOLDER_P 'SWOTPOSTLAUNCH_L2_JPLQC_MOORING-P*5s4s_alltimes.nc']);
cd(FOLDER_P)

TIME_0 = datenum('1950-01-01 00:00:00');
MOORING_NAMES = cell(length(DIR_P),1);

for ii = 1:length(DIR_P)
    MOORING_NAME = DIR_P(ii).name(33:34);
    MOORING_NAMES{ii} = MOORING_NAME;
    eval([MOORING_NAME,'p = ncreadall(DIR_P(ii).name);']);%%%
    eval([MOORING_NAME,'p.TIME_datenum = ',MOORING_NAME,'p.TIME + TIME_0;']);
    disp([MOORING_NAME,'p (profiler) loaded']);
end
wsm %%% custom fucntion "WorkSpaceMemory" to show how much memory is being used

MOORING_NAMES % unsuppressed for validation purposes

%% Variables:

    %                          TIME: [347835×1 double]
    %                      LATITUDE: 36.0930
    %                     LONGITUDE: -125.0900
    %                          TEMP: [347835×1 double]
    %                          PSAL: [347835×1 double]
    %                          CNDC: [347835×1 double]
    %                          PRES: [347835×1 double]
    %                       TEMP_QC: [347835×1 int8]
    %                       PSAL_QC: [347835×1 int8]
    %                       CNDC_QC: [347835×1 int8]
    %                       PRES_QC: [347835×1 int8]
    %                       QC_FLAG: [347835×1 double]
    %                         DEPTH: [347835×1 double]
    %                      PROF_NUM: [347835×1 double]
    %                           RHO: [347835×1 double]
    %         STERIC_HEIGHT_ANOMALY: [1951×1 double]
    %                  RHO_TIMEMEAN: [501×1 double]
    %            STERIC_HEIGHT_MEAN: 0.8550
    %            RHO_TIMEMEAN_DEPTH: [501×1 double]
    %  LATITUDE_GPSSB_STERIC_HEIGHT: [1951×1 double]
    % LONGITUDE_GPSSB_STERIC_HEIGHT: [1951×1 double]
    %            TIME_STERIC_HEIGHT: [1951×1 double]
    %        PROF_GOOD_REALIZATIONS: [1951×1 double]
    %    STERIC_HEIGHT_CUTOFF_DEPTH: 500
    %  STERIC_HEIGHT_MIN_PROF_DEPTH: [1951×1 double]
    %  STERIC_HEIGHT_MAX_PROF_DEPTH: [1951×1 double]
    %                         RHO_0: 1029
    %                  TIME_datenum: [347835×1 double]

%% Times covered:

% datestr(TIME_0 + P1p.TIME(1))
% datestr(TIME_0 + P2p.TIME(1))
% datestr(TIME_0 + P3p.TIME(1))
% datestr(TIME_0 + P4p.TIME(1))
% datestr(TIME_0 + P5p.TIME(1))
% datestr(TIME_0 + P6p.TIME(1))
% datestr(TIME_0 + P7p.TIME(1))
%     '20-Feb-2023 00:00:08'
%     '20-Feb-2023 00:46:56'
%     '20-Feb-2023 16:57:39'
%     '20-Feb-2023 00:00:06'
%     '21-Feb-2023 00:11:01'
%     '20-Feb-2023 02:39:22'
%     '20-Feb-2023 22:56:05'
% 
% datestr(TIME_0 + P1p.TIME(end))
% datestr(TIME_0 + P2p.TIME(end))
% datestr(TIME_0 + P3p.TIME(end))
% datestr(TIME_0 + P4p.TIME(end))
% datestr(TIME_0 + P5p.TIME(end))
% datestr(TIME_0 + P6p.TIME(end))
% datestr(TIME_0 + P7p.TIME(end))
%     '02-Aug-2023 18:27:47'
%     '02-Aug-2023 21:08:06'
%     '02-Aug-2023 18:53:27'
%     '02-Aug-2023 16:21:42'
%     '28-Jun-2023 07:53:46'
%     '02-Aug-2023 19:08:56'
%     '09-Jul-2023 18:08:39'

%% Set gridding parameters

% This is strictly a time-domain gridding, but the code I have to implement
% OI is necessarily 2D, so the gridding coordinates will be linear in time
% and constant in "depth" (or more accurately some fake variable).

t_decor = 0.5/24; % units of days
dt_grid =   0.5/24;

% P1p.TIME_GRID = [datenum('20-Feb-2023 00:00:00') - TIME_0]:dt_grid:...
%                 [datenum('02-Aug-2023 18:30:00') - TIME_0];
% P1p.TIME_GRID = P1p.TIME_GRID';
% 
% P2p.TIME_GRID = [datenum('20-Feb-2023 00:00:00') - TIME_0]:dt_grid:...
%                 [datenum('02-Aug-2023 21:00:00') - TIME_0];
% P2p.TIME_GRID = P2p.TIME_GRID';
% 
% P3p.TIME_GRID = [datenum('20-Feb-2023 17:00:00') - TIME_0]:dt_grid:...
%                 [datenum('02-Aug-2023 18:30:00') - TIME_0];
% P3p.TIME_GRID = P3p.TIME_GRID';
% 
% P4p.TIME_GRID = [datenum('20-Feb-2023 00:00:00') - TIME_0]:dt_grid:...
%                 [datenum('02-Aug-2023 16:30:00') - TIME_0];
% P4p.TIME_GRID = P4p.TIME_GRID';
% 
% P5p.TIME_GRID = [datenum('21-Feb-2023 00:00:00') - TIME_0]:dt_grid:...
%                 [datenum('28-Jun-2023 07:30:00') - TIME_0];
% P5p.TIME_GRID = P5p.TIME_GRID';
% 
% P6p.TIME_GRID = [datenum('20-Feb-2023 03:00:00') - TIME_0]:dt_grid:...
%                 [datenum('02-Aug-2023 19:00:00') - TIME_0];
% P6p.TIME_GRID = P6p.TIME_GRID';
% 
% P7p.TIME_GRID = [datenum('20-Feb-2023 23:00:00') - TIME_0]:dt_grid:...
%                 [datenum('09-Jul-2023 18:00:00') - TIME_0];
% P7p.TIME_GRID = P7p.TIME_GRID';

P1p.TIME_GRID = [datenum('19-Feb-2023 00:00:00') - TIME_0]:dt_grid:...
                [datenum('07-Oct-2023 00:00:00') - TIME_0];
P1p.TIME_GRID = P1p.TIME_GRID';
P2p.TIME_GRID = P1p.TIME_GRID;
P3p.TIME_GRID = P1p.TIME_GRID;
%P4p.TIME_GRID = P1p.TIME_GRID;
P5p.TIME_GRID = P1p.TIME_GRID;
P6p.TIME_GRID = P1p.TIME_GRID;
P7p.TIME_GRID = P1p.TIME_GRID;

%% Perform the optimal interpolation

% S#.STERIC_HEIGHT_ANOMALY_GRID_PROF
% S#.STERIC_HEIGHT_ANOMALY_GRID_FIXED
% S#.STERIC_HEIGHT_ANOMALY_GRID_TOTAL

GOOD_PROF_CUTOFF = 150;
PROF_LENGTH_CUTOFF = 415;

for ii = 1:length(DIR_P)
    PP = MOORING_NAMES{ii};

    eval([PP 'p.STERIC_HEIGHT_ANOMALY_GRID_PROF = nan(size(' PP 'p.TIME_GRID));']);

    for ti = 1:length(eval([PP 'p.TIME_GRID']))
        % % % Grid fixed steric height anomaly (in practice, all matrices in the OI inversion are singular, so instead use bin averaging):
        % 
        % % Grid profiler steric height anomaly:
        % INDp = eval(['[' SS 'p.TIME_MIDPROFILE > ' SS '.TIME_GRID(ti)-dt_grid  &  ' SS 'p.TIME_MIDPROFILE < ' SS '.TIME_GRID(ti)+dt_grid]']); %  & ~' SS 'p.QC_FLAG
        % if sum(INDp)>0 && sum(isfinite(eval([SS 'p.STERIC_HEIGHT_ANOMALY(INDp)'])))>0
        %     STERIC_HEIGHT_ANOMALY_GRID_p = ObjMap(...
        %         eval([SS 'p.STERIC_HEIGHT_ANOMALY(INDp)']), ... % data
        %         eval([SS 'p.TIME_MIDPROFILE(INDp)']), ... % x
        %         eval([SS 'p.TIME_MIDPROFILE(INDp)*0']), ... % y, fake variable
        %         t_decor, ... % xdecor
        %         1, ... % ydecor, fake variable
        %         eval([SS '.TIME_GRID(ti)']), ... % xgrid
        %         0, ... % ygrid, fake variable
        %         0, ... % theta
        %         0.1, ... % noise2sig
        %         0); % plotboolean
        % else
        %     disp([SS ': no valid data at ti = ' num2str(ti)])
        %     STERIC_HEIGHT_ANOMALY_GRID_p = nan;
        % end

        % Grid profiler steric height anomaly:
        INDp = eval(['[' PP 'p.TIME_STERIC_HEIGHT > ' PP 'p.TIME_GRID(ti)-dt_grid  &  ' PP 'p.TIME_STERIC_HEIGHT < ' PP 'p.TIME_GRID(ti)+dt_grid  &  ' ...
                         PP 'p.PROF_GOOD_REALIZATIONS > ' num2str(GOOD_PROF_CUTOFF) ...
                         ' & [' PP 'p.STERIC_HEIGHT_MAX_PROF_DEPTH - ' PP 'p.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' num2str(PROF_LENGTH_CUTOFF) ']']); %  & ~' SS 'p.QC_FLAG
        % sum(INDp)
        if sum(INDp)>0 && sum(isfinite(eval([PP 'p.STERIC_HEIGHT_ANOMALY(INDp)'])))>0
            STERIC_HEIGHT_ANOMALY_GRID_p = mean(eval([PP 'p.STERIC_HEIGHT_ANOMALY(INDp)']),'omitnan');
        else
            disp([PP 'p: no valid data at ti = ' num2str(ti)])
            STERIC_HEIGHT_ANOMALY_GRID_p = nan;
        end

        % Save variables before finally saving them into the netCDT parent files
        eval([PP 'p.STERIC_HEIGHT_ANOMALY_GRID_PROF(ti) = STERIC_HEIGHT_ANOMALY_GRID_p;']);
    end
    eval([PP 'p.STERIC_HEIGHT_MEAN_PROF = ' PP 'p.STERIC_HEIGHT_MEAN;']);
end

%% Save into netCDF files

% error('Make sure that this works by testing manually.')

for ii = 1:length(DIR_P)
    PP = MOORING_NAMES{ii};

    % Save time grid to profiler files
    save_to_nc([DIR_P(ii).folder '/' DIR_P(ii).name],... file name
        'TIME_GRID',... name of new variable
        eval([PP 'p.TIME_GRID']),... data
        {'TIME_GRID',length(eval([PP 'p.TIME_GRID']))},... DIM_METADATA
        {{'long_name',['Time at which steric height is gridded via interval averaging (+-' num2str(dt_grid) ' days).']},...
        {'standard_name','time_grid'},...
        {'units','days since 1950-01-01T00:00:00Z'},...
        {'coordinates','TIME LONGITUDE LATITUDE'}},... METADATA
        true,true);% OVERWRITE, SKIP, MULTIPLE_DIMS

    % save profiler SH to profiler files
    save_to_nc([DIR_P(ii).folder '/' DIR_P(ii).name],... file name
        'STERIC_HEIGHT_ANOMALY_GRID_PROF',... name of new variable
        eval([PP 'p.STERIC_HEIGHT_ANOMALY_GRID_PROF']),... data
        {'TIME_GRID'},... DIM_METADATA
        {{'long_name',['Steric height from profilers (above 500 m) gridded via interval averaging (+-' num2str(dt_grid) ' days).']},...
        {'standard_name','steric_height_grid_prof'},...
        {'units','meters'},...
        {'coordinates','TIME LONGITUDE LATITUDE'}},... METADATA
        true,true);% OVERWRITE, SKIP, MULTIPLE_DIMS

end

% Done seprately in case there is an issue:
for ii = 1:length(DIR_P)
    PP = MOORING_NAMES{ii};
    % save GOOD_PROF_CUTOFF to profiler file. This quantity indicates the
    % minimum number of good data per cast that was required for each
    % steric height estimate used for gridding.
    save_to_nc([DIR_P(ii).folder '/' DIR_P(ii).name],... file name
        'GOOD_PROF_CUTOFF',... name of new variable
        GOOD_PROF_CUTOFF,... data
        {'SCALAR'},... DIM_METADATA
        {{'long_name','Minimum number of good data per profiler cast that was required for each steric height estimate used for gridding.'},...
        {'standard_name','good_prof_cutoff'},...
        {'units','none'}},... METADATA
        true,true);% OVERWRITE, SKIP, MULTIPLE_DIMS

    save_to_nc([DIR_P(ii).folder '/' DIR_P(ii).name],... file name
        'PROF_LENGTH_CUTOFF',... name of new variable
        PROF_LENGTH_CUTOFF,... data
        {'SCALAR'},... DIM_METADATA
        {{'long_name','Minimum length (max-min) of profiler cast that was required for each steric height estimate used for gridding.'},...
        {'standard_name','good_prof_cutoff'},...
        {'units','meters'}},... METADATA
        true,true);% OVERWRITE, SKIP, MULTIPLE_DIMS
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Done with (interval average) gridding steric height and saving to files. ' datestr(now)])


%%


%% (Old format) Variables:


% % % PROFILERS
    %                          TIME: [339608x1 double]
    %                      LATITUDE: 36.0950
    %                     LONGITUDE: -125.0900
    %                          CNDC: [339608x1 double]
    %                          TEMP: [339608x1 double]
    %                          PRES: [339608x1 double]
    %                          PSAL: [339608x1 double]
    %                         DEPTH: [339608x1 double]
    %                      PROF_NUM: [339608x1 double]
    %        STERIC_HEIGHT_ABOVE500: [1964x1 double]
    %           PRESSURE_MIDPROFILE: 247.5000
    %                PRESSURE_RANGE: [2x1 double]
    %               TIME_MIDPROFILE: [1964x1 double]
    %                       QC_FLAG: [339608x1 uint8]
    %                    INSTR_INFO: ' '
    %                    INSTR_MAKE: [20x1 char]
    %                   INSTR_MODEL: [20x1 char]
    %                      INSTR_SN: [20x1 char]
    %         TIME_GPS_SURFACE_BUOY: [162x1 double]
    %    LONGITUDE_GPS_SURFACE_BUOY: [162x1 double]
    %     LATITUDE_GPS_SURFACE_BUOY: [162x1 double]
    %                           RHO: [339608x1 double]
    %         STERIC_HEIGHT_ANOMALY: [1964x1 double]
    %            STERIC_HEIGHT_MEAN: -0.8377
    %                  RHO_TIMEMEAN: [501x1 double]
    %            RHO_TIMEMEAN_DEPTH: [501x1 double]
    %  LATITUDE_GPSSB_STERIC_HEIGHT: [1964x1 double]
    % LONGITUDE_GPSSB_STERIC_HEIGHT: [1964x1 double]
    %            TIME_STERIC_HEIGHT: [1964x1 double]
    %        PROF_GOOD_REALIZATIONS: [1964x1 double]
    %    STERIC_HEIGHT_CUTOFF_DEPTH: 445
    %  STERIC_HEIGHT_MIN_PROF_DEPTH: [1964x1 double]
    %  STERIC_HEIGHT_MAX_PROF_DEPTH: [1964x1 double]
    %                         RHO_0: 1029
    %                  TIME_datenum: [339608x1 double]