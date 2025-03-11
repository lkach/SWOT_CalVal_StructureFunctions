% Grid steric height contribution from both deep and profiler observations

%% Load

% Add necessary file paths before the rest of the script
addpath ~; lkaddpath

% Replace with your paths where the .nc files including density are stored
% FOLDER = '/Users/kachelein/Documents/JPL/work/07.2023_QC_moorings/FIXED_CTD/'; % local
% FOLDER = '/mnt/flow/swot/calval/Data/insitu_mooring/JPL_QC/FIXED_CTD/'; % server
FOLDER = '/mnt/flow/swot/Analysis_Luke/Moorings2/DATA/FIXED/';
FOLDER = '/Users/kachelein/Documents/JPL/work/02.2024/OSM/DATA/FIXED/'; % local
DIR = dir([FOLDER 'SWOTPOSTLAUNCH_L2_JPLQC_MOORING-S*.nc']);
cd(FOLDER)

TIME_0 = datenum('1950-01-01 00:00:00');
MOORING_NAMES = cell(length(DIR),1);

for ii = 1:length(DIR)
    MOORING_NAME = DIR(ii).name(33:34);
    MOORING_NAMES{ii} = MOORING_NAME;
    eval([MOORING_NAME,' = ncreadall(DIR(ii).name);']); % custome function to open .nc files easily
    eval([MOORING_NAME,'.TIME_datenum = ',MOORING_NAME,'.TIME + TIME_0;']);
    disp([MOORING_NAME,' loaded']);
end

MOORING_NAMES % unsuppressed for validation purposes

%% Make a pseudo-CTD at 500 meters depth from profiler data in order to close the gap between 500 and 600 meters for steric height

% % % LOAD PROFILERS
% FOLDER_P = '/Users/kachelein/Documents/JPL/work/07.2023_QC_moorings/PROFILERS/';
% FOLDER_P = '/mnt/flow/swot/calval/Data/insitu_mooring/JPL_QC/PROFILERS/';
FOLDER_P = '/mnt/flow/swot/Analysis_Luke/Moorings2/DATA/PROF/';
FOLDER_P = '/Users/kachelein/Documents/JPL/work/02.2024/OSM/DATA/PROF/'; % local
DIR_P = dir([FOLDER_P 'SWOTPOSTLAUNCH_L2_JPLQC_MOORING-S*.nc']);
cd(FOLDER_P)

MOORING_NAMES = cell(length(DIR_P),1);

for ii = 1:length(DIR_P)
    MOORING_NAME = DIR_P(ii).name(33:34);
    MOORING_NAMES{ii} = MOORING_NAME;
    eval([MOORING_NAME,'p = ncreadall(DIR_P(ii).name);']);%%%
    % ^ This loads everything, so the next step elimiates memory-heavy unneeded variables:
    eval([MOORING_NAME 'p = rmfield(' MOORING_NAME 'p, {''CNDC'',''TEMP'',''PRES'',''PSAL'',''DEPTH'',''PROF_NUM'',''QC_FLAG'',''RHO''});']);
    eval([MOORING_NAME 'p.TIME_datenum = ',MOORING_NAME,'p.TIME + TIME_0;']);
    disp([MOORING_NAME 'p (profiler) loaded']);
end
wsm %%% custom fucntion "WorkSpaceMemory" to show how much memory is being used

MOORING_NAMES

%% Variables:

% FIXED:
    %                   CNDC: [56966×7 double]
    %                CNDC_QC: [56966×7 int8]
    %                  DEPTH: [56966×7 double]
    %          DEPTH_NOMINAL: [7×1 double]
    %     INSTR_HAS_P_SENSOR: [7×1 int8]
    %             INSTR_INFO: -127
    %             INSTR_MAKE: [19×1 char]
    %            INSTR_MODEL: [9×7 char]
    %               INSTR_SN: [7×1 int16]
    %               LATITUDE: 36.1850
    %              LONGITUDE: -125.1250
    %                   PRES: [56966×7 double]
    %                PRES_QC: [56966×7 int8]
    %                   PSAL: [56966×7 double]
    %                PSAL_QC: [56966×7 int8]
    %                QC_FLAG: [56966×7 double]
    %                    RHO: [56966×7 double]
    %                   TEMP: [56966×7 double]
    %                TEMP_QC: [56966×7 int8]
    %                   TIME: [56966×1 double]
    % STERIC_HEIGHT_BY_LAYER: [56966×7 double]
    %     RHO_TIMEMEAN_DEPTH: [397×1 double]
    %           RHO_TIMEMEAN: [397×1 double]
    %     STERIC_HEIGHT_MEAN: -10.9816
    %   TIME_PSEUDO_CTD_500m: [8751×1 double]
    %  DEPTH_PSEUDO_CTD_500m: [8751×1 double]
    %    RHO_PSEUDO_CTD_500m: [8751×1 double]

% PROF:
    %                          TIME: [33970165×1 double]
    %                      LATITUDE: 35.9170
    %                     LONGITUDE: -125.0440
    %*                         CNDC: [33970165×1 double]
    %*                         TEMP: [33970165×1 double]
    %*                         PRES: [33970165×1 double]
    %*                         PSAL: [33970165×1 double]
    %*                        DEPTH: [33970165×1 double]
    %*                     PROF_NUM: [33970165×1 double]
    %        STERIC_HEIGHT_ABOVE500: [8499×1 double]
    %           PRESSURE_MIDPROFILE: 247.5000
    %                PRESSURE_RANGE: [2×1 double]
    %               TIME_MIDPROFILE: [8499×1 double]
    %*                      QC_FLAG: [33970165×1 uint8]
    %                    INSTR_INFO: ' '
    %                    INSTR_MAKE: [10×1 char]
    %                   INSTR_MODEL: [10×1 char]
    %                      INSTR_SN: [10×1 char]
    %*                          RHO: [33970165×1 double]
    %         STERIC_HEIGHT_ANOMALY: [8499×1 double]
    %            STERIC_HEIGHT_MEAN: 0.8550
    %                  RHO_TIMEMEAN: [501×1 double]
    %            RHO_TIMEMEAN_DEPTH: [501×1 double]
    %  LATITUDE_GPSSB_STERIC_HEIGHT: [8499×1 double]
    % LONGITUDE_GPSSB_STERIC_HEIGHT: [8499×1 double]
    %            TIME_STERIC_HEIGHT: [8499×1 double]
    %        PROF_GOOD_REALIZATIONS: [8499×1 double]
    %    STERIC_HEIGHT_CUTOFF_DEPTH: 500
    %  STERIC_HEIGHT_MIN_PROF_DEPTH: [8499×1 double]
    %  STERIC_HEIGHT_MAX_PROF_DEPTH: [8499×1 double]
    %                         RHO_0: 1029


%% Times covered:

TIME_0 = datenum('1950-01-01 00:00:00');

% datestr(TIME_0 + S1p.TIME(1))
% datestr(TIME_0 + S2p.TIME(1))
% datestr(TIME_0 + S3p.TIME(1))
% datestr(TIME_0 + S4p.TIME(1))
    % '03-Mar-2023 23:11:30'
    % '03-Mar-2023 01:18:59'
    % '27-Feb-2023 23:18:57'
    % '27-Feb-2023 00:19:56'

% datestr(TIME_0 + S1.TIME(1))
% datestr(TIME_0 + S2.TIME(1))
% datestr(TIME_0 + S3.TIME(1))
% datestr(TIME_0 + S4.TIME(1))
    % '03-Mar-2023 22:50:01'
    % '03-Mar-2023 00:55:01'
    % '27-Feb-2023 23:00:00'
    % '26-Feb-2023 23:55:01'



% datestr(TIME_0 + S1p.TIME(end))
% datestr(TIME_0 + S2p.TIME(end))
% datestr(TIME_0 + S3p.TIME(end))
% datestr(TIME_0 + S4p.TIME(end))
    % '17-Sep-2023 18:24:12'
    % '15-Sep-2023 15:23:41'
    % '16-Sep-2023 13:52:18'
    % '16-Sep-2023 23:47:42'


% datestr(TIME_0 + S1.TIME(end))
% datestr(TIME_0 + S2.TIME(end))
% datestr(TIME_0 + S3.TIME(end))
% datestr(TIME_0 + S4.TIME(end))
    % '17-Sep-2023 17:55:01'
    % '15-Sep-2023 15:10:01'
    % '16-Sep-2023 13:50:00'
    % '16-Sep-2023 23:25:01'

%% Set gridding parameters

% This is strictly a time-domain gridding, but the code I have to implement
% OI is necessarily 2D, so the gridding coordinates will be linear in time
% and constant in "depth" (or more accurately some fake variable).

t_decor = 0.5/24; % units of days
dt_grid =   0.5/24;

% S1.TIME_GRID = [datenum('04-Mar-2023 00:00:00') - TIME_0]:dt_grid:...
%                [datenum('04-Jul-2023 11:00:00') - TIME_0];
% S1.TIME_GRID = S1.TIME_GRID';
% 
% S2.TIME_GRID = [datenum('03-Mar-2023 02:00:00') - TIME_0]:dt_grid:...
%                [datenum('02-Aug-2023 23:30:00') - TIME_0];
% S2.TIME_GRID = S2.TIME_GRID';
% 
% S3.TIME_GRID = [datenum('28-Feb-2023 00:00:00') - TIME_0]:dt_grid:...
%                [datenum('02-Aug-2023 22:00:00') - TIME_0];
% S3.TIME_GRID = S3.TIME_GRID';
% 
% S4.TIME_GRID = [datenum('27-Feb-2023 00:00:00') - TIME_0]:dt_grid:...
%                [datenum('30-Jul-2023 12:00:00') - TIME_0];
% S4.TIME_GRID = S4.TIME_GRID';

S1.TIME_GRID = [datenum('27-Feb-2023 00:00:00') - TIME_0]:dt_grid:...
               [datenum('17-Sep-2023 00:00:00') - TIME_0];
S1.TIME_GRID = S1.TIME_GRID';
S2.TIME_GRID = S1.TIME_GRID;
S3.TIME_GRID = S1.TIME_GRID;
S4.TIME_GRID = S1.TIME_GRID;

%% Perform the optimal interpolation

% S#.STERIC_HEIGHT_ANOMALY_GRID_PROF
% S#.STERIC_HEIGHT_ANOMALY_GRID_FIXED
% S#.STERIC_HEIGHT_ANOMALY_GRID_TOTAL

GOOD_PROF_CUTOFF = 350;
PROF_LENGTH_CUTOFF = 480; % omits about 1% of casts

for ii = 1:length(DIR)
    SS = MOORING_NAMES{ii};

    eval([SS '.STERIC_HEIGHT_ANOMALY_GRID_FIXED = nan(size(' SS '.TIME_GRID));']);
    eval([SS '.STERIC_HEIGHT_ANOMALY_GRID_PROF = nan(size(' SS '.TIME_GRID));']);
    eval([SS '.STERIC_HEIGHT_ANOMALY_GRID_TOTAL = nan(size(' SS '.TIME_GRID));']);

    for ti = 1:length(eval([SS '.TIME_GRID']))
        % % % Grid fixed steric height anomaly (in practice, all matrices in the OI inversion are singular, so instead use bin averaging):
        % IND = eval(['[' SS '.TIME > ' SS '.TIME_GRID(ti)-dt_grid  &  ' SS '.TIME < ' SS '.TIME_GRID(ti)+dt_grid & ~sum(' SS '.QC_FLAG,1)'']']); % 
        % if sum(IND)>0 && sum(isfinite(eval(['sum(' SS '.STERIC_HEIGHT_BY_LAYER(:,IND),2)'])))>0
        %     STERIC_HEIGHT_ANOMALY_GRID = ObjMap(...
        %         eval(['sum(' SS '.STERIC_HEIGHT_BY_LAYER(:,IND),1)'])', ... % data
        %         eval([SS '.TIME(IND)']), ... % x
        %         eval([SS '.TIME(IND)*0']), ... % y, fake variable
        %         t_decor, ... % xdecor
        %         1, ... % ydecor, fake variable
        %         eval([SS '.TIME_GRID(ti)']), ... % xgrid
        %         0, ... % ygrid, fake variable
        %         0, ... % theta
        %         0.1, ... % noise2sig
        %         0); % plotboolean
        % else
        %     STERIC_HEIGHT_ANOMALY_GRID = nan;
        % end
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


        % Grid fixed steric height anomaly
        % IND = eval(['[' SS '.TIME > ' SS '.TIME_GRID(ti)-dt_grid  &  ' SS '.TIME < ' SS '.TIME_GRID(ti)+dt_grid & ~sum(' SS '.QC_FLAG, 2)]']); %
        % % % ^ This was the old way, where the bottom profiler wasn't
        % % % always broken; the next line accounts for this:
        IND = eval(['[' SS '.TIME > ' SS '.TIME_GRID(ti)-dt_grid  &  ' SS '.TIME < ' SS '.TIME_GRID(ti)+dt_grid &   sum(~~' SS '.QC_FLAG, 2)<=1]']); %
        if sum(IND)>0 && sum(isfinite(eval(['sum(' SS '.STERIC_HEIGHT_BY_LAYER(IND, 1:[end-1]), 2 )'])))>0
            STERIC_HEIGHT_ANOMALY_GRID = mean(eval(['sum(' SS '.STERIC_HEIGHT_BY_LAYER(IND, [end-1]), 2 )']),'omitnan');
        else
            STERIC_HEIGHT_ANOMALY_GRID = nan;
        end

        % Grid profiler steric height anomaly:
        INDp = eval(['[' SS 'p.TIME_MIDPROFILE > ' SS '.TIME_GRID(ti)-dt_grid  &  ' ...
                         SS 'p.TIME_MIDPROFILE < ' SS '.TIME_GRID(ti)+dt_grid  &  ' ...
                         SS 'p.PROF_GOOD_REALIZATIONS > ' num2str(GOOD_PROF_CUTOFF) ...
                  ' & [' SS 'p.STERIC_HEIGHT_MAX_PROF_DEPTH - ' SS 'p.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' num2str(PROF_LENGTH_CUTOFF) ']']); %  & ~' SS 'p.QC_FLAG
        if sum(INDp)>0 && sum(isfinite(eval([SS 'p.STERIC_HEIGHT_ANOMALY(INDp)'])))>0
            STERIC_HEIGHT_ANOMALY_GRID_p = mean(eval([SS 'p.STERIC_HEIGHT_ANOMALY(INDp)']),'omitnan');
        else
            disp([SS ': no valid data at ti = ' num2str(ti)])
            STERIC_HEIGHT_ANOMALY_GRID_p = nan;
        end

        % Save variables before finally saving them into the netCDT parent files
        eval([SS '.STERIC_HEIGHT_ANOMALY_GRID_FIXED(ti) = STERIC_HEIGHT_ANOMALY_GRID;']);
        eval([SS '.STERIC_HEIGHT_ANOMALY_GRID_PROF(ti) = STERIC_HEIGHT_ANOMALY_GRID_p;']);
        eval([SS '.STERIC_HEIGHT_ANOMALY_GRID_TOTAL(ti) = STERIC_HEIGHT_ANOMALY_GRID + STERIC_HEIGHT_ANOMALY_GRID_p;']);

        if ~mod(ti,500); disp(['ti = ' num2str(ti) 'out of ' num2str(length(eval([SS '.TIME_GRID'])))]); else; end
    end
    eval([SS '.STERIC_HEIGHT_MEAN_FIXED = ' SS '.STERIC_HEIGHT_MEAN;']);
    eval([SS '.STERIC_HEIGHT_MEAN_PROF = ' SS 'p.STERIC_HEIGHT_MEAN;']);
    eval([SS '.STERIC_HEIGHT_MEAN_TOTAL = ' SS '.STERIC_HEIGHT_MEAN + ' SS 'p.STERIC_HEIGHT_MEAN;']);
end

%% Save into netCDF files

% error('Make sure that this works by testing manually.')

for ii = 1:length(DIR_P)
    SS = MOORING_NAMES{ii};
    % Reminder: DIR is for fixed, DIR_P is for profilers

    % Save time grid to fixed files
    save_to_nc([DIR(ii).folder '/' DIR(ii).name],... file name
        'TIME_GRID',... name of new variable
        eval([SS '.TIME_GRID']),... data
        {'TIME_GRID',length(eval([SS '.TIME_GRID']))},... DIM_METADATA
        {{'long_name',['Time at which steric height is gridded via interval averaging (+-' num2str(dt_grid) ' days).']},...
        {'standard_name','time_grid'},...
        {'units','days since 1950-01-01T00:00:00Z'},...
        {'coordinates','TIME LONGITUDE LATITUDE'}},... METADATA
        true,true);% OVERWRITE, SKIP, MULTIPLE_DIMS

    % Save time grid to profiler files
    save_to_nc([DIR_P(ii).folder '/' DIR_P(ii).name],... file name
        'TIME_GRID',... name of new variable
        eval([SS '.TIME_GRID']),... data
        {'TIME_GRID',length(eval([SS '.TIME_GRID']))},... DIM_METADATA
        {{'long_name',['Time at which steric height is gridded via interval averaging (+-' num2str(dt_grid) ' days).']},...
        {'standard_name','time_grid'},...
        {'units','days since 1950-01-01T00:00:00Z'},...
        {'coordinates','TIME LONGITUDE LATITUDE'}},... METADATA
        true,true);% OVERWRITE, SKIP, MULTIPLE_DIMS

    % save profiler SH to profiler files
    save_to_nc([DIR_P(ii).folder '/' DIR_P(ii).name],... file name
        'STERIC_HEIGHT_ANOMALY_GRID_PROF',... name of new variable
        eval([SS '.STERIC_HEIGHT_ANOMALY_GRID_PROF']),... data
        {'TIME_GRID'},... DIM_METADATA
        {{'long_name',['Steric height from profilers (above 500 m) gridded via interval averaging (+-' num2str(dt_grid) ' days).']},...
        {'standard_name','steric_height_grid_prof'},...
        {'units','meters'},...
        {'coordinates','TIME LONGITUDE LATITUDE'}},... METADATA
        true,true);% OVERWRITE, SKIP, MULTIPLE_DIMS

    % save fixed SH to fixed files
    save_to_nc([DIR(ii).folder '/' DIR(ii).name],... file name
        'STERIC_HEIGHT_ANOMALY_GRID_FIXED',... name of new variable
        eval([SS '.STERIC_HEIGHT_ANOMALY_GRID_FIXED']),... data
        {'TIME_GRID'},... DIM_METADATA
        {{'long_name',['Steric height from fixed CTDs (below 500 m) gridded via interval averaging (+-' num2str(dt_grid) ' days).']},...
        {'standard_name','steric_height_grid_fixed'},...
        {'units','meters'},...
        {'coordinates','TIME LONGITUDE LATITUDE'}},... METADATA
        true,true);% OVERWRITE, SKIP, MULTIPLE_DIMS

    % save total SH to profiler files
    save_to_nc([DIR_P(ii).folder '/' DIR_P(ii).name],... file name
        'STERIC_HEIGHT_ANOMALY_GRID_TOTAL',... name of new variable
        eval([SS '.STERIC_HEIGHT_ANOMALY_GRID_TOTAL']),... data
        {'TIME_GRID'},... DIM_METADATA
        {{'long_name',['Total steric height from fixed CTDs (0 m to the bottom) gridded via interval averaging (+-' num2str(dt_grid) ' days).']},...
        {'standard_name','steric_height_grid_total'},...
        {'units','meters'},...
        {'coordinates','TIME LONGITUDE LATITUDE'}},... METADATA
        true,true);% OVERWRITE, SKIP, MULTIPLE_DIMS

    % save total SH to fixed files
    save_to_nc([DIR(ii).folder '/' DIR(ii).name],... file name
        'STERIC_HEIGHT_ANOMALY_GRID_TOTAL',... name of new variable
        eval([SS '.STERIC_HEIGHT_ANOMALY_GRID_TOTAL']),... data
        {'TIME_GRID'},... DIM_METADATA
        {{'long_name',['Total steric height from fixed CTDs (0 m to the bottom) gridded via interval averaging (+-' num2str(dt_grid) ' days).']},...
        {'standard_name','steric_height_grid_total'},...
        {'units','meters'},...
        {'coordinates','TIME LONGITUDE LATITUDE'}},... METADATA
        true,true);% OVERWRITE, SKIP, MULTIPLE_DIMS

end

% Done seprately in case there is an issue:
for ii = 1:length(DIR_P)
    SS = MOORING_NAMES{ii};
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

% % % FIXED
  % CNDC                                6x32529            1561392  double              
  % DEPTH                               6x32529            1561392  double              
  % DEPTH_NOMINAL                       6x1                     48  double              
  % DEPTH_PSEUDO_CTD_500m            4816x1                  38528  double              
  % INSTR_INFO                          1x1                      2  char                
  % INSTR_MAKE                         20x6                    240  char                
  % INSTR_MODEL                        20x6                    240  char                
  % INSTR_SN                           20x6                    240  char                
  % LATITUDE                            1x1                      8  double              
  % LATITUDE_GPS_SURFACE_BUOY        5888x1                  47104  double              
  % LONGITUDE                           1x1                      8  double              
  % LONGITUDE_GPS_SURFACE_BUOY       5888x1                  47104  double              
  % PRES                                6x32529            1561392  double              
  % PRESSURE_MIDPROFILE                 1x1                      8  double              
  % PRESSURE_RANGE                      2x1                     16  double              
  % PSAL                                6x32529            1561392  double              
  % QC_FLAG                             6x32529            1561392  double              
  % RHO                                 6x32529            1561392  double              
  % RHO_PSEUDO_CTD_500m              4816x1                  38528  double              
  % RHO_TIMEMEAN                      397x1                   3176  double              
  % RHO_TIMEMEAN_DEPTH                397x1                   3176  double              
  % STERIC_HEIGHT_BELOW500          32529x1                 260232  double              
  % STERIC_HEIGHT_BY_LAYER              6x32529            1561392  double              
  % STERIC_HEIGHT_MEAN                  1x1                      8  double              
  % TEMP                                6x32529            1561392  double              
  % TIME                            32529x1                 260232  double              
  % TIME_GPS_SURFACE_BUOY            5888x1                  47104  double              
  % TIME_PSEUDO_CTD_500m             4816x1                  38528  double 


% % % PROFILERS
  % CNDC                               4858562x1             38868496  double              
  % DEPTH                              4858562x1             38868496  double              
  % INSTR_INFO                               1x1                    2  char                
  % INSTR_MAKE                              10x1                   20  char                
  % INSTR_MODEL                             10x1                   20  char                
  % INSTR_SN                                10x1                   20  char                
  % LATITUDE                                 1x1                    8  double              
  % LATITUDE_GPSSB_STERIC_HEIGHT          5589x1                44712  double              
  % LATITUDE_GPS_SURFACE_BUOY             5888x1                47104  double              
  % LONGITUDE                                1x1                    8  double              
  % LONGITUDE_GPSSB_STERIC_HEIGHT         5589x1                44712  double              
  % LONGITUDE_GPS_SURFACE_BUOY            5888x1                47104  double              
  % PRES                               4858562x1             38868496  double              
  % PRESSURE_MIDPROFILE                      1x1                    8  double              
  % PRESSURE_RANGE                           2x1                   16  double              
  % PROF_GOOD_REALIZATIONS                5589x1                44712  double              
  % PROF_NUM                           4858562x1             38868496  double              
  % PSAL                               4858562x1             38868496  double              
  % QC_FLAG                            4858562x1              4858562  uint8               
  % RHO                                4858562x1             38868496  double              
  % RHO_0                                    1x1                    8  double              
  % RHO_TIMEMEAN                           501x1                 4008  double              
  % RHO_TIMEMEAN_DEPTH                     501x1                 4008  double              
  % STERIC_HEIGHT_ABOVE500                5589x1                44712  double              
  % STERIC_HEIGHT_ANOMALY                 5589x1                44712  double              
  % STERIC_HEIGHT_CUTOFF_DEPTH               1x1                    8  double              
  % STERIC_HEIGHT_MAX_PROF_DEPTH          5589x1                44712  double              
  % STERIC_HEIGHT_MEAN                       1x1                    8  double              
  % STERIC_HEIGHT_MIN_PROF_DEPTH          5589x1                44712  double              
  % TEMP                               4858562x1             38868496  double              
  % TIME                               4858562x1             38868496  double              
  % TIME_GPS_SURFACE_BUOY                 5888x1                47104  double              
  % TIME_MIDPROFILE                       5589x1                44712  double              
  % TIME_STERIC_HEIGHT                    5589x1                44712  double 