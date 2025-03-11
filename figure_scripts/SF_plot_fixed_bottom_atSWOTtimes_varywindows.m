%% Figure and numerical results: interpolated to fixed depth and at densely spaced times
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

SWOT_time_boolean = true;  % if you want to look at SWOT times
SWOT_time_boolean = false; % if you want to look at time steps spaced by averaging window

%% SP.600m.##hr

SWOT_time = nan(length(SWOT),1);
for ti = 1:length(SWOT.time)
    SWOT_time(ti) = median(SWOT.time{ti},'omitnan');
end

% HR_intervals_hr = [1 6 12 24 25 50]';
HR_intervals_hr = [1 2 3 6 12 12.42 24 25 50]';

close all

for ii = 1:length(HR_intervals_hr)
    
    SF_TIME_STEP = HR_intervals_hr(ii)/24; % days
    % TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
    % TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
    TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];
    
    DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
    SEP_VEC = [];
    SF_VEC = [];
    % figure% visualize intervals
    tic
    for ti = 1:[SWOT_time_boolean*length(SWOT_time) + ~SWOT_time_boolean*(length(TIME_BIN_EDGES) - 1)] % ti = 1:[length(TIME_BIN_EDGES) - 1]
        DATA_in_interval = [];
        XY_in_interval = [];
        if SWOT_time_boolean
            SF_TIME_INTERVAL = [SWOT_time(ti) + SF_TIME_STEP*[-0.5 0.5]] - T0;
        else
            SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
        end
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
            % disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
            disp([num2str(ti) '/' num2str(length(SWOT_time)) ' - ' replace(num2str(HR_intervals_hr(ii)),'.','_') ' hr'])
        end
        % plot(ti*[1 1], SF_TIME_INTERVAL, '.-k');hold on% visualize intervals
    end
    toc
    
    eval(['Separation_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr = SEP_VEC;'])
    eval(['StructFunc_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr = SF_VEC;'])
    
    clear SF_VEC SEP_VEC

end

%% ! SP.600m.1hr.NoP1

% close all
% 
% OMITTED_MOORING = 'P1';
% 
% SF_TIME_STEP = 1/24; % days
% % TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
% TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
% TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];
% 
% DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
% SEP_VEC = [];
% SF_VEC = [];
% 
% tic
% for ti = 1:[SWOT_time_boolean*length(SWOT_time) + ~SWOT_time_boolean*(length(TIME_BIN_EDGES) - 1)] % ti = 1:[length(TIME_BIN_EDGES) - 1]
%     DATA_in_interval = [];
%     XY_in_interval = [];
%         if SWOT_time_boolean
%             SF_TIME_INTERVAL = [SWOT_time(ti) + SF_TIME_STEP*[-0.5 0.5]] - T0;
%         else
%             SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
%         end
%     for mi = 1:length(MM_list)
%         MM = MM_list{mi};
%         if strcmp(MM(1),'P') & ~[strcmp(MM,OMITTED_MOORING)]
%             % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
%             %           interp1_custom( T0 + eval([MM '.F.TIME']), ...
%             %                  eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
%             %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])), ...
%             %           '.:', 'MarkerSize',10);
% 
%             P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
%                        '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
%                        '~SPIKE_IND{mi}']); % <-- Spike removal
%             DATA_in_interval = [DATA_in_interval ; ...
%                  mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
%             % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
%             %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
%             XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
%              mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
%         elseif strcmp(MM(1),'S')
%             % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']), ...
%             %      interp1_custom( T0 + eval([MM '.F.TIME']), ...
%             %                  eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
%             %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])) , ...
%             %      '.-')
% 
%             P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
%                         '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
%                         '~SPIKE_IND{mi}']); % <-- Spike removal
%             DATA_in_interval = [DATA_in_interval ; ...
%                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
%             % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
%             %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
%             XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
%                 mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
%         else
%         end
%     end
%     % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
%     SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
%                               [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
%     SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;
% 
%     SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
%     SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];
% 
%     if ~mod(ti,10)
%         disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
%     end
% end
% toc
% 
% % To reiterate, the results from this section are for S an P moorings and
% % combining profilers with the upper-most fixed CTD (~600 m).
% Separation_S_P_600m_1hr_noP1 = SEP_VEC;
% StructFunc_S_P_600m_1hr_noP1 = SF_VEC;
% 
% close all
% figure
% % % % SP.600m.1hr
% LEG_TEXT = {}; LEG_TEXT = {LEG_TEXT{:},'Moorings (600m, 1hr)'};
% BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
% BIN_CENTR = interval_avg(Separation_S_P_600m_1hr,Separation_S_P_600m_1hr,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
% BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
%                [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
%                [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% % plot(Separation_S_P_600m_1hr, StructFunc_S_P_600m_1hr, '.k', 'MarkerSize',10); hold on
% [Mean_SF_Separation_S_P_600m_1hr] = ...
%     errorbar_stdofmean(BIN_CENTR, Separation_S_P_600m_1hr, m2_to_cm2*StructFunc_S_P_600m_1hr, BIN_EDGES, CO(1,:)); hold on
% BIN_CENTR_S_P = BIN_CENTR;
% % % % SP.600m.1hr.noP1
% LEG_TEXT = {LEG_TEXT{:},['Moorings, no ' OMITTED_MOORING ' (600m, 1hr)']};
% BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
% BIN_CENTR = interval_avg(Separation_S_P_600m_1hr_noP1,Separation_S_P_600m_1hr_noP1,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
% BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
%                [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
%                [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% % plot(Separation_S_P_600m_1hr_noP1, StructFunc_S_P_600m_1hr_noP1, '.k', 'MarkerSize',10); hold on
% [Mean_SF_Separation_S_P_600m_1hr_noP1] = ...
%     errorbar_stdofmean(BIN_CENTR, Separation_S_P_600m_1hr_noP1, m2_to_cm2*StructFunc_S_P_600m_1hr_noP1, BIN_EDGES, CO(1,:), '--'); hold on
% LEG = legend(LEG_TEXT,'Location','southeast');
% LEG_TEXT = {};
% 
% 
% 
%% ! SP.600m.24hr.NoP1
% 
% close all
% 
% SF_TIME_STEP = 24/24; % days
% % TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
% TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
% TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];
% 
% DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
% SEP_VEC = [];
% SF_VEC = [];
% 
% tic
% for ti = 1:[SWOT_time_boolean*length(SWOT_time) + ~SWOT_time_boolean*(length(TIME_BIN_EDGES) - 1)] % ti = 1:[length(TIME_BIN_EDGES) - 1]
%     DATA_in_interval = [];
%     XY_in_interval = [];
%     if SWOT_time_boolean
%         SF_TIME_INTERVAL = [SWOT_time(ti) + SF_TIME_STEP*[-0.5 0.5]] - T0;
%     else
%         SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
%     end
%     for mi = 1:length(MM_list)
%         MM = MM_list{mi};
%         if strcmp(MM(1),'P') & ~[strcmp(MM,'P1')]
%             % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']) , ...
%             %           interp1_custom( T0 + eval([MM '.F.TIME']), ...
%             %                  eval(['sum(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED.*[~' MM '.F.QC_FLAG]./[~' MM '.F.QC_FLAG] , 2)']) , ..., ''omitnan''
%             %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])), ...
%             %           '.:', 'MarkerSize',10);
% 
%             P_IND = eval([ MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
%                        '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
%                        '~SPIKE_IND{mi}']); % <-- Spike removal
%             DATA_in_interval = [DATA_in_interval ; ...
%                  mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
%             % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
%             %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
%             XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
%              mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
%         elseif strcmp(MM(1),'S')
%             % plot(T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})']), ...
%             %      interp1_custom( T0 + eval([MM '.F.TIME']), ...
%             %                  eval(['(' MM '.F.STERIC_HEIGHT_BY_LAYER_AUGMENTED(:,1).*[~' MM '.F.QC_FLAG(:,1)]./[~' MM '.F.QC_FLAG(:,1)] )']) , ...
%             %                  T0 + eval([MM '.P.TIME_STERIC_HEIGHT(' MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ~SPIKE_IND{mi})'])) , ...
%             %      '.-')
% 
%             P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
%                         '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
%                         '~SPIKE_IND{mi}']); % <-- Spike removal
%             DATA_in_interval = [DATA_in_interval ; ...
%                 mean(eval([MM '_sh_600(P_IND)']), 'omitnan') ];
%             % XY_in_interval = [XY_in_interval ; ... NOMINAL LOCATION FOR NOW
%             %           eval([MM '.P.LONGITUDE*ones(sum(P_IND),1) + 1i*' MM '.P.LATITUDE*ones(sum(P_IND),1)']) ];
%             XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
%                 mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
%         else
%         end
%     end
%     % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
%     SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
%                               [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
%     SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;
% 
%     SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
%     SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];
% 
%     if ~mod(ti,10)
%         disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
%     end
% end
% toc
% 
% % To reiterate, the results from this section are for S an P moorings and
% % combining profilers with the upper-most fixed CTD (~600 m).
% Separation_S_P_600m_24hr_noP1 = SEP_VEC;
% StructFunc_S_P_600m_24hr_noP1 = SF_VEC;

%% ! S.600m.1hr
% allps=interval_avg(Separation_S_P_600m_24hr,StructFunc_S_P_600m_24hr,BIN_EDGES);


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
for ti = 1:length(SWOT.time)
    % if isfinite(median(SWOT.time{swot_ti},'omitnan'))
    %     ti = dsearchn(TIME_BIN_CENTERS',median(SWOT.time{swot_ti},'omitnan'));
        DATA_in_interval = [];
        XY_in_interval = [];
        if SWOT_time_boolean
            SF_TIME_INTERVAL = [SWOT_time(ti) + SF_TIME_STEP*[-0.5 0.5]] - T0;
        else
            SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
        end
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
    % else
    % end
end
toc

Separation_S_600m_1hr = SEP_VEC;
StructFunc_S_600m_1hr = SF_VEC;

clear SF_VEC SEP_VEC




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
for ti = 1:[SWOT_time_boolean*length(SWOT_time) + ~SWOT_time_boolean*(length(TIME_BIN_EDGES) - 1)] % ti = 1:[length(TIME_BIN_EDGES) - 1]
    DATA_in_interval = [];
    XY_in_interval = [];
        if SWOT_time_boolean
            SF_TIME_INTERVAL = [SWOT_time(ti) + SF_TIME_STEP*[-0.5 0.5]] - T0;
        else
            SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
        end
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

%% ! SPG.600m.24hr - Moorings AND Gliders

close all

% GM_offset = 0;

SF_TIME_STEP = 24/24; % days
TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];

DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
SEP_VEC = [];
SF_VEC = [];

MM_GG_list = {'P1','P2','P3',    'P5','P6','P7',    'S1','S2','S3','S4',     'ru32_600','ru38_600'};
% MM_GG_list = {'ru32_600','ru38_600'}; % for testing purposes

tic
for ti = 1:[SWOT_time_boolean*length(SWOT_time) + ~SWOT_time_boolean*(length(TIME_BIN_EDGES) - 1)] % ti = 1:[length(TIME_BIN_EDGES) - 1]
    DATA_in_interval = [];
    XY_in_interval = [];
        if SWOT_time_boolean
            SF_TIME_INTERVAL = [SWOT_time(ti) + SF_TIME_STEP*[-0.5 0.5]] - T0;
        else
            SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
        end
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

%% S.1800m.24hr

% close all
% 
% SF_TIME_STEP = 24/24; % days
% TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
% TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
% TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];
% 
% DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
% SEP_VEC = [];
% SF_VEC = [];
% 
% % S_MM_list = {'S1','S2','S3','S4'};
% 
% tic
% for ti = 1:[SWOT_time_boolean*length(SWOT_time) + ~SWOT_time_boolean*(length(TIME_BIN_EDGES) - 1)] % ti = 1:[length(TIME_BIN_EDGES) - 1]
%     DATA_in_interval = [];
%     XY_in_interval = [];
%     if SWOT_time_boolean
%         SF_TIME_INTERVAL = [SWOT_time(ti) + SF_TIME_STEP*[-0.5 0.5]] - T0;
%     else
%         SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
%     end
%     for mi = 1:length(MM_list)
%         MM = MM_list{mi};
%         if strcmp(MM(1),'P')
%             % % % P-MOORINGS OMITTED FROM THIS ANALYSIS
%         elseif strcmp(MM(1),'S')
%             P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
%                         '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
%                         '~SPIKE_IND{mi}']); % <-- Spike removal
%             DATA_in_interval = [DATA_in_interval ; ...
%                 mean(eval([MM '_sh_1800(P_IND)']), 'omitnan') ];
%             XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
%                   mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
%         else
%         end
%     end
%     % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
%     SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
%                               [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
%     SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;
% 
%     SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
%     SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];
% 
%     if ~mod(ti,10)
%         disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
%     end
% end
% toc
% 
% Separation_S_deep_24hr = SEP_VEC;
% StructFunc_S_deep_24hr = SF_VEC;

%% S.1800m.1hr

% % S1.F.DEPTH_NOMINAL'
% % S2.F.DEPTH_NOMINAL'
% % S3.F.DEPTH_NOMINAL'
% % S4.F.DEPTH_NOMINAL'
% %          605         804        1205        1814        2523        3432        4427
% %          604         803        1204        1706        2361        3268        4373
% %          605         804        1205        1814        2523        3432        4428
% %          608         798        1199        1813        2522        3432        4327
% i_deepest_F = 4;
% 
% close all
% 
% SF_TIME_STEP = 1/24; % days
% TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
% TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
% TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];
% 
% DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
% SEP_VEC = [];
% SF_VEC = [];
% 
% % S_MM_list = {'S1','S2','S3','S4'};
% 
% tic
% for ti = 1:[SWOT_time_boolean*length(SWOT_time) + ~SWOT_time_boolean*(length(TIME_BIN_EDGES) - 1)] % ti = 1:[length(TIME_BIN_EDGES) - 1]
%     DATA_in_interval = [];
%     XY_in_interval = [];
%     if SWOT_time_boolean
%         SF_TIME_INTERVAL = [SWOT_time(ti) + SF_TIME_STEP*[-0.5 0.5]] - T0;
%     else
%         SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
%     end
%     for mi = 1:length(MM_list)
%         MM = MM_list{mi};
%         if strcmp(MM(1),'P')
%             % % % P-MOORINGS OMITTED FROM THIS ANALYSIS
%         elseif strcmp(MM(1),'S')
%             P_IND = eval([  MM '.P.TIME_STERIC_HEIGHT > SF_TIME_INTERVAL(1) & ' MM '.P.TIME_STERIC_HEIGHT < SF_TIME_INTERVAL(2) & ' ... <-- Time interval
%                         '[' MM '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
%                         '~SPIKE_IND{mi}']); % <-- Spike removal
%             DATA_in_interval = [DATA_in_interval ; ...
%                 mean(eval([MM '_sh_1800(P_IND)']), 'omitnan') ];
%             XY_in_interval = [XY_in_interval ; ... GPS LOCATIONS
%                   mean(eval([MM '.P.LONGITUDE_GPSSB_STERIC_HEIGHT(P_IND) + 1i*' MM '.P.LATITUDE_GPSSB_STERIC_HEIGHT(P_IND)']), 'omitnan') ];
%         else
%         end
%     end
%     % SEP_MAT_in_interval = abs([XY_in_interval - XY_in_interval.']); % Naive cartesian approximation
%     SEP_MAT_in_interval = DD*sqrt(imag(XY_in_interval - XY_in_interval.').^2 + ...
%                               [cosd(mean(imag(XY_in_interval),'omitnan')).^2]*real(XY_in_interval - XY_in_interval.').^2); % better approximation (tangent to sphere)
%     SF_MAT_in_interval  = [DATA_in_interval - DATA_in_interval.'].^2;
% 
%     SEP_VEC = [SEP_VEC ; SEP_MAT_in_interval(:)];
%     SF_VEC  = [SF_VEC  ; SF_MAT_in_interval(:) ];
% 
%     if ~mod(ti,10)
%         disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
%     end
% end
% toc
% 
% Separation_S_deep_1hr = SEP_VEC;
% StructFunc_S_deep_1hr = SF_VEC;

%% S.1800m.##hr

% HR_intervals_hr = [1 6 12 24 25 50]';

close all

for ii = 1:length(HR_intervals_hr)
    
    SF_TIME_STEP = HR_intervals_hr(ii)/24; % days
    % TIME_BIN_EDGES = [datenum('2023-05-11 00:00:00'):SF_TIME_STEP:datenum('2023-05-11 23:00:00')];
    % TIME_BIN_EDGES = [datenum('2023-02-24 00:00:00'):SF_TIME_STEP:datenum('2023-09-15 23:00:00')];
    TIME_BIN_EDGES = [datenum('2023-04-01 00:00:00'):SF_TIME_STEP:datenum('2023-07-10 23:00:00')];
    
    DD = 111.3195; % from m_lldist at any fixed longitude, and a 1deg latitude shift
    SEP_VEC = [];
    SF_VEC = [];
    % figure% visualize intervals
    tic
    for ti = 1:[SWOT_time_boolean*length(SWOT_time) + ~SWOT_time_boolean*(length(TIME_BIN_EDGES) - 1)] % ti = 1:[length(TIME_BIN_EDGES) - 1]
        DATA_in_interval = [];
        XY_in_interval = [];
        if SWOT_time_boolean
            SF_TIME_INTERVAL = [SWOT_time(ti) + SF_TIME_STEP*[-0.5 0.5]] - T0;
        else
            SF_TIME_INTERVAL = [TIME_BIN_EDGES(ti) TIME_BIN_EDGES(ti+1)] - T0;
        end
        for mi = 1:length(MM_list)
            MM = MM_list{mi};
            if strcmp(MM(1),'P')
                % Do nothing (no 1800m CTD's on P moorings)
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
            % disp([num2str(ti) '/' num2str(length(TIME_BIN_EDGES) - 1)])
            disp([num2str(ti) '/' num2str(length(SWOT_time)) ' - ' replace(num2str(HR_intervals_hr(ii)),'.','_') ' hr'])
        end
        % plot(ti*[1 1], SF_TIME_INTERVAL, '.-k');hold on% visualize intervals
    end
    toc
    
    eval(['Separation_S_1800m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr = SEP_VEC;'])
    eval(['StructFunc_S_1800m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr = SF_VEC;'])
    
    clear SF_VEC SEP_VEC

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Updated Publication Figure %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1-PANEL SF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
figure('Color',[1 1 1])
dBIN = 0.1;
LEG_TEXT = {};
CO = colororder;
m2_to_cm2 = 100*100; % =10^4 if you are to plot cm^2; =1 if m^2
XSCALE = 'lin'; YSCALE = 'lin';

% tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');
% 
% nexttile;%subplot(121) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


% % % % SP.600m.1hr.noP1
% LEG_TEXT = {LEG_TEXT{:},'Moorings, no P1 (600m, 1hr)'};
% BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
% BIN_CENTR = interval_avg(Separation_S_P_600m_1hr_noP1,Separation_S_P_600m_1hr_noP1,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
% BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
%                [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
%                [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% % plot(Separation_S_P_600m_1hr_noP1, StructFunc_S_P_600m_1hr_noP1, '.k', 'MarkerSize',10); hold on
% [Mean_SF_Separation_S_P_600m_1hr_noP1] = ...
%     errorbar_stdofmean(BIN_CENTR, Separation_S_P_600m_1hr_noP1, m2_to_cm2*StructFunc_S_P_600m_1hr_noP1, BIN_EDGES, CO(1,:), '--'); hold on
% 
% 
% % % SP.600m.24hr.noP1
% LEG_TEXT = {LEG_TEXT{:},'Moorings, no P1 (600m, 24hr)'};
% BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
% BIN_CENTR = interval_avg(Separation_S_P_600m_24hr_noP1,Separation_S_P_600m_24hr_noP1,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
% BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
%                [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
%                [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% % plot(Separation_S_P_600m_24hr_noP1, StructFunc_S_P_600m_24hr_noP1, 'k.', 'MarkerSize',10); hold on
% [Mean_SF_Separation_S_P_600m_24hr_noP1] = ...
%     errorbar_stdofmean(BIN_CENTR, Separation_S_P_600m_24hr_noP1, m2_to_cm2*StructFunc_S_P_600m_24hr_noP1, BIN_EDGES, CO(2,:), '--');

set(gca,'FontSize',18); grid on
set(gca,'yscale',YSCALE,'xscale',XSCALE)

xlabel('Separation (km)')
ylabel('$$\overline{\big[\eta''(x) - \eta''(x+s)\big]^2}$$ \quad (cm$^2$)','Interpreter','latex','FontSize',16)

set(gca,'YLim',[0 50],'YScale','lin')
set(gca,'XLim',[0 103],'XScale','lin')
set(gca,'FontSize',20); grid on
LEG = legend(LEG_TEXT,'Location','southeast');
LEG.Position = [0.1404    0.6490    0.4699    0.2549];
LEG_TEXT = {};

xticks(0:10:100)
yticks(0:5:50)

set(gcf,'Position',[-1439         419         680         392])

INPUT = input(['Do you want to save this figure? Enter any number for "yes" or\n push ' ...
    '"Enter" with a blank input for "No".\n']);

if isempty(INPUT)
    disp(['Image not saved'])
else
    figure(1)
    exportgraphics(gcf,...
        '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/F3.pdf',...
        'BackgroundColor','none','ContentType','vector')
    disp(['Image saved'])
end

%% Figure showing different time windows for SF averaging:

close all

figure('Color','w')
CO = turbo(length(HR_intervals_hr));

for ii = 1:length(HR_intervals_hr)
    LEG_TEXT = {LEG_TEXT{:},['Moorings (600m, ' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr)']};
    BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
    BIN_CENTR = interval_avg(eval(['Separation_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']), ...
                             eval(['Separation_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']),BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
    BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
        [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
        [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
    [Mean_SF_Separation_S_P_600m_NNhr] = ...
        errorbar_stdofmean(BIN_CENTR, ...
                           eval(['Separation_S_P_600m_'           replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']),...
                           eval(['m2_to_cm2*StructFunc_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']), ...
                           BIN_EDGES, CO(ii,:));
    eval(['Mean_SF_Separation_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr = Mean_SF_Separation_S_P_600m_NNhr;']);
    set(gca,'FontSize',18); grid on; hold on
end

set(gca,'FontSize',18); grid on
set(gca,'yscale',YSCALE,'xscale',XSCALE)

xlabel('Separation (km)')
ylabel('$$\overline{\big[\eta''(x) - \eta''(x+s)\big]^2}$$ \quad (cm$^2$)','Interpreter','latex','FontSize',16)

set(gca,'YLim',[0 50],'YScale','lin')
set(gca,'XLim',[0 103],'XScale','lin')
set(gca,'FontSize',20); grid on
LEG = legend(LEG_TEXT,'Location','southeast');
LEG.Position = [0.1404    0.6490    0.4699    0.2549];
LEG_TEXT = {};

xticks(0:10:100)
yticks(0:5:50)

set(gcf,'Position',[-1439         419         680         392])

% INPUT = input(['Do you want to save this figure? Enter any number for "yes" or\n push ' ...
%     '"Enter" with a blank input for "No".\n']);

%% Calculate slopes of SFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% for ii = 1:length(HR_intervals_hr) % including the 50km bump
%     LEG_TEXT = {LEG_TEXT{:},['Moorings (600m, ' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr)']};
%     BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
%     BIN_CENTR = interval_avg(eval(['Separation_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']),...
%                              eval(['Separation_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']),BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
%     BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
%                    [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
%                    [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
%     SF_avg = interval_avg(eval(['Separation_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']),...
%                           eval(['m2_to_cm2*StructFunc_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']), BIN_EDGES);
%     [NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
%                 nonlinear_lsqf(SF_avg(2:[end-1]),BIN_CENTR(2:[end-1]),...
%                 FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
%     NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:[end-1])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
%     disp(['600m/' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr/M to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
%     [NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
%                 nonlinear_lsqf(SF_avg(2:5),BIN_CENTR(2:5),...
%                 FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
%     NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:5)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
%     disp(['600m/' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr/M to ' num2str(BIN_CENTR(5)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
%     disp(' ')
% end

for ii = 1:length(HR_intervals_hr) % Omitting the 50km bump in the 600m SF
    LEG_TEXT = {LEG_TEXT{:},['Moorings (600m, ' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr)']};
    BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
    BIN_CENTR = interval_avg(eval(['Separation_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']),...
                             eval(['Separation_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']),BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
    BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
                   [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
                   [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
    SF_avg = interval_avg(eval(['Separation_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']),...
                          eval(['m2_to_cm2*StructFunc_S_P_600m_' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr']), BIN_EDGES);
    [NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
                nonlinear_lsqf(SF_avg([2:5,7:[end-1]]),BIN_CENTR([2:5,7:[end-1]]),...
                FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
    NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg([2:5,7:[end-1]])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
    disp(['600m/' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr/M to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
    [NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
                nonlinear_lsqf(SF_avg(2:5),BIN_CENTR(2:5),...
                FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
    NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:5)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
    disp(['600m/' replace(num2str(HR_intervals_hr(ii)),'.','_') 'hr/M to ' num2str(BIN_CENTR(5)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
    disp(' ')
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

set(gca,'YLim',[0 7])
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

