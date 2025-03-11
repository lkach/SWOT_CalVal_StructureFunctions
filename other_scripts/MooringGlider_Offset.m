%% Supplemental Material
%% Load data from all gliders and some data from all moorings

% Note: this only compares profiles to 500 meters depth; assume that the
% offset is due to this layer
G_DIR = dir('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/GLIDERS/*Processed_2023*.nc');
M_DIR = dir('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/PROF/SWOT*.nc');

T0 = datenum('1950-01-01 00:00:00');

% % % For now, use these files for GPS mooring data:
M_GPS_DIR = dir('/Users/kachelein/Documents/JPL/work/07.2023_QC_moorings/PROFILERS/SWOT*.nc');
warning('Before publication, be sure to package the GPS data with the other data (or ideally SIO will send us an updated version).')

% % % Variables from source:
% LATITUDE_GPS_SURFACE_BUOY,    LONGITUDE_GPS_SURFACE_BUOY,    TIME_GPS_SURFACE_BUOY
% % % Variables by us (interpolated):
% LATITUDE_GPSSB_STERIC_HEIGHT,    LONGITUDE_GPSSB_STERIC_HEIGHT

% LOAD GLIDERS
cd(G_DIR(1).folder)
GG_list = {};
for ii=1:length(G_DIR)
    GG = G_DIR(ii).name(1:5);
    GG_list{ii} = GG;
    eval([GG '.LONGITUDE_PROFILE = ncread(''' G_DIR(ii).name ''',''LONGITUDE_PROFILE'');']);
    eval([GG '.LATITUDE_PROFILE = ncread(''' G_DIR(ii).name ''',''LATITUDE_PROFILE'');']);
    eval([GG '.TIME_STERIC_HEIGHT = ncread(''' G_DIR(ii).name ''',''TIME_STERIC_HEIGHT'');']);
    eval([GG '.STERIC_HEIGHT_ANOMALY = ncread(''' G_DIR(ii).name ''',''STERIC_HEIGHT_ANOMALY'');']);
    eval([GG '.STERIC_HEIGHT_MAX_PROF_DEPTH = ncread(''' G_DIR(ii).name ''',''STERIC_HEIGHT_MAX_PROF_DEPTH'');']);
    eval([GG '.STERIC_HEIGHT_MIN_PROF_DEPTH = ncread(''' G_DIR(ii).name ''',''STERIC_HEIGHT_MIN_PROF_DEPTH'');']);

    eval([GG '.TIME = ncread(''' G_DIR(ii).name ''',''TIME'');']);
    eval([GG '.DEPTH = single(ncread(''' G_DIR(ii).name ''',''DEPTH''));']);
    eval([GG '.PROF_NUM = int16(ncread(''' G_DIR(ii).name ''',''PROFILE_NUM''));']);
    eval([GG '.QC_FLAG = logical(ncread(''' G_DIR(ii).name ''',''QC_FLAG''));']);
    % eval([MM '.TIME_MIDPROFILE = ncread(''' M_DIR(ii).name ''',''TIME_MIDPROFILE'');']); % provided only for S moorings
end

%% LOAD MOORINGS
cd(M_DIR(1).folder)
MM_list = {};
for ii=1:length(M_DIR)
    MM = M_DIR(ii).name(33:34);
    MM_list{ii} = MM;
    eval([MM '.LATITUDE = ncread(''' M_DIR(ii).name ''',''LATITUDE'');']);
    eval([MM '.LONGITUDE = ncread(''' M_DIR(ii).name ''',''LONGITUDE'');']);

    % eval([MM '.LATITUDE_GPSSB_STERIC_HEIGHT = ncread(''' M_DIR(ii).name ''',''LATITUDE_GPSSB_STERIC_HEIGHT'');']);
    % eval([MM '.LONGITUDE_GPSSB_STERIC_HEIGHT = ncread(''' M_DIR(ii).name ''',''LONGITUDE_GPSSB_STERIC_HEIGHT'');']);

    eval([MM '.TIME_STERIC_HEIGHT = ncread(''' M_DIR(ii).name ''',''TIME_STERIC_HEIGHT'');']); % JPL added
    eval([MM '.STERIC_HEIGHT_ANOMALY = ncread(''' M_DIR(ii).name ''',''STERIC_HEIGHT_ANOMALY'');']); % JPL added
    % eval([MM '.TIME_MIDPROFILE = ncread(''' M_DIR(ii).name ''',''TIME_MIDPROFILE'');']); % provided only for S moorings

    eval([MM '.STERIC_HEIGHT_MIN_PROF_DEPTH = ncread(''' M_DIR(ii).name ''',''STERIC_HEIGHT_MIN_PROF_DEPTH'');']); % JPL added
    eval([MM '.STERIC_HEIGHT_MAX_PROF_DEPTH = ncread(''' M_DIR(ii).name ''',''STERIC_HEIGHT_MAX_PROF_DEPTH'');']); % JPL added

    % These larger variables are used for profile examination
    eval([MM '.TIME = ncread(''' M_DIR(ii).name ''',''TIME'');']);
    eval([MM '.DEPTH = single(ncread(''' M_DIR(ii).name ''',''DEPTH''));']);
    eval([MM '.PROF_NUM = int16(ncread(''' M_DIR(ii).name ''',''PROF_NUM''));']);
    eval([MM '.QC_FLAG = logical(ncread(''' M_DIR(ii).name ''',''QC_FLAG''));']);

    % GPS variables as given by source (we will interpolate to steric height time here on our own):
    eval([MM '.GPS.TIME_GPS_SURFACE_BUOY = ncread(''' M_GPS_DIR(ii).folder '/' M_GPS_DIR(ii).name ''',''TIME_GPS_SURFACE_BUOY'');']);
    eval([MM '.GPS.LATITUDE_GPS_SURFACE_BUOY = ncread(''' M_GPS_DIR(ii).folder '/' M_GPS_DIR(ii).name ''',''LATITUDE_GPS_SURFACE_BUOY'');']);
    eval([MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY = ncread(''' M_GPS_DIR(ii).folder '/' M_GPS_DIR(ii).name ''',''LONGITUDE_GPS_SURFACE_BUOY'');']);
end

%% Interpolate GPS locations to steric height times:
for ii=1:length(M_DIR)
    MM = MM_list{ii};
    % eval(['intlat = interp1(timegps,latgps,timestericheight,''linear'');']);

    eval([MM '.LATITUDE_GPSSB_STERIC_HEIGHT =  interp1(' MM '.GPS.TIME_GPS_SURFACE_BUOY, ' ...
                                                         MM '.GPS.LATITUDE_GPS_SURFACE_BUOY,  ' ...
                                                         MM '.TIME_STERIC_HEIGHT,''linear'');']);
    eval([MM '.LONGITUDE_GPSSB_STERIC_HEIGHT = interp1(' MM '.GPS.TIME_GPS_SURFACE_BUOY, ' ...
                                                         MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY, ' ...
                                                         MM '.TIME_STERIC_HEIGHT,''linear'');']);

    % eval(['MM_LATITUDE_GPS_GLIDERTIME  = interp1(' MM '.GPS.TIME_GPS_SURFACE_BUOY, ' ...
    %                                                MM '.GPS.LATITUDE_GPS_SURFACE_BUOY,  ' ...
    %                                                GG '.TIME_STERIC_HEIGHT,''linear'');']);
    % eval(['MM_LONGITUDE_GPS_GLIDERTIME = interp1(' MM '.GPS.TIME_GPS_SURFACE_BUOY, ' ...
    %                                                MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY, ' ...
    %                                                GG '.TIME_STERIC_HEIGHT,''linear'');']);
end

%% Calculate offsets in time

TOffset_g = cell(length(GG_list),1);
TOffset_m = cell(length(MM_list),1);
for gi = 1:length(GG_list)
    GG = GG_list{gi};
    GG_time = eval([GG '.TIME_STERIC_HEIGHT']);
    TOffset_g{gi} = nan(length(GG_time),length(MM_list));
    for mi = 1:length(MM_list)
        MM = MM_list{mi};
        MM_time = eval([MM '.TIME_STERIC_HEIGHT']);
        for ti = 1:length(GG_time)
            TOffset_ti = abs(GG_time(ti) - MM_time);
            TOffset_g{gi}(ti,mi) = min(TOffset_ti); % days
        end
        TOffset_m_ti = nan(size(MM_time));
        for ti = 1:length(MM_time)
            TOffset_ti = abs(GG_time - MM_time(ti));
            TOffset_m_ti(ti) = min(TOffset_ti); % days
        end
        TOffset_m{mi}(:,gi) = TOffset_m_ti;
    end
end

%% Go through all times and check how far it is from each mooring,
% noting when it is within 1 hr and 2 km of a mooring and note which
% mooring that is. Just check steric height at this point, maybe do T and S
% or RHO later.

% For GPS locations:
% TBD
% Dist_g = [];
% Dist_m = {};
Dist_g = cell(length(GG_list),1);
for gi = 1:length(GG_list)
    GG = GG_list{gi};
    for mi = 1:length(MM_list)
        MM = MM_list{mi};

        eval(['MM_LATITUDE_GPS_GLIDERTIME  = interp1(' MM '.GPS.TIME_GPS_SURFACE_BUOY, ' ...
                                                       MM '.GPS.LATITUDE_GPS_SURFACE_BUOY,  ' ...
                                                       GG '.TIME_STERIC_HEIGHT,''linear'');']);
        eval(['MM_LONGITUDE_GPS_GLIDERTIME = interp1(' MM '.GPS.TIME_GPS_SURFACE_BUOY, ' ...
                                                       MM '.GPS.LONGITUDE_GPS_SURFACE_BUOY, ' ...
                                                       GG '.TIME_STERIC_HEIGHT,''linear'');']);

        MM_lon_vec = MM_LONGITUDE_GPS_GLIDERTIME;
        LON_zippered = [MM_lon_vec' ; eval([GG '.LONGITUDE_PROFILE'])']; LON_zippered = LON_zippered(:);

        MM_lat_vec = MM_LATITUDE_GPS_GLIDERTIME;
        LAT_zippered = [MM_lat_vec' ; eval([GG '.LATITUDE_PROFILE'])']; LAT_zippered = LAT_zippered(:);

        Dist_mi = m_lldist(LON_zippered, LAT_zippered); % M_MAP package required
        Dist_g{gi}(:,mi) = Dist_mi(1:2:[end]); % kilometers
    end
end
Dist_m = cell(length(MM_list),1);
for mi = 1:length(MM_list) % pre-allocate
    MM = MM_list{mi};
    MM_time = eval([MM '.TIME_STERIC_HEIGHT']);
    for gi = 1:length(GG_list)
        GG = GG_list{gi};
        Dist_m{mi} = nan(length(MM_time),length(GG_list));
    end
end

for gi = 1:length(GG_list)
    tic
    GG = GG_list{gi};
    for mi = 1:length(MM_list)
        MM = MM_list{mi};
        MM_time = eval([MM '.TIME_STERIC_HEIGHT']);
        Dist_mi = [];
        for ti = 1:length(MM_time)
            if ~isfinite(MM_time(ti))
                Dist_mi(ti) = NaN;
            else
                Dist_mi(ti) = m_lldist([eval([MM '.LONGITUDE_GPSSB_STERIC_HEIGHT(ti)']), ...
                                            eval([GG '.LONGITUDE_PROFILE(abs(' GG '.TIME_STERIC_HEIGHT-MM_time(ti)) == min(abs(' GG '.TIME_STERIC_HEIGHT-MM_time(ti))))'])],...
                                       [eval([MM '.LATITUDE_GPSSB_STERIC_HEIGHT(ti)'] ), ...
                                            eval([GG '.LATITUDE_PROFILE( abs(' GG '.TIME_STERIC_HEIGHT-MM_time(ti)) == min(abs(' GG '.TIME_STERIC_HEIGHT-MM_time(ti))))'])]);
            end
        end
        Dist_m{mi}(:,gi) = Dist_mi';
    end
    toc
    disp(['Done with glider #' num2str(gi) ' out of ' num2str(length(GG_list))])
end

%% Calculate a boolean vector for each mooring and the glider
% These vectors are true when either instrument is within BOTH 2 km AND 1
% hour of the other type.

KM_cutoff = 2;
HR_cutoff = 1;

Mooring_2km_1hr = cell(length(Dist_m),1);
for mi = 1:length(MM_list)
    Dist_m_mi = Dist_m{mi};
    Dist_m_mi(~isfinite(Dist_m_mi)) = Inf;
    Mooring_2km_1hr{mi} = [Dist_m_mi < KM_cutoff] & [TOffset_m{mi} < HR_cutoff/24];
end

Glider_2km_1hr = cell(length(GG_list),1);
for gi = 1:length(GG_list)
    Dist_g_nan2inf = Dist_g{gi}; Dist_g_nan2inf(~isfinite(Dist_g_nan2inf)) = Inf;
    Glider_2km_1hr{gi} = [Dist_g{gi} < KM_cutoff] & [TOffset_g{gi}(:,mi) < HR_cutoff/24];
end

error('Forced error to stop the full running of the script.')

%% VISUAL CHECK: Compare steric height from both instruments in time

% {'ru32_'}
% {'ru38_'}
gi = 2;
GL = eval(GG_list{gi});

close all

COLORS = rand(length(MM_list),3);
% COLORS = colorcube(11);

figure
subplot(211) % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    plot(eval([MM '.TIME_STERIC_HEIGHT']) + T0,...
         eval([MM '.STERIC_HEIGHT_ANOMALY']),'.-','Color',COLORS(mi,:)); hold on
end
plot(GL.TIME_STERIC_HEIGHT + T0,...
     GL.STERIC_HEIGHT_ANOMALY,'k.--'); hold on
datetick
title(GG_list{gi})
subplot(212) % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    plot(eval([MM '.TIME_STERIC_HEIGHT(Mooring_2km_1hr{mi}(:,gi))']) + T0,...
         eval([MM '.STERIC_HEIGHT_ANOMALY(Mooring_2km_1hr{mi}(:,gi))']),'.-','Color',COLORS(mi,:)); hold on
end
for mi = 1:length(MM_list)
    plot(   GL.TIME_STERIC_HEIGHT(Glider_2km_1hr{gi}(:,mi)) + T0,...
         GL.STERIC_HEIGHT_ANOMALY(Glider_2km_1hr{gi}(:,mi)),'o--','Color',COLORS(mi,:),'HandleVisibility','off'); hold on
end
legend(MM_list)
XLIM_212 = get(gca,'XLim');
datetick
title('o = glider, . = mooring')
% % % 
subplot(211)
set(gca,'XLim',XLIM_212);

%% Scatter the depths and see if there is anything awry

close all
clc
disp(' '); disp(' '); disp('% % % % % % % % % % % % % % % % %')
% Mooring data slightly denser, so scatter glider data with nearest-time
% mooring steric height:
for gi = 1:length(GG_list)
    ColocatedGliderDEPTH_atMtime = {};
    LEG_list = {}; II = 1;
    ColocatedMooringDEPTH_all = [];
    ColocatedGliderDEPTH_atMtime_all = [];

    GL = eval(GG_list{gi});
    figure(gi)
    for mi = 1:length(MM_list)
        MM = MM_list{mi};
        ColocatedMooringDEPTH_mi = eval([MM '.STERIC_HEIGHT_MAX_PROF_DEPTH(Mooring_2km_1hr{mi}(:,gi))']) - ...
                                   eval([MM '.STERIC_HEIGHT_MIN_PROF_DEPTH(Mooring_2km_1hr{mi}(:,gi))']);
        ColocatedMooringT_mi = T0 + eval([MM '.TIME_STERIC_HEIGHT(Mooring_2km_1hr{mi}(:,gi))']);
        ColocatedMooringDEPTH_all = [ColocatedMooringDEPTH_all; ColocatedMooringDEPTH_mi];
        %
        ColocatedGliderDEPTH_mi = GL.STERIC_HEIGHT_MAX_PROF_DEPTH(Glider_2km_1hr{gi}(:,mi)) - ...
                                  GL.STERIC_HEIGHT_MIN_PROF_DEPTH(Glider_2km_1hr{gi}(:,mi));
        ColocatedGliderT_mi = T0 + GL.TIME_STERIC_HEIGHT(Glider_2km_1hr{gi}(:,mi));
        %
        ColocatedGliderDEPTH_atMtime_mi = nan(size(ColocatedMooringT_mi));
        for ti = 1:length(ColocatedMooringT_mi)
            ColocatedGliderDEPTH_atMtime_mi(ti) = ColocatedGliderDEPTH_mi(dsearchn(ColocatedGliderT_mi,ColocatedMooringT_mi(ti)));
        end
        ColocatedGliderDEPTH_atMtime{mi} = ColocatedGliderDEPTH_atMtime_mi;
        ColocatedGliderDEPTH_atMtime_all = [ColocatedGliderDEPTH_atMtime_all; ColocatedGliderDEPTH_atMtime_mi];
        if isempty(ColocatedGliderDEPTH_atMtime{mi}); HandVis = 'off'; else; HandVis = 'on'; LEG_list{II} = MM_list{mi}; II = II + 1; end
        %
        SCATTER = scatter(ColocatedMooringDEPTH_mi,ColocatedGliderDEPTH_atMtime{mi},100,'.','HandleVisibility',HandVis); hold on
        %
    end
    xlabel('Mooring max-min depth (m)'); ylabel('Glider max-min depth (m)'); axis equal; grid on

    title(['Glider ' GG_list{gi}])
    legend(LEG_list)

end


%% OFFSET and compare steric height from both instruments when they are approximately colocated and contemporary

close all
clc
FigTextColor = [1 1 1];

CUTOFF = 475;
disp(['Limit to ' num2str(CUTOFF) ' meters of profiled column'])
disp(' '); disp(' '); disp('% % % % % % % % % % % % % % % % %')
COEF_all_all = {};
% Mooring data slightly denser, so scatter glider data with nearest-time
% mooring steric height:
for gi = 1:length(GG_list)
    ColocatedGliderSH_atMtime = {};
    ColocatedGliderIND = {};
    LEG_list = {}; II = 1;
    ColocatedMooringSH_all = [];
    ColocatedGliderSH_atMtime_all = [];

    GL = eval(GG_list{gi});

    % figure(gi)
    if strcmp(GG_list{gi}(1:2),'uw')
        figure('Color',[50,0,111]/255)%0.7 0.1 0.7
    elseif strcmp(GG_list{gi}(1:2),'ru')
        figure('Color',[204 0 51]/255)%0.7 0.2 0.2
    elseif strcmp(GG_list{gi}(1:2),'ng')
        figure('Color',[2 42 58]/255)%0.2 0.5 0.2
    else
        warning('Unexpected result, what glider is this?')
    end
    for mi = 1:length(MM_list)
        MM = MM_list{mi};
        % CUTTOF_ind = eval(['[' MM '.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.STERIC_HEIGHT_MIN_PROF_DEPTH] > CUTOFF']);
        CUTTOF_ind = eval(['[' MM '.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM '.STERIC_HEIGHT_MIN_PROF_DEPTH] > CUTOFF & ' ...
                           '[' MM '.STERIC_HEIGHT_MIN_PROF_DEPTH < 5]']);
        ColocatedMooringSH_mi =   eval([MM '.STERIC_HEIGHT_ANOMALY(Mooring_2km_1hr{mi}(:,gi) & CUTTOF_ind)']);
        ColocatedMooringT_mi = T0 +  eval([MM '.TIME_STERIC_HEIGHT(Mooring_2km_1hr{mi}(:,gi) & CUTTOF_ind)']);
        ColocatedMooringIND_mi = [1:length(eval([MM '.STERIC_HEIGHT_ANOMALY']))]';
            ColocatedMooringIND_mi = ColocatedMooringIND_mi(Mooring_2km_1hr{mi}(:,gi) & CUTTOF_ind);
        ColocatedMooringSH_all = [ColocatedMooringSH_all; ColocatedMooringSH_mi];
        ColocatedMooringSeparation = Dist_m{mi}(Mooring_2km_1hr{mi}(:,gi) & CUTTOF_ind,gi);
        ColocatedMooringTimeLag = TOffset_m{mi}(Mooring_2km_1hr{mi}(:,gi) & CUTTOF_ind,gi);
        %
        ColocatedGliderSH_mi = GL.STERIC_HEIGHT_ANOMALY( Glider_2km_1hr{gi}(:,mi));
        ColocatedGliderT_mi = T0 + GL.TIME_STERIC_HEIGHT(Glider_2km_1hr{gi}(:,mi));
        ColocatedGliderIND_mi = [1:length(GL.STERIC_HEIGHT_ANOMALY)]';
            ColocatedGliderIND_mi = ColocatedGliderIND_mi(Glider_2km_1hr{gi}(:,mi));
        %
        ColocatedGliderSH_atMtime_mi  = nan(size(ColocatedMooringT_mi));
        ColocatedGliderIND_atMtime_mi = nan(size(ColocatedMooringT_mi));
        for ti = 1:length(ColocatedMooringT_mi)
            ColocatedGliderSH_atMtime_mi(ti)  = ColocatedGliderSH_mi(dsearchn(ColocatedGliderT_mi,ColocatedMooringT_mi(ti)));
            ColocatedGliderIND_atMtime_mi(ti) = ColocatedGliderIND_mi(dsearchn(ColocatedGliderT_mi,ColocatedMooringT_mi(ti)));
        end
        ColocatedGliderSH_atMtime{mi} = ColocatedGliderSH_atMtime_mi;
        ColocatedGliderIND{mi} = ColocatedGliderIND_atMtime_mi;
        ColocatedGliderSH_atMtime_all = [ColocatedGliderSH_atMtime_all; ColocatedGliderSH_atMtime_mi];
        if isempty(ColocatedGliderSH_atMtime{mi}); HandVis = 'off'; else; HandVis = 'on'; LEG_list{II} = MM_list{mi}; II = II + 1; end
        %
        SCATTER = scatter(ColocatedMooringSH_mi,ColocatedGliderSH_atMtime{mi},100,'.','HandleVisibility',HandVis); hold on
        % SCATTER = scatter(ColocatedMooringSH_mi,ColocatedGliderSH_atMtime{mi},100,ColocatedMooringSeparation,'.','HandleVisibility',HandVis); hold on
        %           CB = colorbar; CB.Label.String = 'Separation (km)';
        % SCATTER = scatter(ColocatedMooringSH_mi,ColocatedGliderSH_atMtime{mi},100,24*ColocatedMooringTimeLag,'.','HandleVisibility',HandVis); hold on
        %           CB = colorbar; CB.Label.String = 'Lag (hr)';
        SCATTER.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Mooring',repmat({MM},length(ColocatedMooringIND_mi),1));
        SCATTER.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Mooring Profile Number',ColocatedMooringIND_mi);
        SCATTER.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Glider',repmat({GG_list{gi}},length(ColocatedMooringIND_mi),1));
        SCATTER.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Glider Profile Number',ColocatedGliderIND{mi});
        SCATTER.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Profile Separation (km)',ColocatedMooringSeparation);
        SCATTER.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Profile time lag (hr)',24*ColocatedMooringTimeLag);
        %
        COEF =         [ones(size(ColocatedMooringSH_mi)), ColocatedMooringSH_mi]\ColocatedGliderSH_atMtime_mi;
        GliderFit_mi = [ones(size(ColocatedMooringSH_mi)), ColocatedMooringSH_mi]*COEF;
        
        % plot(ColocatedMooringSH_mi, GliderFit_mi, '-', 'Color',SCATTER.CData,'HandleVisibility','off');
        plot(ColocatedMooringSH_mi, GliderFit_mi, 'k-', 'HandleVisibility','off');
    end
    set(gca,'XColor',FigTextColor,'YColor',FigTextColor,'GridColor',[1 1 1]*0.15)
    xlabel('Mooring Steric Height (m)','Color',FigTextColor); ylabel('Glider Steric Height (m)','Color',FigTextColor); axis equal; grid on
    % Fit of all scattered points together, regardless of mooring
    SH_ISFINITE = [isfinite(ColocatedGliderSH_atMtime_all) & isfinite(ColocatedMooringSH_all)];
    ColocatedGliderSH_atMtime_all = ColocatedGliderSH_atMtime_all(SH_ISFINITE);
    ColocatedMooringSH_all = ColocatedMooringSH_all(SH_ISFINITE);
    [ColocatedMooringSH_all, SortOrder] = sort(ColocatedMooringSH_all);
    ColocatedGliderSH_atMtime_all = ColocatedGliderSH_atMtime_all(SortOrder);
    % Sorting is not necessary, it's only done here so that the plotted line
    % doesn't go over itself.
    COEF_all = [ones(size(ColocatedMooringSH_all)), ColocatedMooringSH_all]\ColocatedGliderSH_atMtime_all
    GliderFit_all = [ones(size(ColocatedMooringSH_all)), ColocatedMooringSH_all]*COEF_all;
    plot(ColocatedMooringSH_all, GliderFit_all, 'k--', 'LineWidth',2,'HandleVisibility','off')

    if COEF_all(2) >= 0
        title({['Glider ' GG_list{gi} ' , STD = ' num2str(100*std(ColocatedMooringSH_all - ColocatedGliderSH_atMtime_all,'omitnan')) ' cm'];...
               ['Glider offset = ' num2str(100*COEF_all(1)*sign(COEF_all(2))) ' cm']},'FontSize',24,'Color',FigTextColor)
    elseif COEF_all(2) < 0
        title({['Glider ' GG_list{gi} ' , STD = ' num2str(100*std(-ColocatedMooringSH_all - ColocatedGliderSH_atMtime_all,'omitnan')) ' cm'];...
               ['Glider offset = ' num2str(100*COEF_all(1)*sign(COEF_all(2))) ' cm']},'FontSize',24,'Color',FigTextColor)
    else
        title({['Glider ' GG_list{gi} ' , STD = ' num2str(100*std(ColocatedMooringSH_all - ColocatedGliderSH_atMtime_all,'omitnan')) ' cm'];...
               ['Glider offset = ' num2str(100*COEF_all(1)*sign(COEF_all(2))) ' cm']},'FontSize',24,'Color',FigTextColor)
    end
    datacursormode(gcf,'on')

    legend(LEG_list)

    disp(' ')
    disp(['Add ' num2str(-100*COEF_all(1)) ' cm to ' GG_list{gi} ' glider steric height to match mooring steric height'])
    disp('% % % % % % % % % % % % % % % % %')
    disp(' '); disp(' ');

    COEF_all_all{gi} = COEF_all;
end
% legend(LEG_list)


% disp(' ')
% disp(['Add ' num2str(-100*COEF_all(1)) ' cm to glider steric height to match mooring steric height'])
% disp(' ')

%% Examine individual profiles to examine outliers (use tooltips on figures above)

Cursors = who('cursor_info*');
CC = eval(Cursors{end}); % CC = "Current Cursor"
% CC.Target.DataTipTemplate.DataTipRows.Label:
%     'X'
%     'Y'
%     'Mooring'
%     'Mooring Profile Number'
%     'Glider'
%     'Glider Profile Number'

T0 = datenum('1950-01-01 00:00:00');

figure(13)
% TIME, DEPTH, PROF_NUM
IND = [CC.Target.XData == CC.Position(1)] & [CC.Target.YData == CC.Position(2)];
MM = CC.Target.DataTipTemplate.DataTipRows(3).Value{IND}; % all elements will be the same here, but use IND anyway
MM_pnum = CC.Target.DataTipTemplate.DataTipRows(4).Value(IND);
GG = CC.Target.DataTipTemplate.DataTipRows(5).Value{IND}; % all elements will be the same here, but use IND anyway
GG_pnum = CC.Target.DataTipTemplate.DataTipRows(6).Value(IND);
MG_dist = CC.Target.DataTipTemplate.DataTipRows(7).Value(IND);
MG_tlag = CC.Target.DataTipTemplate.DataTipRows(8).Value(IND);
% 
plot(T0 + eval([MM  '.TIME(' MM '.PROF_NUM == MM_pnum & ~' MM '.QC_FLAG)']),...
          eval([MM '.DEPTH(' MM '.PROF_NUM == MM_pnum & ~' MM '.QC_FLAG)']), '.-'); hold on
plot(T0 + eval([GG  '.TIME(' GG '.PROF_NUM == GG_pnum & ~' GG '.QC_FLAG)']),...
          eval([GG '.DEPTH(' GG '.PROF_NUM == GG_pnum & ~' GG '.QC_FLAG)']), '.-');
datetick('x','yyyy-mm-dd HH:MM')
legend(MM,GG)
xlabel('Time')
ylabel('Depth (m)')
title({['Mooring (' MM ' ) profile = ' num2str(MM_pnum) ',     Glider (' GG ' ) profile = ' num2str(GG_pnum)] ; ...
       ['Separation = ' num2str(MG_dist) ' km,     Lag = ' num2str(MG_tlag) ' hr']})


%% VISUAL CHECK: Compare steric height from both instruments in time WITH OFFSET

% {'ng252'}
% {'ng351'}
% {'ng644'}
% {'ng738'}
% {'ng782'}
% {'ru32_'}
% {'ru38_'}
% {'uw180'}
% {'uw219'}
% {'uw220'}
% {'uw247'}
% {'uw248'}
gi = 2;
GL = eval(GG_list{gi});

close all

COLORS = rand(length(MM_list),3);
% COLORS = colorcube(11);

figure
subplot(211) % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    plot(eval([MM '.TIME_STERIC_HEIGHT']) + T0,...
         eval([MM '.STERIC_HEIGHT_ANOMALY']),'.-','Color',COLORS(mi,:)); hold on
end
plot(GL.TIME_STERIC_HEIGHT + T0,...
     GL.STERIC_HEIGHT_ANOMALY,'k.-'); hold on
datetick
title(GG_list{gi})
subplot(212) % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
for mi = 1:length(MM_list)
    MM = MM_list{mi};
    plot(eval([MM '.TIME_STERIC_HEIGHT(Mooring_2km_1hr{mi}(:,gi))']) + T0,...
         eval([MM '.STERIC_HEIGHT_ANOMALY(Mooring_2km_1hr{mi}(:,gi))']),'.-','Color',COLORS(mi,:)); hold on
end
for mi = 1:length(MM_list)
    plot(   GL.TIME_STERIC_HEIGHT(Glider_2km_1hr{gi}(:,mi)) + T0,...
         GL.STERIC_HEIGHT_ANOMALY(Glider_2km_1hr{gi}(:,mi)) - COEF_all_all{gi}(1),'o--','Color',COLORS(mi,:),'HandleVisibility','off'); hold on
end
% legend(MM_list)
XLIM_212 = get(gca,'XLim');
datetick
title('o = glider, . = mooring')
% % % 
subplot(211)
set(gca,'XLim',XLIM_212);


