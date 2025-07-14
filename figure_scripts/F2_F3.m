%% Figure for Manuscript
% Look at the relatively small California SWOT Cal/Val diamond for SWOT data:

% % Core analysis, L2
SWOT = load('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/SWOTdata_CCS_33_385.mat',...
            'time','ssha_karin_2','ssha_karin_2_qual','name','lat','lon', ...
            'internal_tide_hret','height_cor_xover','height_cor_xover_qual','cross_track_distance');
LEVEL = 2;

% % Alternative, L3
% SWOT = load('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/SWOTdata_L3_CCS_33_385.mat',...
%             'time','ssha_karin_2','ssha_karin_2_qual','name','lat','lon', ...
%             'internal_tide_hret','height_cor_xover','height_cor_xover_qual','cross_track_distance');
% LEVEL = 3;

% L3: ssha_filtered: 'Height of the sea surface anomaly with all
                    % corrections applied and with calibration, data
                    % selection and noise reduction (using Unet model)
                    % applied; see the product user manual for details.'

% DIR_P = dir('/Users/kachelein/Documents/JPL/work/07.2023_QC_moorings/PROFILERS/SWOTPOSTLAUNCH_L2_JPLQC*.nc');
DIR_P = dir('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/PROF/SWOTPOSTLAUNCH_L2_JPLQC_MOORING*.nc');
for ii = 1:length(DIR_P)
    tic
    filename = DIR_P(ii).name;
    Mname = filename(33:34);
    MOORINGS{ii} =  Mname;
    eval([Mname '.LAT = ncread([DIR_P(ii).folder ''/'' filename],''LATITUDE'');']);
    eval([Mname '.LON = ncread([DIR_P(ii).folder ''/'' filename],''LONGITUDE'');']);
    toc
end

wsm


ALONGTRACK_RANGE = [1:807];
TRACK = 50;
SLOPE = -5/3; SLOPE = -2;
SEGS = 3;

IND = 1;
ALONGTRACK_DISTANCE = abs([SWOT.lon{IND}(TRACK,:) + 1i*SWOT.lat{IND}(TRACK,:)] - [SWOT.lon{IND}(TRACK,1) + 1i*SWOT.lat{IND}(TRACK,1)]);

dx = 111*median(diff(ALONGTRACK_DISTANCE));

%% Calculate mean ascending and descending tracks in order to remove the geoid per JW's advice

warning(['Also plot the spectra resulting from the removed and unremoved ' ...
         'alongside each other for the supplementary data'])

INDS = 1:length(SWOT.ssha_karin_2);
JJ1 = 1; JJ2 = 1;
MeanAscendingTrack = [];
MeanDescendingTrack = [];
for jj = 1:length(INDS)
    IND = INDS(jj)
    DATA = [SWOT.ssha_karin_2{IND} + 0*SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
           [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];
    ASCENDING_DESCENDING = sign(diff((SWOT.lat{jj}(50,1:2))));
    if ASCENDING_DESCENDING == 1 % ascending
        MeanAscendingTrack(:,:,JJ1) = DATA;
        JJ1 = JJ1 + 1;
    elseif ASCENDING_DESCENDING == -1 % descending
        MeanDescendingTrack(:,:,JJ2) = DATA;
        JJ2 = JJ2 + 1;
    else
        error(' ')
    end
end
MeanAscendingTrack = mean(MeanAscendingTrack,3,'omitnan');
MeanDescendingTrack = mean(MeanDescendingTrack,3,'omitnan');

%% (FIG 2) PREP: Plot all tracks and spectra of all tracks over multiple passes, ascending and descending (non-tiled, segmented)

close all

warning('off','all')

% Set the set of wavenumbers to probe in advance:
SEGS = 1;
WINDOW_FUNCTION = 'hanning';
IND = 180;
DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
        [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];
[~,Freq_predefined,~] = nanspectrum1(DATA(20,:),dx,'km',SEGS,'-k',false,0,WINDOW_FUNCTION);
% % % % % % % % % % % % % % % % % % % % % % % 

TRACK_LENGTH = nan(length(SWOT.ssha_karin_2),1);
for ii=1:length(TRACK_LENGTH)
    TRACK_LENGTH(ii) = size(SWOT.ssha_karin_2{ii},2);
end

clear SPEC_MEDIAN FREQ
% ASCENDING_DESCENDING = [1142]; % 1142 = ascending, 807 = descending
II = 1;

INDS = 1:length(SWOT.ssha_karin_2);
% INDS = INDS(TRACK_LENGTH == A_D); % ascending (1142) or descending (807)
SPEC_MEAN_track = [];
WhichMethod = []; % to keep track to whether or not the pass was discarded

for jj = 1:length(INDS)
    IND = INDS(jj);
    DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
        [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];

    ALONGTRACK_RANGE = 1:size(DATA,2);
    ASCENDING_DESCENDING = sign(diff((SWOT.lat{jj}(50,1:2))));
    if ASCENDING_DESCENDING == 1 % ascending
        % % Shorten only up to coast
        % SEGS = 2;
        % ALONGTRACK_RANGE = ALONGTRACK_RANGE(1:570);
        DATA = DATA - MeanAscendingTrack;
        % figure(2);m_pcolor_centered(SWOT.lon{jj}(:,ALONGTRACK_RANGE),SWOT.lat{jj}(:,ALONGTRACK_RANGE),DATA(:,ALONGTRACK_RANGE));figure(1)%$
    elseif ASCENDING_DESCENDING == -1 % descending
        % SEGS = 4;
        % ALONGTRACK_RANGE = ALONGTRACK_RANGE(1:1140);
        DATA = DATA - MeanDescendingTrack;
        % figure(2);m_pcolor_centered(SWOT.lon{jj}(:,ALONGTRACK_RANGE),SWOT.lat{jj}(:,ALONGTRACK_RANGE),DATA(:,ALONGTRACK_RANGE));figure(1)%$
    else
        error(' ')
    end
    TRACK = 50;

    ALONGTRACK_DISTANCE = abs([SWOT.lon{IND}(TRACK,:) + 1i*SWOT.lat{IND}(TRACK,:)] - [SWOT.lon{IND}(TRACK,1) + 1i*SWOT.lat{IND}(TRACK,1)]);

    TS = []; SPEC = []; SPEC_NFRemoved = []; DATA_COVERAGE = [];
    % for ii = 3:67
    %     TRACK = ii;
    %     [Spec,Freq,~] = nanspectrum1(DATA(TRACK,ALONGTRACK_RANGE),dx,'km',SEGS,'-',false,0,'hanning');
    %     SPEC(:,ii) = Spec;
    % end
    for ii = 1:size(DATA,1)% 3:67
        TRACK = ii;
        [Spec,Freq,~] = nanspectrum1(DATA(TRACK,ALONGTRACK_RANGE),dx,'km',SEGS,'-k',false,0,WINDOW_FUNCTION);

        if sum(isfinite(DATA(TRACK,ALONGTRACK_RANGE))) > 100
            [Spec,Freq,~] = nunanspectrum(DATA(TRACK,ALONGTRACK_RANGE), ALONGTRACK_RANGE*dx, 'km', ...
                'Segments',SEGS,'Plot_option','-','Plot',false,'Window',WINDOW_FUNCTION,'Freq',Freq_predefined);
            WhichMethod = [WhichMethod;0];
        else % all nan, which is hard to deal with using nunanspectrum
            [Spec,Freq,~] = nanspectrum1(DATA(TRACK,ALONGTRACK_RANGE),dx,'km',SEGS,'-k',false,0,WINDOW_FUNCTION);
            Spec = nan*Spec; % i.e. discard (only used nanspectrum1 to get it to be the right shape with minimal computation and no errors)
            WhichMethod = [WhichMethod;1];
        end

        SPEC(:,ii) = Spec;

        NOISE_FLOOR = mean(Spec(dsearchn(Freq,1/8):end),'omitnan');
        SPEC_NFRemoved(:,ii) = Spec - NOISE_FLOOR;

        DATA_COVERAGE(ii,jj) = sum(isfinite(DATA(TRACK,ALONGTRACK_RANGE)))/length(DATA(TRACK,ALONGTRACK_RANGE));
    end
    SPEC(SPEC==0) = NaN; % This will only apply to the first two columns, which default to all zeros
    SPEC_ = mean(SPEC',"omitmissing")';
    SPEC__ = median(SPEC',"omitmissing")';
    SPEC_MEAN_track(:,jj) = SPEC_;
    SPEC_MEDIAN_track(:,jj) = SPEC__;

    % Just the edges. Along the cross-track index, [1:4 32:38 66:69] are
    % nans. Do 12 total.
    SPEC_EDGES_ = mean(SPEC(:,[5:7 29:31 39:41 63:65])',"omitmissing")';
    SPEC_EDGES__ = median(SPEC(:,[5:7 29:31 39:41 63:65])',"omitmissing")';
    SPEC_MEAN_track_edges(:,jj) = SPEC_EDGES_;
    SPEC_MEDIAN_track_edges(:,jj) = SPEC_EDGES__;

    % The middle. Also do 12 total.
    SPEC_CENTER_ = mean(SPEC(:,[16:21 49:54])',"omitmissing")';
    SPEC_CENTER__ = median(SPEC(:,[16:21 49:54])',"omitmissing")';
    SPEC_MEAN_track_center(:,jj) = SPEC_CENTER_;
    SPEC_MEDIAN_track_center(:,jj) = SPEC_CENTER__;

    % Noise floor removed (all tracks):
    SPEC_NFRemoved(SPEC_NFRemoved==0) = NaN; % This will only apply to the first two columns, which default to all zeros
    SPEC_NFRemoved_ = mean(SPEC_NFRemoved',"omitmissing")';
    SPEC_NFRemoved__ = median(SPEC_NFRemoved',"omitmissing")';
    SPEC_NFRemoved_MEAN_track(:,jj) = SPEC_NFRemoved_;
    SPEC_NFRemoved_MEDIAN_track(:,jj) = SPEC_NFRemoved__;
    % Edges
    SPEC_NFRemoved_EDGES_ = mean(SPEC_NFRemoved(:,[5:7 29:31 39:41 63:65])',"omitmissing")';
    SPEC_NFRemoved_EDGES__ = median(SPEC_NFRemoved(:,[5:7 29:31 39:41 63:65])',"omitmissing")';
    SPEC_NFRemoved_MEAN_track_edges(:,jj) = SPEC_NFRemoved_EDGES_;
    SPEC_NFRemoved_MEDIAN_track_edges(:,jj) = SPEC_NFRemoved_EDGES__;
    % Middle
    SPEC_NFRemoved_CENTER_ = mean(SPEC_NFRemoved(:,[16:21 49:54])',"omitmissing")';
    SPEC_NFRemoved_CENTER__ = median(SPEC_NFRemoved(:,[16:21 49:54])',"omitmissing")';
    SPEC_NFRemoved_MEAN_track_center(:,jj) = SPEC_NFRemoved_CENTER_;
    SPEC_NFRemoved_MEDIAN_track_center(:,jj) = SPEC_NFRemoved_CENTER__;

    % figure(2);m_pcolor_centered(SWOT.lon{jj},SWOT.lat{jj},DATA);pause(0.2);figure(1)%$
    disp(jj)
end

%% FIG 2 PLOT

close all

figure('Color',[1 1 1])

fill([1/100 1/100 1/20 1/20],10.^[-6 -4.5 -4.5 -6],   [1 1 1]*0.9,'EdgeColor','none'); hold on % mooring range

SPEC_MEAN_track          = SPEC_MEAN_track .* SPEC_MEAN_track./SPEC_MEAN_track; % turn 0 to nan
SPEC_MEAN_track_edges    = SPEC_MEAN_track_edges .* SPEC_MEAN_track_edges./SPEC_MEAN_track_edges; % turn 0 to nan
SPEC_MEAN_track_center   = SPEC_MEAN_track_center .* SPEC_MEAN_track_center./SPEC_MEAN_track_center; % turn 0 to nan
SPEC_MEAN                = mean(  SPEC_MEAN_track,2,'omitnan');
SPEC_MEAN_edges          = mean(  SPEC_MEAN_track_edges,2,'omitnan');
SPEC_MEAN_center         = mean(  SPEC_MEAN_track_center,2,'omitnan');

SPEC_MEDIAN_track        = SPEC_MEDIAN_track .* SPEC_MEDIAN_track./SPEC_MEDIAN_track; % turn 0 to nan
SPEC_MEDIAN_track_edges  = SPEC_MEDIAN_track_edges .* SPEC_MEDIAN_track_edges./SPEC_MEDIAN_track_edges; % turn 0 to nan
SPEC_MEDIAN_track_center = SPEC_MEDIAN_track_center .* SPEC_MEDIAN_track_center./SPEC_MEDIAN_track_center; % turn 0 to nan
SPEC_MEDIAN              = mean(  SPEC_MEDIAN_track,2,'omitnan');
SPEC_MEDIAN_edges        = mean(  SPEC_MEDIAN_track_edges,2,'omitnan');
SPEC_MEDIAN_center       = mean(  SPEC_MEDIAN_track_center,2,'omitnan');

SPEC_NFRemoved_MEAN_track          = SPEC_NFRemoved_MEAN_track .* SPEC_NFRemoved_MEAN_track./SPEC_NFRemoved_MEAN_track; % turn 0 to nan
SPEC_NFRemoved_MEAN_track_edges    = SPEC_NFRemoved_MEAN_track_edges .* SPEC_NFRemoved_MEAN_track_edges./SPEC_NFRemoved_MEAN_track_edges; % turn 0 to nan
SPEC_NFRemoved_MEAN_track_center   = SPEC_NFRemoved_MEAN_track_center .* SPEC_NFRemoved_MEAN_track_center./SPEC_NFRemoved_MEAN_track_center; % turn 0 to nan
SPEC_NFRemoved_MEAN                = mean(  SPEC_NFRemoved_MEAN_track,2,'omitnan');
SPEC_NFRemoved_MEAN_edges          = mean(  SPEC_NFRemoved_MEAN_track_edges,2,'omitnan');
SPEC_NFRemoved_MEAN_center         = mean(  SPEC_NFRemoved_MEAN_track_center,2,'omitnan');

SPEC_NFRemoved_MEDIAN_track        = SPEC_NFRemoved_MEDIAN_track .* SPEC_NFRemoved_MEDIAN_track./SPEC_NFRemoved_MEDIAN_track; % turn 0 to nan
SPEC_NFRemoved_MEDIAN_track_edges  = SPEC_NFRemoved_MEDIAN_track_edges .* SPEC_NFRemoved_MEDIAN_track_edges./SPEC_NFRemoved_MEDIAN_track_edges; % turn 0 to nan
SPEC_NFRemoved_MEDIAN_track_center = SPEC_NFRemoved_MEDIAN_track_center .* SPEC_NFRemoved_MEDIAN_track_center./SPEC_NFRemoved_MEDIAN_track_center; % turn 0 to nan
SPEC_NFRemoved_MEDIAN              = mean(  SPEC_NFRemoved_MEDIAN_track,2,'omitnan');
SPEC_NFRemoved_MEDIAN_edges        = mean(  SPEC_NFRemoved_MEDIAN_track_edges,2,'omitnan');
SPEC_NFRemoved_MEDIAN_center       = mean(  SPEC_NFRemoved_MEDIAN_track_center,2,'omitnan');

FREQ{II} = Freq;

% % Presentation Figure

% % % Range in which the moorings observe (subject to modification)
% fill([1/100 1/100 1/10 1/10],10.^[-6 1 1 -6],   [1 1 1]*0.9,'EdgeColor','none','HandleVisibility','off'); hold on
set(gca,'YScale','log','XScale','log')

% % % Average spectrum for each track
loglog(Freq,SPEC_MEAN_track,'-','Color',[1 1 1]*0.7,'HandleVisibility','off'); hold on


% % % Median spectrum, slope line, text, and the rest
SLOPE = -2;
xlim(Freq([1 end-1]))
ylim(10.^[-6 1])
set(gca,'FontSize',14)
xlabel('Wavenumber (cpkm)') % xlabel('Wavenumber (km^{-1})')
ylabel('Spectral power density (m^2 cpkm^{-1})')



loglog(nan,nan,'Color',[1 1 1]*0.5)

PLOT_SLOPE_1 = -2;
PLOT_SLOPE_2 = -11/3;

loglog( Freq(dsearchn(Freq,0.01):end), ...
        10^-6.7 * Freq(dsearchn(Freq,0.01):end).^PLOT_SLOPE_1, 'k--', 'LineWidth',3,'HandleVisibility','off')
text(Freq(1)*10,10^-3.85,['\propto k^{' num2str(PLOT_SLOPE_1) '}'],'FontSize',16,'HandleVisibility','off')

loglog( Freq(1:dsearchn(Freq,0.02)), ...
        10^-8.3 * Freq(1:dsearchn(Freq,0.02)).^[PLOT_SLOPE_2], 'k--', 'LineWidth',3,'HandleVisibility','off')
% text(Freq(1)*12,10^-1.5,['\propto k^{' num2str(PLOT_SLOPE_2) '}'],'FontSize',16,'HandleVisibility','off')
text(Freq(1)*7,10^-0.85,['\propto k^{-11/3}'],'FontSize',16,'HandleVisibility','off')

% MEAN
loglog(FREQ{1},SPEC_MEAN,'LineWidth',3,'Color',[24;116;205]'/255)
loglog(FREQ{1},SPEC_NFRemoved_MEAN,':','LineWidth',3,'Color',[24;116;205]'/255,'HandleVisibility','off')

loglog(FREQ{1},SPEC_MEAN_edges,'LineWidth',3,'Color',[217.6000   83.2000   25.0880]/256)
loglog(FREQ{1},SPEC_NFRemoved_MEAN_edges,':','LineWidth',3,'Color',[217.6000   83.2000   25.0880]/256,'HandleVisibility','off')

loglog(FREQ{1},SPEC_MEAN_center,'LineWidth',3,'Color',[77.0560  190.7200  238.8480]/256)
loglog(FREQ{1},SPEC_NFRemoved_MEAN_center,':','LineWidth',3,'Color',[77.0560  190.7200  238.8480]/256,'HandleVisibility','off')

loglog(FREQ{1}(1),NaN,':','LineWidth',3,'Color',[128 128 128]/256,'HandleVisibility','on') % for the legend

% % MEDIAN
% loglog(FREQ{1},SPEC_MEDIAN,'LineWidth',3,'Color',[24;116;205]'/255)
% loglog(FREQ{1},SPEC_MEDIAN_edges,'LineWidth',3,'Color',[217.6000   83.2000   25.0880]/256)
% loglog(FREQ{1},SPEC_MEDIAN_center,'LineWidth',3,'Color',[77.0560  190.7200  238.8480]/256)

% legend('Glider observation range','Mooring observation range','Track mean spectra','Median spectrum (ascending)','Median spectrum (descending)')
legend(['Mooring wavenumber' char(10) 'resolution'], ...
    'Track-mean spectra',...
    'Mean spectrum', ...
    ['Mean spectrum' char(10) '(edges of tracks)'], ...
    ['Mean spectrum' char(10) '(centers of tracks)'], ...
    ['{\it k} > 1/8 cpkm' char(10) 'noise floor removed'])
% title({'Along-track Spectra: Cal/Val Orbit';'US West Coast Crossover'},'FontWeight','bold','FontSize',16)
set(gcf,'Position',[-1776         228         813         730])
set(gca,'TickLength',[1 1]/50)
set(gca,'LineWidth',2)
if LEVEL == 2
    ylim(10.^[-5.5 0])
elseif LEVEL == 3
    ylim(10.^[-7.5 0])
else
end

fontsize(22,"point")

warning('on','all')

grid on

% % 70 km line label:
% text(1/67,10^-4.8,'(70 km)^{-1}','FontSize',20,'HandleVisibility','off')

%% Save

% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/F2.pdf',...
% 'BackgroundColor','none','ContentType','vector')

% % % Save a second figure for revised manuscript:
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/REV1/F2.pdf',...
% 'BackgroundColor','none','ContentType','vector')

% % % Save level 3 version for supplemental material:
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/REV1/F2_L3.pdf',...
% 'BackgroundColor','none','ContentType','vector')

%% Structure functions from full tracks, every cell

close all

% error(['This calculation takes a while'])

SF_ARRAY_swot = nan(size(SWOT.ssha_karin_2{1},1), size(SWOT.ssha_karin_2{1},2) - 1, length(SWOT.ssha_karin_2));
for jj = 1:length(SWOT.ssha_karin_2)
    DATA = [  SWOT.ssha_karin_2{jj}      + SWOT.internal_tide_hret{jj} + SWOT.height_cor_xover{jj}].*...
           [~(SWOT.ssha_karin_2_qual{jj} + SWOT.height_cor_xover_qual{jj})./...
            ~(SWOT.ssha_karin_2_qual{jj} + SWOT.height_cor_xover_qual{jj})];

    ASCENDING_DESCENDING = sign(diff((SWOT.lat{jj}(50,1:2))));
    if ASCENDING_DESCENDING == 1 % ascending (long)
        DATA = DATA - MeanAscendingTrack;
    elseif ASCENDING_DESCENDING == -1 % descending (short)
        DATA = DATA - MeanDescendingTrack;
    else
        error(' ')
    end

    % Structure function analysis for this track

    SF_MAT_swot_i = nan(size(DATA) - [0 1]);
    for track_i = 1:size(DATA,1)
        TRACK_i = DATA(track_i,:);
        LAT_i = SWOT.lat{jj}(track_i,:);
        LON_i = SWOT.lon{jj}(track_i,:);
        dx = median(m_lldist(LON_i,LAT_i),'omitnan');
        DIST_VEC = nan(length(LAT_i)-1,1);
        SF_VEC =   nan(length(LAT_i)-1,1);
        for di = 1:length(DIST_VEC)
            DIST_VEC(di) = di*dx;
            SF_VEC(di) = mean([TRACK_i(1:[end-di]) - TRACK_i([1+di]:end)].^2,'omitnan');
        end
        SF_MAT_swot_i(track_i,:) = SF_VEC;
    end
    SF_ARRAY_swot(:,:,jj) = SF_MAT_swot_i;
    disp(jj)
end
% %
SF_MAT_swot = nan(size(SF_ARRAY_swot,2),size(SF_ARRAY_swot,1)*size(SF_ARRAY_swot,3));
for kk = 1:size(SF_ARRAY_swot,3)
    SF_MAT_swot(:,[1:size(SF_ARRAY_swot,1)] + size(SF_ARRAY_swot,1)*[kk-1]) = SF_ARRAY_swot(:,:,kk)';
end


%% SWOT STRUCTURE FUNCTIONS

load(['/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/' ...
      'SF_plot_fixed_bottom_windowed_SWOTt_data.mat'],...
      'DD','P1','P2','P3','P5','P6','P7','S1','S2','S3','S4','ru32_600','ru38_600');

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

T0 = datenum('01-Jan-1950');

% Set vector of nominal locations of moorings:
MM_list = {'P1','P2','P3','P5','P6','P7','S1','S2','S3','S4'};
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

        if LEVEL == 2
            SWOT_M_rollerror_ti = SWOT.height_cor_xover{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
            SWOT_M_rollerror = [SWOT_M_rollerror ; SWOT_M_rollerror_ti];

            SWOT_M_rollerror_qc_ti = SWOT.height_cor_xover_qual{ti}( Mooring_Loc_Mat==min(Mooring_Loc_Mat(:)) );
            SWOT_M_rollerror_qc = [SWOT_M_rollerror_qc ; SWOT_M_rollerror_qc_ti];
        elseif LEVEL == 3
            SWOT_M_rollerror = 0;
            SWOT_M_rollerror_qc = 0;
        else
        end


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


%%
%%
%%
%%
%% Look at slopes from along-track data and Welch's method

close all

II = 1;

INDS = 1:length(SWOT.ssha_karin_2);
% INDS = INDS(TRACK_LENGTH == A_D); % ascending (1142) or descending (807)
WINDOW_FUNCTION = 'rectwin';
SEGS = 1;
SPEC_MEAN_track = [];
SPEC = [];
DATA_COVERAGE = [];

% Set the set of wavenumbers to probe in advance:
IND = 180;
DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
        [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];
[~,Freq_predefined,~] = nanspectrum1(DATA(20,:),dx,'km',SEGS,'-k',false,0,WINDOW_FUNCTION);
% % % % % % % % % % % % % % % % % % % % % % % 

for jj = 1:length(INDS)
    IND = INDS(jj);
    DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
        [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];

    ALONGTRACK_RANGE = 1:size(DATA,2);
    ASCENDING_DESCENDING = sign(diff((SWOT.lat{jj}(50,1:2))));
    if ASCENDING_DESCENDING == 1 % ascending
        % % Shorten only up to coast
        % SEGS = 2;
        % ALONGTRACK_RANGE = ALONGTRACK_RANGE(1:570);
        DATA = DATA - MeanAscendingTrack;
        % figure(2);m_pcolor_centered(SWOT.lon{jj}(:,ALONGTRACK_RANGE),SWOT.lat{jj}(:,ALONGTRACK_RANGE),DATA(:,ALONGTRACK_RANGE));figure(1)%$
    elseif ASCENDING_DESCENDING == -1 % descending
        % SEGS = 4;
        % ALONGTRACK_RANGE = ALONGTRACK_RANGE(1:1140);
        DATA = DATA - MeanDescendingTrack;
        % figure(2);m_pcolor_centered(SWOT.lon{jj}(:,ALONGTRACK_RANGE),SWOT.lat{jj}(:,ALONGTRACK_RANGE),DATA(:,ALONGTRACK_RANGE));figure(1)%$
    else
        error(' ')
    end
    TRACK = 50;

    ALONGTRACK_DISTANCE = abs([SWOT.lon{IND}(TRACK,:) + 1i*SWOT.lat{IND}(TRACK,:)] - [SWOT.lon{IND}(TRACK,1) + 1i*SWOT.lat{IND}(TRACK,1)]);
    FFT_vs_NUFFT = []; % 0 for FFT spectrum, 1 for NUFFT spectrum

    for ii = 1:69
        TRACK = ii;
        [Spec,Freq,~] = nanspectrum1(DATA(TRACK,ALONGTRACK_RANGE),dx,'km',SEGS,'-k',false,0,WINDOW_FUNCTION);
        
        if sum(isfinite(DATA(TRACK,ALONGTRACK_RANGE))) > 100
            FFT_vs_NUFFT = [FFT_vs_NUFFT; 1];
            [Spec,Freq,~] = nunanspectrum(DATA(TRACK,ALONGTRACK_RANGE), ALONGTRACK_RANGE*dx, 'km', ...
                'Segments',SEGS,'Plot_option','-','Plot',false,'Window',WINDOW_FUNCTION,'Freq',Freq_predefined);
        else % all nan, which is hard to deal with using nunanspectrum
            FFT_vs_NUFFT = [FFT_vs_NUFFT; 0];
            [Spec,Freq,~] = nanspectrum1(DATA(TRACK,ALONGTRACK_RANGE),dx,'km',SEGS,'-k',false,0,WINDOW_FUNCTION);
            Spec = nan*Spec;
        end

        SPEC(:,ii,jj) = Spec;
        DATA_COVERAGE(ii,jj) = sum(isfinite(DATA(TRACK,ALONGTRACK_RANGE)))/length(DATA(TRACK,ALONGTRACK_RANGE));
    end
    % figure(2);m_pcolor_centered(SWOT.lon{jj},SWOT.lat{jj},DATA);pause(0.2);figure(1)%$
    disp([num2str(jj) ' - ' num2str(100*sum(FFT_vs_NUFFT)/length(FFT_vs_NUFFT)) '% tracks using NUFFT'])
end

SPEC_ = mean(SPEC,3,"omitmissing");

%% Figure of the mean spectrum for each cross-track distance:

close all

ApproxCrossTrackDist = nan(size(SWOT.cross_track_distance{1},1),length(SWOT.cross_track_distance));
for ii = 1:length(SWOT.cross_track_distance)
    ApproxCrossTrackDist(:,ii) = mean(SWOT.cross_track_distance{ii},2);
end
ApproxCrossTrackDist = mean(ApproxCrossTrackDist,2,'omitnan');

% Slope of each spectrum past some frequency (Freq is in cpkm)
KM_LONGEST = 100;
ind_min = dsearchn(Freq,1/[KM_LONGEST]);
ind_max = dsearchn(Freq,1/[20]);
% ind_max = length(Freq); % maximum
HH = [ones(length(Freq(ind_min:ind_max)),1), log10(Freq(ind_min:ind_max))];
SPEC_COEF = [];


%% Figure with multiple length scale regimes

close all

TURBO = [flip(turbo(34)) ; 0 0 0 ; turbo(34)];
% figure % verify that we are cutting it off as expected
% for ii = 1:size(SPEC_,2)
%     loglog(Freq,SPEC_(:,ii),'-','Color',TURBO(ii,:),'LineWidth',1); hold on
% end
% KM_LONGEST  = [200 100 60 10];
% KM_SHORTEST = [60  20  4  4];
KM_LONGEST  = [200 100 60];
KM_SHORTEST = [60  20  8 ];
CROSSTRACK_SLOPES = nan(size(SPEC_,2),length(KM_LONGEST));
for jj = 1:length(KM_LONGEST)
    ind_min = dsearchn(Freq,1/KM_LONGEST(jj));
    ind_max = dsearchn(Freq,1/KM_SHORTEST(jj));
    % ind_max = length(Freq); % maximum
    HH = [ones(length(Freq(ind_min:ind_max)),1), log10(Freq(ind_min:ind_max))];
    SPEC_COEF = [];
    for ii = 1:size(SPEC_,2)
        SPEC_COEF(ii,:) = [HH\log10(SPEC_(ind_min:ind_max,ii))]';
        % loglog(Freq(ind_min:ind_max), [10.^SPEC_COEF(ii,1) * Freq(ind_min:ind_max).^SPEC_COEF(ii,2)], 'k');hold on
    end
    CROSSTRACK_SLOPES(:,jj) = SPEC_COEF(:,2);
end

DX = 5;

KM_LONGEST  = [200 100 60];
KM_SHORTEST = [60  20  8 ];
CROSSTRACK_SLOPES = nan(size(SPEC_,2),length(KM_LONGEST));
for jj = 1:length(KM_LONGEST)
    ind_min = dsearchn(Freq,1/KM_LONGEST(jj));
    ind_max = dsearchn(Freq,1/KM_SHORTEST(jj));
    % ind_max = length(Freq); % maximum
    HH = [ones(length(Freq(ind_min:ind_max)),1), log10(Freq(ind_min:ind_max))];
    SPEC_COEF = [];
    for ii = 1:size(SPEC_,2)
        SPEC_COEF(ii,:) = [HH\log10(SPEC_(ind_min:ind_max,ii))]';
        % loglog(Freq(ind_min:ind_max), [10.^SPEC_COEF(ii,1) * Freq(ind_min:ind_max).^SPEC_COEF(ii,2)], 'k');hold on
    end
    CROSSTRACK_SLOPES(:,jj) = SPEC_COEF(:,2);
end

%% Better approach: Slope(TimeMean(SPEC - N.F.))

ApproxCrossTrackDist = nan(size(SWOT.cross_track_distance{1},1),length(SWOT.cross_track_distance));
for ii = 1:length(SWOT.cross_track_distance)
    ApproxCrossTrackDist(:,ii) = mean(SWOT.cross_track_distance{ii},2);
end
ApproxCrossTrackDist = mean(ApproxCrossTrackDist,2,'omitnan');
if LEVEL==2
elseif LEVEL==3
    ApproxCrossTrackDist = ApproxCrossTrackDist*1000;
else
    error('LEVEL must be 2 or 3')
end

close all

figure % verify that we are cutting it off as expected
KM_LONGEST  = [200 100 60];
KM_SHORTEST = [60  20  8 ];
CROSSTRACK_SLOPES_NFR = nan(size(SPEC,2), length(KM_LONGEST));
% ^ NFR = "noise floor removed"
for nn = 1:length(KM_LONGEST)
    ind_min = dsearchn(Freq,1/KM_LONGEST(nn));
    ind_max = dsearchn(Freq,1/KM_SHORTEST(nn));
    % ind_max = length(Freq); % maximum
    HH = [ones(length(Freq(ind_min:ind_max)),1), log10(Freq(ind_min:ind_max))];
    SPEC_COEF = [];
    for jj = 1:size(SPEC,2)
        NOISE_FLOOR = mean(SPEC(dsearchn(Freq,1/8):end,jj),'omitnan');
        SPEC_jj = SPEC(ind_min:ind_max,jj) - NOISE_FLOOR;
        COEF_jj = [HH(SPEC_jj>0,:)\log10(SPEC_jj(SPEC_jj>0))]';
        SPEC_COEF(1,jj) = COEF_jj(1); SPEC_COEF(2,jj) = COEF_jj(2);
    end
    CROSSTRACK_SLOPES_NFR(:,nn) = SPEC_COEF(2,:);
end
CROSSTRACK_SLOPES_NFR(CROSSTRACK_SLOPES_NFR==0) = NaN;

figure('Color','w')
plot(ApproxCrossTrackDist, mean(squeeze(CROSSTRACK_SLOPES_NFR),2,'omitnan') , '.-', 'LineWidth',1,'MarkerSize',20); hold on
LEG_NFR = {};
for nn = 1:length(KM_LONGEST)
% for nn = [1 3]
    LEG_NFR{nn} = [num2str(KM_LONGEST(nn)) ' to ' num2str(KM_SHORTEST(nn)) ' km - noise floor removed'];
end
legend(LEG_NFR)
xlabel('Approx. distance from nadir')
ylabel(['Spectral slope of average spectrum'])

% With the previous results (before removing the noise floor):
figure('Color','w')
% LINE_COLOR = colororder;
LINE_COLOR = [.1 .1 .8;   .1 .8 .1;   .8 .1 .1];

% % Omit the edges manually
% plot(ApproxCrossTrackDist,-2*ApproxCrossTrackDist./ApproxCrossTrackDist,'.-k');hold on % For index counting
CrossTrackPlotIndices = 1:65;%[1:29 32:(length(ApproxCrossTrackDist)-4)];
% CrossTrackPlotIndices = [1:length(ApproxCrossTrackDist)];
for nn = 1:size(CROSSTRACK_SLOPES,2)
    plot([10^-3]*ApproxCrossTrackDist(CrossTrackPlotIndices), CROSSTRACK_SLOPES(CrossTrackPlotIndices,nn), ...
        '.-', 'Color', LINE_COLOR(nn,:), 'LineWidth',1,'MarkerSize',20); hold on
end
for nn = 1:size(CROSSTRACK_SLOPES_NFR,2)
    plot([10^-3]*ApproxCrossTrackDist(CrossTrackPlotIndices), CROSSTRACK_SLOPES_NFR(CrossTrackPlotIndices,nn), ...
        'o-', 'Color', LINE_COLOR(nn,:), 'LineWidth',1,'MarkerSize',10); hold on
end
grid on
grid minor

LEG = {};
LEG_NFR = {};
for nn = 1:length(KM_LONGEST)
    LEG{nn} = [num2str(KM_LONGEST(nn)) ' to ' num2str(KM_SHORTEST(nn)) ' km'];
    LEG_NFR{nn} = [num2str(KM_LONGEST(nn)) ' to ' num2str(KM_SHORTEST(nn)) ' km - noise floor removed'];
end
LEG_HANDLE = legend([LEG]);
set(gca,'XLim',[-70 70])
xlabel('Distance from nadir (km)')
ylabel(['Spectral slope of average spectrum'])

fontsize(24,"point")
LEG_HANDLE.FontSize = 16;
% LEG_HANDLE.Position = [0.6451 0.7633 0.1875 0.1378];
% LEG_HANDLE.Position = [0.5816    0.7726    0.1875    0.1378];
LEG_HANDLE.Position = [0.5983    0.7726    0.1875    0.1378];

%% Export the figure

% figure(3)
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/F3.pdf',...
% 'BackgroundColor','none','ContentType','vector')

% % % Revision 1 redone figure:
% figure(3)
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/REV1/F3.pdf',...
% 'BackgroundColor','none','ContentType','vector')


% % % Supplemental figure for L3 data:
% figure(3)
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/REV1/F3_L3.pdf',...
% 'BackgroundColor','none','ContentType','vector')

%% Evaluate the effects of varying the cutoff wavenumber and noise-floor as a function of cross-track-distance:

% ApproxCrossTrackDist = nan(size(SWOT.cross_track_distance{1},1),length(SWOT.cross_track_distance));
% for ii = 1:length(SWOT.cross_track_distance)
%     ApproxCrossTrackDist(:,ii) = mean(SWOT.cross_track_distance{ii},2);
% end
% ApproxCrossTrackDist = mean(ApproxCrossTrackDist,2,'omitnan');

% k_cutoff = [1/8 3/16];
k_cutoff = 1./[8 7 6 5];

NOISE_FLOOR_MAT = nan(length(k_cutoff),size(SPEC,2),size(SPEC,3));
% ^ NFR = "noise floor removed"
for ki = 1:length(k_cutoff)
    for jj = 1:size(SPEC,2)
        for kk = 1:size(SPEC,3)
            NOISE_FLOOR = mean(SPEC(dsearchn(Freq,k_cutoff(ki)):end,jj,kk),'omitnan');
            NOISE_FLOOR_MAT(ki,jj,kk) = NOISE_FLOOR;
        end
    end
end

k_lowcutoff = 1./[100 50 25 16];
LOW_K_MAT = nan(length(k_lowcutoff),size(SPEC,2),size(SPEC,3));
% ^ 
for ki = 1:length(k_lowcutoff)
    for jj = 1:size(SPEC,2)
        for kk = 1:size(SPEC,3)
            LOW_K_SPEC = Freq(1)*sum(SPEC(1:dsearchn(Freq,k_lowcutoff(ki)),jj,kk),'omitnan');
            LOW_K_MAT(ki,jj,kk) = LOW_K_SPEC;
        end
    end
end

%% Plot noise floor vs. xtd
close all
figure('Color','w')
plot([10^-3]*ApproxCrossTrackDist,mean(squeeze(NOISE_FLOOR_MAT(1,:,:)),2,'omitnan'),'.-','LineWidth',1,'MarkerSize',15); hold on
plot([10^-3]*ApproxCrossTrackDist,mean(squeeze(NOISE_FLOOR_MAT(2,:,:)),2,'omitnan'),'.-','LineWidth',1,'MarkerSize',15)
plot([10^-3]*ApproxCrossTrackDist,mean(squeeze(NOISE_FLOOR_MAT(3,:,:)),2,'omitnan'),'.-','LineWidth',1,'MarkerSize',15)
plot([10^-3]*ApproxCrossTrackDist,mean(squeeze(NOISE_FLOOR_MAT(4,:,:)),2,'omitnan'),'.-','LineWidth',1,'MarkerSize',15)
LEG = ...
legend(['k_{cutoff} = 1/' num2str(1/k_cutoff(1)) ' cpkm'],...
       ['k_{cutoff} = 1/' num2str(1/k_cutoff(2)) ' cpkm'],...
       ['k_{cutoff} = 1/' num2str(1/k_cutoff(3)) ' cpkm'],...
       ['k_{cutoff} = 1/' num2str(1/k_cutoff(4)) ' cpkm']);
LEG.Position = [0.6133    0.6124    0.1787    0.2545];

set(gca,'XLim',[-70 70])
xlabel('Distance from nadir (km)')
ylabel(['Mean of S(k>k_{cutoff})'])
set(gca,'FOntSize',16)
set(gcf,'Position',[-1439         420         873         391])

figure('Color','w')
plot([10^-3]*ApproxCrossTrackDist,mean(squeeze(LOW_K_MAT(1,:,:)),2,'omitnan'),'.-','LineWidth',1,'MarkerSize',15); hold on
plot([10^-3]*ApproxCrossTrackDist,mean(squeeze(LOW_K_MAT(2,:,:)),2,'omitnan'),'.-','LineWidth',1,'MarkerSize',15)
plot([10^-3]*ApproxCrossTrackDist,mean(squeeze(LOW_K_MAT(3,:,:)),2,'omitnan'),'.-','LineWidth',1,'MarkerSize',15)
plot([10^-3]*ApproxCrossTrackDist,mean(squeeze(LOW_K_MAT(4,:,:)),2,'omitnan'),'.-','LineWidth',1,'MarkerSize',15)
LEG = ...
legend(['k_{max} = 1/' num2str(1/k_lowcutoff(1)) ' cpkm'],...
       ['k_{max} = 1/' num2str(1/k_lowcutoff(2)) ' cpkm'],...
       ['k_{max} = 1/' num2str(1/k_lowcutoff(3)) ' cpkm'],...
       ['k_{max} = 1/' num2str(1/k_lowcutoff(4)) ' cpkm']);
LEG.Position = [0.6202    0.2160    0.1879    0.2545];

set(gca,'XLim',[-70 70])
xlabel('Distance from nadir (km)')
ylabel(['Sum of \Deltak\timesS(k<k_{max})'])
set(gca,'FOntSize',16)
set(gcf,'Position',[-1439         420         873         391])
