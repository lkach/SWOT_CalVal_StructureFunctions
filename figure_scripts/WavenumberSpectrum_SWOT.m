%% Figure for Manuscript and Data Analysis
% Compare to what MA did

SWOT = load('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/SWOTdata_CCS.mat',...
            'time','ssha_karin_2','ssha_karin_2_qual','name','lat','lon','internal_tide_hret','height_cor_xover','height_cor_xover_qual');
% SWOT = load('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/SWOTdata_CCS_33_385.mat',... More limited latitude range than above
%             'time','ssha_karin_2','ssha_karin_2_qual','name','lat','lon','internal_tide_hret','height_cor_xover','height_cor_xover_qual');

LL = load('/Users/kachelein/Documents/JPL/work/09.2023/SWOT_STM/MATLAB/LatLonLines.mat');

DIR_P = dir('/Users/kachelein/Documents/JPL/work/07.2023_QC_moorings/PROFILERS/SWOTPOSTLAUNCH_L2_JPLQC*.nc');
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

%% Visualize any pass

close all

figure

for IND = 1:length(SWOT.time)
    DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
           [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];

    pcolor(SWOT.lon{IND}, SWOT.lat{IND}, DATA); shading flat; hold on
    title(num2str(IND))
    set(gca,'XLim',[-132 -119],'YLim',[29.5 50.5])
    xlabel('LON'); ylabel('LAT')

    pause(.1)
end


%% Plot all tracks and spectra of all tracks over multiple passes (SEGMENTED)

close all
TRACK_LENGTH = nan(length(SWOT.ssha_karin_2),1);
for ii=1:length(TRACK_LENGTH)
    TRACK_LENGTH(ii) = size(SWOT.ssha_karin_2{ii},2);
end

INDS = 1:length(SWOT.ssha_karin_2);
% INDS = INDS(TRACK_LENGTH == 807); % descending
INDS = INDS(TRACK_LENGTH == 1142); % ascending
SEGS = 3;
% SPEC_MEAN = 0;
SPEC_MEAN_track = [];


TRACK = 35;
IND = 1;
ALONGTRACK_DISTANCE = abs([SWOT.lon{IND}(TRACK,:) + 1i*SWOT.lat{IND}(TRACK,:)] - [SWOT.lon{IND}(TRACK,1) + 1i*SWOT.lat{IND}(TRACK,1)]);
dx = 111*median(diff(ALONGTRACK_DISTANCE));

figure

for jj = 1:length(INDS)
    IND = INDS(jj);
    DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
           [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];


    ALONGTRACK_RANGE = 1:size(DATA,2);
    if size(DATA,2) == 807
        % Do not update the along-track data (already the right size)
        SEGS = 3;
    elseif size(DATA,2) == 1142
        ALONGTRACK_RANGE = ALONGTRACK_RANGE(1:570);
        SEGS = 3;
    else
        error(' ')
    end
    TRACK = 50;
    
    ALONGTRACK_DISTANCE = abs([SWOT.lon{IND}(TRACK,:) + 1i*SWOT.lat{IND}(TRACK,:)] - [SWOT.lon{IND}(TRACK,1) + 1i*SWOT.lat{IND}(TRACK,1)]);

    % figure
    % plot(LL.LineLon, LL.LineLat, 'k'); hold on
    % pcolor_centered(WC.lon{IND}, WC.lat{IND}, DATA); shading flat
    % colormap(bwr)
    % axis equal
    % set(gca,'YLim',[32 46],'XLim',[-135 -115])
    % clim([-1 1]*0.15); CB = colorbar; CB.Label.String = 'SSHA';

    TS = []; SPEC = [];
    % figure
    % subplot(211)
    % for ii = 3:67
    %     TRACK = ii;
    %     plot(111*ALONGTRACK_DISTANCE(ALONGTRACK_RANGE), DATA(TRACK,ALONGTRACK_RANGE)); hold on
    %     % plot(DATA(TRACK,ALONGTRACK_RANGE)); hold on
    %     TS(:,ii) = DATA(TRACK,ALONGTRACK_RANGE);
    % end
    % plot(111*ALONGTRACK_DISTANCE(ALONGTRACK_RANGE), mean(TS',"omitmissing")', 'k','LineWidth',5 )
    % xlabel('alongtrack distance (km)'); ylabel('SSHA')
    % subplot(212)
    for ii = 3:67
        TRACK = ii;
        [Spec,Freq,~] = nanspectrum(DATA(TRACK,ALONGTRACK_RANGE),dx,'km',SEGS,'-',false,0,'hanning');
        SPEC(:,ii) = Spec;
        % loglog(Freq,Spec,'-','Color',[1 1 1]*0.7); hold on
    end
    SPEC_ = mean(SPEC',"omitmissing")';
    SPEC_MEAN_track(:,jj) = SPEC_;
    
    % loglog(Freq,SPEC_,'k','LineWidth',5)
    % % loglog(Freq, 0.1 * [Spec(1)/[Freq(1)^SLOPE]] * Freq.^SLOPE, 'k--')
    % loglog( Freq, 0.1 * [SPEC_(1)/[Freq(1)^SLOPE]] * Freq.^SLOPE, 'r', 'LineWidth',5)
    % text(Freq(1), 0.1 * [SPEC_(1)/[Freq(1)^SLOPE]] * Freq(1).^SLOPE, ['\propto k^{' num2str(SLOPE) '}'], 'FontSize',20)

    % SPEC_MEAN = SPEC_MEAN + SPEC_/length(INDS);

    disp(jj)
end

SPEC_MEAN   = mean(  SPEC_MEAN_track,2);
SPEC_MEDIAN = median(SPEC_MEAN_track,2);

% figure
loglog(Freq,SPEC_MEAN_track,'-','Color',[1 1 1]*0.7); hold on
loglog(Freq,SPEC_MEAN,'k','LineWidth',3); hold on
loglog(Freq,SPEC_MEDIAN,'k:','LineWidth',3); hold on
loglog( Freq, 0.1 * [SPEC_(1)/[Freq(1)^SLOPE]] * Freq.^SLOPE, 'r', 'LineWidth',5)

%% Half of a Presentation Figure

close all

figure('Color',[1 1 1])
% % % Range in which the moorings observe (subject to modification)
fill([1/100 1/100 1/10 1/10],10.^[-6 0 0 -6],   [1 1 1]*0.9,'EdgeColor','none'); hold on
set(gca,'YScale','log','XScale','log')

% % % Average spectrum for each track
loglog(Freq,SPEC_MEAN_track,'-','Color',[1 1 1]*0.7,'HandleVisibility','off'); hold on
% for ii=1:size(SPEC_MEAN_track,2) % find which tracks are problematic
%     loglog(Freq,SPEC_MEAN_track(:,ii),'-','Color',[1 1 1]*0.7,'HandleVisibility','off'); hold on
%     text(Freq(end),SPEC_MEAN_track(end,ii),num2str(INDS(ii)))
% end

% % % Median spectrum, slope line, text, and the rest
loglog(nan,nan,'Color',[1 1 1]*0.5)
loglog(Freq,SPEC_MEDIAN,'k','LineWidth',3)
SLOPE = -2;
loglog( Freq(dsearchn(Freq,0.01):end), ...
        10^-6.7 * Freq(dsearchn(Freq,0.01):end).^SLOPE, 'k--', 'LineWidth',1)
text(Freq(1)*10,10^-4.2,'\propto k^{-2}','FontSize',16)
xlim(Freq([1 end-1]))
set(gca,'FontSize',14)
xlabel('Wavenumber (km^{-1})')
ylabel('Spectral power density (m^2 cpkm^{-1})')
legend('Mooring observation range','Track mean spectra','Median spectrum')
title('Along-track Spectra: Cal/Val Orbit, US West Coast Crossover')
set(gcf,'Position',[-1439 271 654 540])

%% Calculate mean ascending and descending tracks in order to remove the geoid per JW's advice

INDS = 1:length(SWOT.ssha_karin_2);
JJ1 = 1; JJ2 = 1;
MeanAscendingTrack = [];
MeanDescendingTrack = [];
for jj = 1:length(INDS)
    IND = INDS(jj)
    DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
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

%% FIG 2: Plot all tracks and spectra of all tracks over multiple passes, ascending and descending (non-tiled, segmented)

close all

% Set the set of wavenumbers to probe in advance:
IND = 180;
DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
        [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];
[~,Freq_predefined,~] = nanspectrum1(DATA(20,1:1140),dx,'km',4,'-k',false,0,'hanning');
% % % % % % % % % % % % % % % % % % % % % % % 

TRACK_LENGTH = nan(length(SWOT.ssha_karin_2),1);
for ii=1:length(TRACK_LENGTH)
    TRACK_LENGTH(ii) = size(SWOT.ssha_karin_2{ii},2);
end

clear SPEC_MEDIAN FREQ
% ASCENDING_DESCENDING = [1142]; % 1142 = ascending, 807 = descending
II = 1;

figure('Color',[1 1 1])
% fill([1/300 1/300 1/2 1/2],10.^[-6 1 1 -6],   [1 1 1]*0.95,'EdgeColor','none'); hold on % glider range
fill([1/100 1/100 1/20 1/20],10.^[-6 1 1 -6],   [1 1 1]*0.9,'EdgeColor','none'); hold on % mooring range
% fill([1/71 1/71 1/69 1/69],10.^[-6 1 1 -6],   [1 1 1]*0.6,'EdgeColor','none','HandleVisibility','off'); hold on % ~70km

% figure('Color',[1 1 1]) % test to be sure correct tracks are being analyzed %$


INDS = 1:length(SWOT.ssha_karin_2);
% INDS = INDS(TRACK_LENGTH == A_D); % ascending (1142) or descending (807)
SEGS = 3;
SPEC_MEAN_track = [];

for jj = 1:length(INDS)
    IND = INDS(jj);
    DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
        [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];

    ALONGTRACK_RANGE = 1:size(DATA,2);
    ASCENDING_DESCENDING = sign(diff((SWOT.lat{jj}(50,1:2))));
    if ASCENDING_DESCENDING == 1 % ascending
        % Shorten only up to coast
        SEGS = 2;
        ALONGTRACK_RANGE = ALONGTRACK_RANGE(1:570);
        DATA = DATA - MeanAscendingTrack;
        % figure(2);m_pcolor_centered(SWOT.lon{jj}(:,ALONGTRACK_RANGE),SWOT.lat{jj}(:,ALONGTRACK_RANGE),DATA(:,ALONGTRACK_RANGE));figure(1)%$
    elseif ASCENDING_DESCENDING == -1 % descending
        SEGS = 4;
        ALONGTRACK_RANGE = ALONGTRACK_RANGE(1:1140);
        DATA = DATA - MeanDescendingTrack;
        % figure(2);m_pcolor_centered(SWOT.lon{jj}(:,ALONGTRACK_RANGE),SWOT.lat{jj}(:,ALONGTRACK_RANGE),DATA(:,ALONGTRACK_RANGE));figure(1)%$
    else
        error(' ')
    end
    TRACK = 50;

    ALONGTRACK_DISTANCE = abs([SWOT.lon{IND}(TRACK,:) + 1i*SWOT.lat{IND}(TRACK,:)] - [SWOT.lon{IND}(TRACK,1) + 1i*SWOT.lat{IND}(TRACK,1)]);

    TS = []; SPEC = []; DATA_COVERAGE = [];
    % for ii = 3:67
    %     TRACK = ii;
    %     [Spec,Freq,~] = nanspectrum1(DATA(TRACK,ALONGTRACK_RANGE),dx,'km',SEGS,'-',false,0,'hanning');
    %     SPEC(:,ii) = Spec;
    % end
    for ii = 3:67
        TRACK = ii;
        [Spec,Freq,~] = nanspectrum1(DATA(TRACK,ALONGTRACK_RANGE),dx,'km',SEGS,'-k',false,0,'hanning');

        if sum(isfinite(DATA(TRACK,ALONGTRACK_RANGE))) > 100
            [Spec,Freq,~] = nunanspectrum(DATA(TRACK,ALONGTRACK_RANGE), ALONGTRACK_RANGE*dx, 'km', ...
                'Segments',SEGS,'Plot_option','-','Plot',false,'Window','hanning','Freq',Freq_predefined);
        else % all nan, which is hard to deal with using nunanspectrum
            [Spec,Freq,~] = nanspectrum1(DATA(TRACK,ALONGTRACK_RANGE),dx,'km',SEGS,'-k',false,0,'hanning');
            Spec = nan*Spec;
        end

        SPEC(:,ii) = Spec;
        DATA_COVERAGE(ii,jj) = sum(isfinite(DATA(TRACK,ALONGTRACK_RANGE)))/length(DATA(TRACK,ALONGTRACK_RANGE));
    end
    SPEC_ = mean(SPEC',"omitmissing")';
    SPEC_MEAN_track(:,jj) = SPEC_;

    % figure(2);m_pcolor_centered(SWOT.lon{jj},SWOT.lat{jj},DATA);pause(0.2);figure(1)%$
    disp(jj)
end

SPEC_MEAN_track = SPEC_MEAN_track .* SPEC_MEAN_track./SPEC_MEAN_track; % turn 0 to nan

SPEC_MEAN   = mean(  SPEC_MEAN_track,2,'omitnan');
SPEC_MEDIAN{II} = median(SPEC_MEAN_track,2,'omitnan');
FREQ{II} = Freq;

% % Presentation Figure

% % % Range in which the moorings observe (subject to modification)
% fill([1/100 1/100 1/10 1/10],10.^[-6 1 1 -6],   [1 1 1]*0.9,'EdgeColor','none','HandleVisibility','off'); hold on
set(gca,'YScale','log','XScale','log')

% % % Average spectrum for each track
loglog(Freq,SPEC_MEAN_track,'-','Color',[1 1 1]*0.7,'HandleVisibility','off'); hold on
% for ii=1:size(SPEC_MEAN_track,2) % find which tracks are problematic
%     loglog(Freq,SPEC_MEAN_track(:,ii),'-','Color',[1 1 1]*0.7,'HandleVisibility','off'); hold on
%     text(Freq(end),SPEC_MEAN_track(end,ii),num2str(INDS(ii)))
% end

% % % Median spectrum, slope line, text, and the rest
SLOPE = -2;
xlim(Freq([1 end-1]))
ylim(10.^[-6 1])
set(gca,'FontSize',14)
xlabel('Wavenumber (cpkm)') % xlabel('Wavenumber (km^{-1})')
ylabel('Spectral power density (m^2 cpkm^{-1})')



loglog(nan,nan,'Color',[1 1 1]*0.5)

loglog( Freq(dsearchn(Freq,0.01):end), ...
        10^-6.7 * Freq(dsearchn(Freq,0.01):end).^-2, 'k--', 'LineWidth',3,'HandleVisibility','off')
text(Freq(1)*10,10^-3.9,'\propto k^{-2}','FontSize',16,'HandleVisibility','off')


loglog( Freq(1:dsearchn(Freq,0.02)), ...
        10^-7.5 * Freq(1:dsearchn(Freq,0.02)).^[-3], 'k--', 'LineWidth',3,'HandleVisibility','off')
text(Freq(1)*12,10^-2.1,'\propto k^{-3}','FontSize',16,'HandleVisibility','off')

loglog(FREQ{1},SPEC_MEAN,'LineWidth',3,'Color',[24;116;205]'/255)
% loglog(FREQ{1},SPEC_MEDIAN{1},'LineWidth',3,'Color',[24;116;205]'/255)

% legend('Glider observation range','Mooring observation range','Track mean spectra','Median spectrum (ascending)','Median spectrum (descending)')
legend('Mooring observation scales','Track-mean spectra','Mean spectrum')
% title({'Along-track Spectra: Cal/Val Orbit';'US West Coast Crossover'},'FontWeight','bold','FontSize',16)
set(gcf,'Position',[-1776         228         813         730])
ylim(10.^[-5 0])

fontsize(24,"point")

% % 70 km line label:
% text(1/67,10^-4.8,'(70 km)^{-1}','FontSize',20,'HandleVisibility','off')

%% Save

% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/SWOT_PowerSpectrum_MeanRemoved.pdf',...
% 'BackgroundColor','none','ContentType','vector')

% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/F2.pdf',...
% 'BackgroundColor','none','ContentType','vector')

%% 2D Wavenumber spectra (not for publication, just investigation)


KL_SPEC = [];
JJ = 1;
for jj = 1:length(INDS)
    IND = INDS(jj);
    DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
        [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];

    ALONGTRACK_RANGE = 1:size(DATA,2);
    ASCENDING_DESCENDING = sign(diff((SWOT.lat{jj}(50,1:2))));
    if ASCENDING_DESCENDING == 1 % ascending
        % Shorten only up to coast
        SEGS = 2;
        ALONGTRACK_RANGE = ALONGTRACK_RANGE(1:570);
        DATA = DATA - MeanAscendingTrack;
        % figure(2);m_pcolor_centered(SWOT.lon{jj}(:,ALONGTRACK_RANGE),SWOT.lat{jj}(:,ALONGTRACK_RANGE),DATA(:,ALONGTRACK_RANGE));figure(1)%$
    elseif ASCENDING_DESCENDING == -1 % descending
        SEGS = 4;
        ALONGTRACK_RANGE = ALONGTRACK_RANGE(1:1140);
        DATA = DATA - MeanDescendingTrack;
        % figure(2);m_pcolor_centered(SWOT.lon{jj}(:,ALONGTRACK_RANGE),SWOT.lat{jj}(:,ALONGTRACK_RANGE),DATA(:,ALONGTRACK_RANGE));figure(1)%$
    else
        error(' ')
    end

    TS = []; SPEC = [];
    %%% All:
    % [KL_Spec,KL_Freq,~] = nanspectrum2(DATA(:,ALONGTRACK_RANGE),...
    %                                    [dx,dx],...
    %                                    [SEGS,3],...
    %                                    {'alongtrack km','crosstrack km'},...
    %                                    false, 'hanning');

    %%% Top swath
    [KL_Spec,KL_Freq,~] = nanspectrum2(DATA(5:30,ALONGTRACK_RANGE),...
                                       [dx,dx],...
                                       [SEGS,1],...
                                       {'alongtrack km','crosstrack km'},...
                                       false, 'hanning');

    % %%% Top swath
    % [KL_Spec,KL_Freq,~] = nanspectrum2(DATA(39:65,ALONGTRACK_RANGE),...
    %                                    [dx,dx],...
    %                                    [SEGS,1],...
    %                                    {'alongtrack km','crosstrack km'},...
    %                                    false, 'hanning');

    KL_SPEC(:,:,JJ) = KL_Spec;
    disp(jj)
    JJ = JJ + 1;
end
%%
close all
figure('Color','w')
pcolor_centered(KL_Freq(:,:,1), KL_Freq(:,:,2), log10( mean(KL_SPEC,3,'omitnan') ) )
% pcolor_centered(KL_Freq(:,:,1), KL_Freq(:,:,2), log10( KL_SPEC ) )
clim([prctile(log10(KL_SPEC(KL_SPEC>0)),1) max(log10(KL_SPEC(:)))])
colorbar
shading flat
xlabel('Alongtrack (1/km)')
ylabel('Crosstrack (1/km)')

%% Structure functions from full tracks, every cell

close all

error(['This calculation takes a while'])

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
%%
SF_MAT_swot = nan(size(SF_ARRAY_swot,2),size(SF_ARRAY_swot,1)*size(SF_ARRAY_swot,3));
for kk = 1:size(SF_ARRAY_swot,3)
    SF_MAT_swot(:,[1:size(SF_ARRAY_swot,1)] + size(SF_ARRAY_swot,1)*[kk-1]) = SF_ARRAY_swot(:,:,kk)';
end

%% Plot full structure functions

close all

figure('Color',[1 1 1])
for jj = 1:round(0.1*size(SF_MAT_swot,2)) % 1:size(SF_MAT_swot,2)
    plot(DIST_VEC,SF_MAT_swot(:,jj)*100^2,'-','Color',[1 1 1]*0.7); hold on
end

fill([flip(DIST_VEC); DIST_VEC], ...
     [flip(mean(SF_MAT_swot,2,'omitnan') + std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))); ...
           mean(SF_MAT_swot,2,'omitnan') - std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))]*100^2, ...
     [1 1 1]*0.7); hold on
plot(DIST_VEC, mean(SF_MAT_swot,2,'omitnan')*100^2, 'k.-'); hold on
xlabel('Distance (km)')
ylabel('SF (cm^2)')
set(gca,'FontSize',24)

% set(gca,'xscale','log','yscale','log')

% plot(DIST_VEC, mean(SF_MAT_swot,2,'omitnan') + std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2)), '.-')
% plot(DIST_VEC, mean(SF_MAT_swot,2,'omitnan') - std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2)), '.-')
% plot(DIST_VEC, ones(size(mean(SF_MAT_swot,2,'omitnan') + std(SF_MAT_swot,0,2,'omitnan')))./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2)), '.-')
% plot(DIST_VEC, ones(size(mean(SF_MAT_swot,2,'omitnan') + std(SF_MAT_swot,0,2,'omitnan'))).*(sum(isfinite(SF_MAT_swot(:,:)), 2)), '.-')


%% NLLSF for polynomials

close all

Mean_SF_MAT_swot = mean(SF_MAT_swot,2,'omitnan')*100^2;

FitBasisFunction = {'a*t.^b','a','b'};
FitParameters0 = [1,1];
FitBasisDerivatives = {'t.^b','a*log(t).*t.^b'};
ToleranceLevel = 0.001; % NLLSF parameter, tells change in parameter magnitude, if small then the result was stable
N_nllsf_iterations = 10; % number of iterations before we give up

[NLLSF_COEF,~,~,~] = ...
            nonlinear_lsqf(Mean_SF_MAT_swot(1:20), DIST_VEC(1:20), ...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);


%%
%%
%%
%%
%% Look at slopes from along-track data and Welch's method

close all

II = 1;

INDS = 1:length(SWOT.ssha_karin_2);
% INDS = INDS(TRACK_LENGTH == A_D); % ascending (1142) or descending (807)
SEGS = 3;
SPEC_MEAN_track = [];
SPEC = [];
DATA_COVERAGE = [];

% Set the set of wavenumbers to probe in advance:
IND = 180;
DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
        [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];
[~,Freq_predefined,~] = nanspectrum1(DATA(20,1:1140),dx,'km',4,'-k',false,0,'hanning');
% % % % % % % % % % % % % % % % % % % % % % % 

for jj = 1:length(INDS)
    IND = INDS(jj);
    DATA = [SWOT.ssha_karin_2{IND} + SWOT.internal_tide_hret{IND} + SWOT.height_cor_xover{IND}].*...
        [~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})./~(SWOT.ssha_karin_2_qual{IND} + SWOT.height_cor_xover_qual{IND})];

    ALONGTRACK_RANGE = 1:size(DATA,2);
    ASCENDING_DESCENDING = sign(diff((SWOT.lat{jj}(50,1:2))));
    if ASCENDING_DESCENDING == 1 % ascending
        % Shorten only up to coast
        SEGS = 2;
        ALONGTRACK_RANGE = ALONGTRACK_RANGE(1:570);
        DATA = DATA - MeanAscendingTrack;
        % figure(2);m_pcolor_centered(SWOT.lon{jj}(:,ALONGTRACK_RANGE),SWOT.lat{jj}(:,ALONGTRACK_RANGE),DATA(:,ALONGTRACK_RANGE));figure(1)%$
    elseif ASCENDING_DESCENDING == -1 % descending
        SEGS = 4;
        ALONGTRACK_RANGE = ALONGTRACK_RANGE(1:1140);
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
        [Spec,Freq,~] = nanspectrum1(DATA(TRACK,ALONGTRACK_RANGE),dx,'km',SEGS,'-k',false,0,'hanning');
        
        if sum(isfinite(DATA(TRACK,ALONGTRACK_RANGE))) > 100
            FFT_vs_NUFFT = [FFT_vs_NUFFT; 1];
            [Spec,Freq,~] = nunanspectrum(DATA(TRACK,ALONGTRACK_RANGE), ALONGTRACK_RANGE*dx, 'km', ...
                'Segments',SEGS,'Plot_option','-','Plot',false,'Window','hanning','Freq',Freq_predefined);
        else % all nan, which is hard to deal with using nunanspectrum
            FFT_vs_NUFFT = [FFT_vs_NUFFT; 0];
            [Spec,Freq,~] = nanspectrum1(DATA(TRACK,ALONGTRACK_RANGE),dx,'km',SEGS,'-k',false,0,'hanning');
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

ApproxCrossTrackDist = 111*cosd(SWOT.lat{1}(1,500))*abs([SWOT.lon{1}(:,500) + 1i*SWOT.lat{1}(:,500)] - [SWOT.lon{1}(35,500) + 1i*SWOT.lat{1}(35,500)]) ...
                                                .* sign([SWOT.lon{1}(:,500)] -  [SWOT.lon{1}(35,500)]);

figure('Color',[1 1 1])
fill([1/100 1/100 1/20 1/20],10.^[-6 1 1 -6],   [1 1 1]*0.9,'EdgeColor','none'); hold on % mooring range

% % % Range in which the moorings observe (subject to modification)
set(gca,'YScale','log','XScale','log')

% % % Median spectrum, slope line, text, and the rest
SLOPE = -2;
xlim(Freq([1 end-1]))
ylim(10.^[-6 1])
set(gca,'FontSize',14)
xlabel('Wavenumber (km^{-1})')
ylabel('Spectral power density (m^2 cpkm^{-1})')

% loglog(nan,nan,'Color',[1 1 1]*0.5)

% % % 
loglog( Freq(dsearchn(Freq,0.01):end), ...
        10^-6.7 * Freq(dsearchn(Freq,0.01):end).^-2, 'k--', 'LineWidth',3,'HandleVisibility','off')
text(Freq(1)*10,10^-3.9,'\propto k^{-2}','FontSize',16,'HandleVisibility','off')
loglog( Freq(1:dsearchn(Freq,0.02)), ...
        10^-7.5 * Freq(1:dsearchn(Freq,0.02)).^[-3], 'k--', 'LineWidth',3,'HandleVisibility','off')
text(Freq(1)*12,10^-2.1,'\propto k^{-3}','FontSize',16,'HandleVisibility','off')
% % % 

% loglog(Freq,SPEC_,'-','LineWidth',3,'Color',[24;116;205]'/255)
% loglog(Freq,SPEC_,'-','Color',[1 1 1]*.5)
% TURBO = turbo(size(SPEC_,2));
TURBO = [flip(turbo(34)) ; 0 0 0 ; turbo(34)];
for ii = 1:size(SPEC_,2)
    loglog(Freq,SPEC_(:,ii),'-','Color',TURBO(ii,:),'LineWidth',1)
end

% legend('Mooring observation scales','Track-mean spectra','Mean spectrum')
% title({'Along-track Spectra: Cal/Val Orbit';'US West Coast Crossover'},'FontWeight','bold','FontSize',16)
set(gcf,'Position',[-1776         228         813         730])
ylim(10.^[-5 0])

fontsize(24,"point")

figure
for ii = 1:size(SPEC_,2)
    scatter(ApproxCrossTrackDist(ii),0,20,TURBO(ii,:),'filled'); hold on
end

% Slope of each spectrum past some frequency (Freq is in cpkm)
KM_LONGEST = 100;
ind_min = dsearchn(Freq,1/[KM_LONGEST]);
ind_max = dsearchn(Freq,1/[20]);
% ind_max = length(Freq); % maximum
HH = [ones(length(Freq(ind_min:ind_max)),1), log10(Freq(ind_min:ind_max))];
SPEC_COEF = [];
figure % verify that we are cutting it off as expected
for ii = 1:size(SPEC_,2)
    SPEC_COEF(ii,:) = [HH\log10(SPEC_(ind_min:ind_max,ii))]';
    % FITTED_SPEC = 10.^( HH*SPEC_COEF(ii,:)' );
    % loglog(Freq(ind_min:ind_max),FITTED_SPEC, 'k'); hold on
    loglog(Freq(ind_min:ind_max), [10.^SPEC_COEF(ii,1) * Freq(ind_min:ind_max).^SPEC_COEF(ii,2)], 'k');hold on
    loglog(Freq(ind_min:ind_max),SPEC_(ind_min:ind_max,ii),'-','Color',TURBO(ii,:),'LineWidth',1)
end
% %%
% figure
% for ii = 1:size(SPEC_,2)
%     % loglog(Freq(ind_min:ind_max), 10.^SPEC_COEF(ii,1) + Freq(ind_min:ind_max).^SPEC_COEF(ii,2));hold on
%     loglog(Freq(ind_min:ind_max), [10.^SPEC_COEF(ii,1) * Freq(ind_min:ind_max).^SPEC_COEF(ii,2)], 'k');hold on
% end
% %
figure('Color','w')
plot(ApproxCrossTrackDist, SPEC_COEF(:,2) , '.-k', 'LineWidth',1,'MarkerSize',20); hold on
xlabel('Approx. distance from nadir')
ylabel(['Spectral slope of average spectrum'])

fontsize(24,"point")

%% Figure with multiple length scale regimes

close all

figure % verify that we are cutting it off as expected
for ii = 1:size(SPEC_,2)
    loglog(Freq,SPEC_(:,ii),'-','Color',TURBO(ii,:),'LineWidth',1); hold on
end
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
        loglog(Freq(ind_min:ind_max), [10.^SPEC_COEF(ii,1) * Freq(ind_min:ind_max).^SPEC_COEF(ii,2)], 'k');hold on
    end
    CROSSTRACK_SLOPES(:,jj) = SPEC_COEF(:,2);
end
figure('Color','w')

plot(ApproxCrossTrackDist, CROSSTRACK_SLOPES , '.-', 'LineWidth',1,'MarkerSize',20); hold on
LEG = {};
for ii = 1:length(KM_LONGEST)
    LEG{ii} = [num2str(KM_LONGEST(ii)) ' to ' num2str(KM_SHORTEST(ii)) ' km'];
end
legend(LEG)
xlabel('Approx. distance from nadir')
ylabel(['Spectral slope of average spectrum'])

fontsize(24,"point")

% %% Integrated spectrum
% 
% close all
% 
% figure('Color','w')
% plot(ApproxCrossTrackDist, sum(SPEC_,1,'omitnan'),'.-'); hold on
% % plot(ApproxCrossTrackDist, sum(SPEC_(1,:),1,'omitnan'),'.-')
% % plot(ApproxCrossTrackDist, sum(SPEC_(2,:),1,'omitnan'),'.-')
% % plot(ApproxCrossTrackDist, sum(SPEC_(3,:),1,'omitnan'),'.-')
% % plot(ApproxCrossTrackDist, sum(SPEC_(4,:),1,'omitnan'),'.-')
% plot(ApproxCrossTrackDist, sum(SPEC_(10:end,:),1,'omitnan'),'.-')
% 
% 
% xlabel('Approximate distance from center of swath')
% ylabel('\Sigma S(k)')
% fontsize(18,"point")

%% Subtract the "noise spectrum" by subtracting the PSD (averaged over 1/8 to 1/4 cpkm)
% from the spectra as white noise and then getting the slope. As JW wrote:

% You can use the psd averaged over 1/8 and 1/4 cpkm as the noise floor,
% subtract that from the spectrum, you will get a steeper spectrum, use
% which and redo the spectrum slope as a function of crossswath distance.

% % Demonstration of the assumption:
% figure
% loglog(Freq, Freq.^-3, '.-'); hold on
% loglog(Freq, 50*Freq.^-2, '.-')
% loglog(Freq, Freq.^-3 + 50*Freq.^-2, '.-')

close all

figure % verify that we are cutting it off as expected
KM_LONGEST  = [200 100 60];
KM_SHORTEST = [60  20  8 ];
CROSSTRACK_SLOPES_NFR = nan(size(SPEC_,2),length(KM_LONGEST));
% ^ NFR = "noise floor removed"
for jj = 1:length(KM_LONGEST)
    ind_min = dsearchn(Freq,1/KM_LONGEST(jj));
    ind_max = dsearchn(Freq,1/KM_SHORTEST(jj));
    % ind_max = length(Freq); % maximum
    HH = [ones(length(Freq(ind_min:ind_max)),1), log10(Freq(ind_min:ind_max))];
    SPEC_COEF = [];
    for ii = 1:size(SPEC_,2)
        NOISE_FLOOR = mean(SPEC_(dsearchn(Freq,1/8):end,ii));
        SPEC_COEF(ii,:) = [HH\log10(SPEC_(ind_min:ind_max,ii) - NOISE_FLOOR)]';

        loglog(Freq,SPEC_(:,ii) - NOISE_FLOOR,'-','Color',TURBO(ii,:),'LineWidth',1); hold on
        loglog(Freq(ind_min:ind_max), [10.^SPEC_COEF(ii,1) * Freq(ind_min:ind_max).^SPEC_COEF(ii,2)], 'k');hold on
    end
    CROSSTRACK_SLOPES_NFR(:,jj) = SPEC_COEF(:,2);
end

figure('Color','w')
plot(ApproxCrossTrackDist, CROSSTRACK_SLOPES_NFR , '.-', 'LineWidth',1,'MarkerSize',20); hold on
LEG_NFR = {};
for ii = 1:length(KM_LONGEST)
    LEG_NFR{ii} = [num2str(KM_LONGEST(ii)) ' to ' num2str(KM_SHORTEST(ii)) ' km - noise floor removed'];
end
legend(LEG_NFR)
xlabel('Approx. distance from nadir')
ylabel(['Spectral slope of average spectrum'])

% With the previous results (before removing the noise floor):
figure('Color','w')
% LINE_COLOR = colororder;
LINE_COLOR = [.1 .1 .8;   .1 .8 .1;   .8 .1 .1];
for ii = 1:size(CROSSTRACK_SLOPES,2)
    plot(ApproxCrossTrackDist, CROSSTRACK_SLOPES(:,ii)     , '.-', 'Color', LINE_COLOR(ii,:), 'LineWidth',1,'MarkerSize',20); hold on
end
for ii = 1:size(CROSSTRACK_SLOPES_NFR,2)
    plot(ApproxCrossTrackDist, CROSSTRACK_SLOPES_NFR(:,ii) , 'o-', 'Color', LINE_COLOR(ii,:), 'LineWidth',1,'MarkerSize',10); hold on
end
% for ii = 1:size(CROSSTRACK_SLOPES,2)
%     plot(ApproxCrossTrackDist, CROSSTRACK_SLOPES_NFR(:,ii) - CROSSTRACK_SLOPES(:,ii) , 's-', 'Color', LINE_COLOR(ii,:), 'LineWidth',1,'MarkerSize',10); hold on
% end
LEG_NFR = {};
for ii = 1:length(KM_LONGEST)
    LEG_NFR{ii} = [num2str(KM_LONGEST(ii)) ' to ' num2str(KM_SHORTEST(ii)) ' km - noise floor removed'];
end
% legend([LEG, LEG_NFR])
legend([LEG])
xlabel('Approx. distance from nadir')
ylabel(['Spectral slope of average spectrum'])

fontsize(24,"point")

%% Better approach: Slope(TimeMean(SPEC - N.F.))

close all

ApproxCrossTrackDist = flip(ApproxCrossTrackDist);

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
        SPEC_jj_ = nan(size(squeeze(SPEC(ind_min:ind_max,jj,:))));
        for kk = 1:size(SPEC,3)
                NOISE_FLOOR = mean(SPEC(dsearchn(Freq,1/8):end,jj,kk),'omitnan');
                SPEC_jj_(:,kk) = SPEC(ind_min:ind_max,jj,kk) - NOISE_FLOOR;
        end
        SPEC_jj = mean(SPEC_jj_,2,'omitnan');
        COEF_jj = [HH(SPEC_jj>0,:)\log10(SPEC_jj(SPEC_jj>0))]';
        SPEC_COEF(1,jj) = COEF_jj(1); SPEC_COEF(2,jj) = COEF_jj(2);
    end
    CROSSTRACK_SLOPES_NFR(:,nn) = SPEC_COEF(2,:);
end
CROSSTRACK_SLOPES_NFR(CROSSTRACK_SLOPES_NFR==0) = NaN;

figure('Color','w')
plot(ApproxCrossTrackDist, mean(squeeze(CROSSTRACK_SLOPES_NFR(1,:,:)),2,'omitnan') , '.-', 'LineWidth',1,'MarkerSize',20); hold on
LEG_NFR = {};
for nn = 1:length(KM_LONGEST)
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
CrossTrackPlotIndices = [1:29 32:(length(ApproxCrossTrackDist)-4)];
% CrossTrackPlotIndices = [1:length(ApproxCrossTrackDist)];
for nn = 1:size(CROSSTRACK_SLOPES,2)
    % plot(ApproxCrossTrackDist(CrossTrackPlotIndices), CROSSTRACK_SLOPES(CrossTrackPlotIndices,nn).*[CROSSTRACK_SLOPES(CrossTrackPlotIndices,nn)./CROSSTRACK_SLOPES(CrossTrackPlotIndices,nn)], ...
    %     '.-', 'Color', LINE_COLOR(nn,:), 'LineWidth',1,'MarkerSize',20); hold on
    plot(ApproxCrossTrackDist(CrossTrackPlotIndices), CROSSTRACK_SLOPES(CrossTrackPlotIndices,nn), ...
        '.-', 'Color', LINE_COLOR(nn,:), 'LineWidth',1,'MarkerSize',20); hold on
end
for nn = 1:size(CROSSTRACK_SLOPES_NFR,2)
    % plot(ApproxCrossTrackDist(CrossTrackPlotIndices), CROSSTRACK_SLOPES_NFR(CrossTrackPlotIndices,nn).*[CROSSTRACK_SLOPES_NFR(CrossTrackPlotIndices,nn)./CROSSTRACK_SLOPES_NFR(CrossTrackPlotIndices,nn)], ...
    %     'o-', 'Color', LINE_COLOR(nn,:), 'LineWidth',1,'MarkerSize',10); hold on
    plot(ApproxCrossTrackDist(CrossTrackPlotIndices), CROSSTRACK_SLOPES_NFR(CrossTrackPlotIndices,nn), ...
        'o-', 'Color', LINE_COLOR(nn,:), 'LineWidth',1,'MarkerSize',10); hold on
end
grid on
grid minor

LEG_NFR = {};
for nn = 1:length(KM_LONGEST)
    LEG_NFR{nn} = [num2str(KM_LONGEST(nn)) ' to ' num2str(KM_SHORTEST(nn)) ' km - noise floor removed'];
end
% legend([LEG, LEG_NFR])
LEG_HANDLE = legend([LEG]);
set(gca,'XLim',[-70 70])
xlabel('Distance from nadir (km)')
ylabel(['Spectral slope of average spectrum'])

fontsize(24,"point")
LEG_HANDLE.FontSize = 16;
% LEG_HANDLE.Position = [0.6451 0.7633 0.1875 0.1378];
% LEG_HANDLE.Position = [0.5816    0.7726    0.1875    0.1378];
LEG_HANDLE.Position = [0.5983    0.7726    0.1875    0.1378];


%% Discarded work

% close all
% 
% figure % verify that we are cutting it off as expected
% for ii = 10%1:size(SPEC_,2)
%     loglog(Freq,SPEC_(:,ii),'-','Color',TURBO(ii,:),'LineWidth',1); hold on
% end
% KM_LONGEST  = [8 100];
% KM_SHORTEST = [4 20 ];
% CROSSTRACK_COEFS = nan(size(SPEC_,2),length(KM_LONGEST),2);
% for jj = 1:length(KM_LONGEST)
%     ind_min = dsearchn(Freq,1/KM_LONGEST(jj));
%     ind_max = dsearchn(Freq,1/KM_SHORTEST(jj));
%     % ind_max = length(Freq); % maximum
%     HH = [ones(length(Freq(ind_min:ind_max)),1), log10(Freq(ind_min:ind_max))];
%     SPEC_COEF = [];
%     for ii = 10%1:size(SPEC_,2)
%         SPEC_COEF(ii,:) = [HH\log10(SPEC_(ind_min:ind_max,ii))]';
%         loglog(Freq(ind_min:ind_max), [10.^SPEC_COEF(ii,1) * Freq(ind_min:ind_max).^SPEC_COEF(ii,2)], 'k');hold on
%     end
%     % CROSSTRACK_COEFS(:,jj,1) = SPEC_COEF(:,1);
%     % CROSSTRACK_COEFS(:,jj,2) = SPEC_COEF(:,2);
% end
% figure('Color','w')
% 
% plot(ApproxCrossTrackDist, CROSSTRACK_SLOPES , '.-', 'LineWidth',1,'MarkerSize',20); hold on
% LEG = {};
% for ii = 1:length(KM_LONGEST)
%     LEG{ii} = [num2str(KM_LONGEST(ii)) ' to ' num2str(KM_SHORTEST(ii)) ' km'];
% end
% legend(LEG)
% xlabel('Approx. distance from nadir')
% ylabel(['Spectral slope of average spectrum'])
% 
% fontsize(24,"point")
