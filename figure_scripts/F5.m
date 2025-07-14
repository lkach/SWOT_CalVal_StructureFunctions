%% Figure for SWOT structure functions (FIGURE 5)
%% Load the workspace if needed

% error('Only load if this is needed.')

INPUT = input(['Do you want to load ~530 MB into memory? Enter any number for "yes" or push\n' ...
    '"Enter" with a blank input for "No".\n']);

if isempty(INPUT)
    disp(['No data loaded'])
else
    load(['/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/' ...
        'SF_plot_fixed_bottom_windowed_SWOTt_data.mat'])
    disp('Loaded data')
end

% % % For analyzing Level 3 data, run F2_F3.m up until and including the
% % % section "SWOT_at_S_P_G".

%% Calculate mean ascending and descending tracks in order to remove the geoid per JW's advice

warning(['Also plot the spectra resulting from the removed and unremoved ' ...
         'alongside each other for the supplementary data'])

SWOT.time_mean = nan(length(SWOT.time),1);
for jj = 1:length(SWOT.time_mean)
    SWOT.time_mean(jj) = mean(SWOT.time{jj},'omitnan');
end

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


%% CROSS-TRACK Structure functions from full tracks, every cell

TIME_START = datenum('01-Apr-2023');
TIME_END   = datenum('11-Jul-2023');

TIME_START = datenum('08-May-2023');
TIME_END   = datenum('21-Jun-2023');

xSF_ARRAY_swot = nan(size(SWOT.ssha_karin_2{1},1) - 1, size(SWOT.ssha_karin_2{1},2), length(SWOT.ssha_karin_2));
for jj = dsearchn(SWOT.time_mean,TIME_START):dsearchn(SWOT.time_mean,TIME_END) % 1:length(SWOT.ssha_karin_2)
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

    xSF_MAT_swot_i = nan(size(DATA) - [1 0]);
    for track_i = 1:size(DATA,2)
        TRACK_i = DATA(:,track_i);
        % Consider always omitting some near-nadir/near-edges data at this line.
        LAT_i = SWOT.lat{jj}(:,track_i);
        LON_i = SWOT.lon{jj}(:,track_i);
        xdx = median(m_lldist(LON_i,LAT_i),'omitnan');
        xDIST_VEC = nan(length(LAT_i)-1,1);
        xSF_VEC =   nan(length(LAT_i)-1,1);
        for di = 1:length(xDIST_VEC)
            xDIST_VEC(di) = di*xdx;
            xSF_VEC(di) = mean([TRACK_i(1:[end-di]) - TRACK_i([1+di]:end)].^2,'omitnan');
        end
        xSF_MAT_swot_i(:,track_i) = xSF_VEC;
    end
    xSF_ARRAY_swot(:,:,jj) = xSF_MAT_swot_i;
    disp(jj)
end
xSF_MAT_swot = nan(size(xSF_ARRAY_swot,1),size(xSF_ARRAY_swot,2)*size(xSF_ARRAY_swot,3));
for kk = 1:size(xSF_ARRAY_swot,3)
    xSF_MAT_swot(:,[1:size(xSF_ARRAY_swot,2)] + size(xSF_ARRAY_swot,2)*[kk-1]) = xSF_ARRAY_swot(:,:,kk);
end

%% Along-track structure functions from full tracks, every cell

close all

% error(['This calculation takes a while'])

SF_ARRAY_swot = nan(size(SWOT.ssha_karin_2{1},1), size(SWOT.ssha_karin_2{1},2) - 1, length(SWOT.ssha_karin_2));
for jj = dsearchn(SWOT.time_mean,TIME_START):dsearchn(SWOT.time_mean,TIME_END) % 1:length(SWOT.ssha_karin_2)
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

%% Along-track structure functions from central tracks, every cell

SFcent_ARRAY_swot = nan(size(SWOT.ssha_karin_2{1},1), size(SWOT.ssha_karin_2{1},2) - 1, length(SWOT.ssha_karin_2));
for jj = dsearchn(SWOT.time_mean,TIME_START):dsearchn(SWOT.time_mean,TIME_END) % 1:length(SWOT.ssha_karin_2)
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

    SFcent_MAT_swot_i = nan(size(DATA) - [0 1]);
    for track_i = [16:21 49:54]
        TRACK_i = DATA(track_i,:);
        LAT_i = SWOT.lat{jj}(track_i,:);
        LON_i = SWOT.lon{jj}(track_i,:);
        dx = median(m_lldist(LON_i,LAT_i),'omitnan');
        DIST_VECcent = nan(length(LAT_i)-1,1);
        SFcent_VEC =   nan(length(LAT_i)-1,1);
        for di = 1:length(DIST_VECcent)
            DIST_VECcent(di) = di*dx;
            SFcent_VEC(di) = mean([TRACK_i(1:[end-di]) - TRACK_i([1+di]:end)].^2,'omitnan');
        end
        SFcent_MAT_swot_i(track_i,:) = SFcent_VEC;
    end
    SFcent_ARRAY_swot(:,:,jj) = SFcent_MAT_swot_i;
    disp(jj)
end
% %
SFcent_MAT_swot = nan(size(SFcent_ARRAY_swot,2),size(SFcent_ARRAY_swot,1)*size(SFcent_ARRAY_swot,3));
for kk = 1:size(SFcent_ARRAY_swot,3)
    SFcent_MAT_swot(:,[1:size(SFcent_ARRAY_swot,1)] + size(SFcent_ARRAY_swot,1)*[kk-1]) = SFcent_ARRAY_swot(:,:,kk)';
end

%% Along-track structure functions from edge tracks, every cell

% Copy the above but replace "cent" with "edge" and redefine
% "for track_i = ..."
% Then add the necessary lines below to both panels, copying and
% reformatting the blocks for "cent".

SFedge_ARRAY_swot = nan(size(SWOT.ssha_karin_2{1},1), size(SWOT.ssha_karin_2{1},2) - 1, length(SWOT.ssha_karin_2));
for jj = dsearchn(SWOT.time_mean,TIME_START):dsearchn(SWOT.time_mean,TIME_END) % 1:length(SWOT.ssha_karin_2)
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

    SFedge_MAT_swot_i = nan(size(DATA) - [0 1]);
    for track_i = [5:7 29:31 39:41 63:65]
        TRACK_i = DATA(track_i,:);
        LAT_i = SWOT.lat{jj}(track_i,:);
        LON_i = SWOT.lon{jj}(track_i,:);
        dx = median(m_lldist(LON_i,LAT_i),'omitnan');
        DIST_VECedge = nan(length(LAT_i)-1,1);
        SFedge_VEC =   nan(length(LAT_i)-1,1);
        for di = 1:length(DIST_VECedge)
            DIST_VECedge(di) = di*dx;
            SFedge_VEC(di) = mean([TRACK_i(1:[end-di]) - TRACK_i([1+di]:end)].^2,'omitnan');
        end
        SFedge_MAT_swot_i(track_i,:) = SFedge_VEC;
    end
    SFedge_ARRAY_swot(:,:,jj) = SFedge_MAT_swot_i;
    disp(jj)
end
% %
SFedge_MAT_swot = nan(size(SFedge_ARRAY_swot,2),size(SFedge_ARRAY_swot,1)*size(SFedge_ARRAY_swot,3));
for kk = 1:size(SFedge_ARRAY_swot,3)
    SFedge_MAT_swot(:,[1:size(SFedge_ARRAY_swot,1)] + size(SFedge_ARRAY_swot,1)*[kk-1]) = SFedge_ARRAY_swot(:,:,kk)';
end

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
CENTER_COLOR = [77  191  239]/255; % [ 41, 203, 255 ]./255; % Custom color
EDGE_COLOR   = [218   83   25]/255; % [ 244, 179, 0  ]./255; % Custom color
CROSS_COLOR  = [185, 19, 209 ]/255; % [ 244, 179, 0  ]./255; % Custom color

tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');

nexttile;%subplot(121) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % SWOT_ALL
% LEG_TEXT = {LEG_TEXT{:},'SWOT cal/val all data (29.9N-50.2N)'};
LEG_TEXT = {LEG_TEXT{:},['Along-track SWOT cal/val ' newline 'all data (' num2str(round(min([SWOT.lat{end}(:) ; SWOT.lat{end-1}(:)]),1)) 'N-' ...
                                                   num2str(round(max([SWOT.lat{end}(:) ; SWOT.lat{end-1}(:)]),1)) 'N)']};
fill([flip(DIST_VEC); DIST_VEC], ...
     [flip(mean(SF_MAT_swot,2,'omitnan') + std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))); ...
           mean(SF_MAT_swot,2,'omitnan') - std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))], ...
     [1 1 1]*0.7, 'EdgeAlpha',0,'HandleVisibility','off'); hold on
Mean_SF_MAT_swot = mean(SF_MAT_swot,2,'omitnan')*m2_to_cm2;
STD_SF_MAT_swot = std(SF_MAT_swot,[],2,'omitnan')*m2_to_cm2./sqrt(sum(isfinite(SF_MAT_swot),2));
% plot(DIST_VEC, Mean_SF_MAT_swot, 'k.-', 'MarkerSize',10); hold on
errorbar(DIST_VEC, Mean_SF_MAT_swot,STD_SF_MAT_swot, 'k.-', 'MarkerSize',15); hold on


% % % % SWOT 12 most central tracks (6 per swath)
LEG_TEXT = {LEG_TEXT{:},['Along-track SWOT cal/val ' newline '(12 central tracks)']};
fill([flip(DIST_VECcent); DIST_VECcent], ...
     [flip(mean(SFcent_MAT_swot,2,'omitnan') + std(SFcent_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SFcent_MAT_swot(:,:)), 2))); ...
           mean(SFcent_MAT_swot,2,'omitnan') - std(SFcent_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SFcent_MAT_swot(:,:)), 2))], ...
     [1 1 1]*0.7, 'EdgeAlpha',0,'HandleVisibility','off'); hold on
Mean_SFcent_MAT_swot = mean(SFcent_MAT_swot,2,'omitnan')*m2_to_cm2;
STD_SFcent_MAT_swot = std(SFcent_MAT_swot,[],2,'omitnan')*m2_to_cm2./sqrt(sum(isfinite(SFcent_MAT_swot),2));
% plot(DIST_VECcent, Mean_SFcent_MAT_swot, '.-', 'Color',BLUE); hold on
errorbar(DIST_VECcent, Mean_SFcent_MAT_swot, STD_SFcent_MAT_swot, '.-', 'Color',CENTER_COLOR); hold on

% % % % SWOT 12 edge-most tracks (6 per swath)
LEG_TEXT = {LEG_TEXT{:},['Along-track SWOT cal/val ' newline '(12 edge-most tracks)']};
fill([flip(DIST_VECedge); DIST_VECedge], ...
     [flip(mean(SFedge_MAT_swot,2,'omitnan') + std(SFedge_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SFedge_MAT_swot(:,:)), 2))); ...
           mean(SFedge_MAT_swot,2,'omitnan') - std(SFedge_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SFedge_MAT_swot(:,:)), 2))], ...
     [1 1 1]*0.7, 'EdgeAlpha',0,'HandleVisibility','off'); hold on
Mean_SFedge_MAT_swot = mean(SFedge_MAT_swot,2,'omitnan')*m2_to_cm2;
STD_SFedge_MAT_swot = std(SFedge_MAT_swot,[],2,'omitnan')*m2_to_cm2./sqrt(sum(isfinite(SFedge_MAT_swot),2));
% plot(DIST_VECedge, Mean_SFedge_MAT_swot, '.-', 'Color',ORANGE); hold on
errorbar(DIST_VECedge, Mean_SFedge_MAT_swot, STD_SFedge_MAT_swot, '.-', 'Color',EDGE_COLOR); hold on


% Plot cross-track structure function, but with the last 8 elements omitted
% because the last 7 are NaN and the one before is clearly spurious:
LEG_TEXT = {LEG_TEXT{:},['Cross-track']};
Mean_xSF_MAT_swot = mean(xSF_MAT_swot,2,'omitnan')*m2_to_cm2;
STD_xSF_MAT_swot = std(xSF_MAT_swot,[],2,'omitnan')*m2_to_cm2./sqrt(sum(isfinite(xSF_MAT_swot),2));
% loglog(xDIST_VEC(1:[end-8]), Mean_xSF_MAT_swot(1:[end-8],:), 'r.-'); hold on
errorbar(xDIST_VEC(1:[end-8]), Mean_xSF_MAT_swot(1:[end-8],:), STD_xSF_MAT_swot(1:[end-8],:), '.-','Color',CROSS_COLOR); hold on


% % % SWOT_at_S_P_G
LEG_TEXT = {LEG_TEXT{:},['SWOT cal/val at ' newline 'Moorings+Gliders']};
BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_SWOT_moor_12hr,Separation_SWOT_moor_12hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_SWOT_moor_12hr, StructFunc_SWOT_moor_12hr, '.', 'MarkerSize',10); hold on
errorbar_stdofmean(BIN_CENTR, Separation_SWOT_moor_12hr, m2_to_cm2*StructFunc_SWOT_moor_12hr, BIN_EDGES, [1 1 1]*0.7); hold on
Mean_SWOT_at_S_P_G = interval_avg(Separation_SWOT_moor_12hr, m2_to_cm2*StructFunc_SWOT_moor_12hr, BIN_EDGES);
% plot(BIN_CENTR,Mean_SWOT_at_S_P_G,'r') % verify that this is the quantity of interest

% %% Plot full CROSS-TRACK structure functions
% 
% close all
% 
% figure('Color',[1 1 1])
% for jj = 1:round(0.1*size(xSF_MAT_swot,2)) % 1:size(SF_MAT_swot,2)
%     % plot(DIST_VEC,SF_MAT_swot(:,jj)*100^2,'-','Color',[1 1 1]*0.7); hold on
%     loglog(xDIST_VEC,xSF_MAT_swot(:,jj)*100^2,'-','Color',[1 1 1]*0.7); hold on
% end
% 
% fill([flip(xDIST_VEC); xDIST_VEC], ...
%      [flip(mean(xSF_MAT_swot,2,'omitnan') + std(xSF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(xSF_MAT_swot(:,:)), 2))); ...
%            mean(xSF_MAT_swot,2,'omitnan') - std(xSF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(xSF_MAT_swot(:,:)), 2))]*100^2, ...
%      [1 1 1]*0.7); hold on
% % plot(DIST_VEC, mean(SF_MAT_swot,2,'omitnan')*100^2, 'k.-'); hold on
% loglog(xDIST_VEC, mean(xSF_MAT_swot,2,'omitnan')*100^2, 'k.-'); hold on
% xlabel('Distance (km)')
% ylabel('SF (cm^2)')
% set(gca,'FontSize',24)


xlabel('Separation (km)')
ylabel('$\overline{\big(\mathrm{SSHA}(x) - \mathrm{SSHA}(x+s)\big)^2}$ \quad (cm$^2$)','Interpreter','latex','FontSize',16)

set(gca,'XLim',DIST_VEC([1 end]))
set(gca,'YLim',[0.1 100]); set(gca,'XLim',[0 1000])
set(gca,'yscale',YSCALE,'xscale',XSCALE)
set(gca,'yscale','log','xscale','log')
set(gca,'FontSize',20); grid on



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
    
    % Along-track
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

    % Along-track, 12 central tracks
    % IND_fitting_lags = [4:2:10 12:4:20 25:5:50 60:10:100 150]';
    SWOTcent_fitted_SpecSlopes            = nan(size(IND_fitting_lags));
    SWOTcent_fitted_SpecSlopes_covariance = nan(size(IND_fitting_lags));
    for xi = 1:length(IND_fitting_lags)
        XI = IND_fitting_lags(xi);
        [NLLSF_COEF,JACOBIAN,SF_FIT,ConvergenceRecord] = ...
        nonlinear_lsqf(mean(SFcent_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2, DIST_VECcent(1:XI), ...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations);
        NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ [mean(SFcent_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2] - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
        SWOTcent_fitted_SpecSlopes(xi) = -NLLSF_COEF(2) - 1;
        SWOTcent_fitted_SpecSlopes_covariance(xi) = NLLSF_COEF_cov(2,2);
        disp(ConvergenceRecord(end))
    end

    % Along-track, 12 edge-most tracks
    % IND_fitting_lags = [4:2:10 12:4:20 25:5:50 60:10:100 150]';
    SWOTedge_fitted_SpecSlopes            = nan(size(IND_fitting_lags));
    SWOTedge_fitted_SpecSlopes_covariance = nan(size(IND_fitting_lags));
    for xi = 1:length(IND_fitting_lags)
        XI = IND_fitting_lags(xi);
        [NLLSF_COEF,JACOBIAN,SF_FIT,ConvergenceRecord] = ...
        nonlinear_lsqf(mean(SFedge_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2, DIST_VECedge(1:XI), ...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations);
        NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ [mean(SFedge_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2] - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
        SWOTedge_fitted_SpecSlopes(xi) = -NLLSF_COEF(2) - 1;
        SWOTedge_fitted_SpecSlopes_covariance(xi) = NLLSF_COEF_cov(2,2);
        disp(ConvergenceRecord(end))
    end

    % Cross-track
    xIND_fitting_lags = [4:2:10 12:4:20 25:5:50]';
    xSWOT_fitted_SpecSlopes            = nan(size(xIND_fitting_lags));
    xSWOT_fitted_SpecSlopes_covariance = nan(size(xIND_fitting_lags));
    for xi = 1:length(xIND_fitting_lags)
        XI = xIND_fitting_lags(xi);
        [NLLSF_COEF,JACOBIAN,SF_FIT,ConvergenceRecord] = ...
        nonlinear_lsqf(mean(xSF_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2, xDIST_VEC(1:XI), ...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations);
        NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ [mean(xSF_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2] - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
        xSWOT_fitted_SpecSlopes(xi) = -NLLSF_COEF(2) - 1;
        xSWOT_fitted_SpecSlopes_covariance(xi) = NLLSF_COEF_cov(2,2);
        disp(ConvergenceRecord(end))
    end
else
end; disp(' ')

plot(DIST_VEC(IND_fitting_lags), SWOT_fitted_SpecSlopes, 'k.-', 'MarkerSize',20,'LineWidth',1.5); hold on
errorbar(DIST_VEC(IND_fitting_lags), SWOT_fitted_SpecSlopes, SWOT_fitted_SpecSlopes_covariance, 'MarkerSize',20,'LineWidth',1.5,'Color','k','HandleVisibility','off')

plot(DIST_VEC(IND_fitting_lags), SWOTcent_fitted_SpecSlopes, '.-', 'MarkerSize',20,'LineWidth',1.5,'Color',CENTER_COLOR); hold on
errorbar(DIST_VEC(IND_fitting_lags), SWOTcent_fitted_SpecSlopes, SWOTcent_fitted_SpecSlopes_covariance, 'MarkerSize',20,'LineWidth',1.5,'Color',CENTER_COLOR,'HandleVisibility','off')

plot(DIST_VEC(IND_fitting_lags), SWOTedge_fitted_SpecSlopes, '.-', 'MarkerSize',20,'LineWidth',1.5,'Color',EDGE_COLOR); hold on
errorbar(DIST_VEC(IND_fitting_lags), SWOTedge_fitted_SpecSlopes, SWOTedge_fitted_SpecSlopes_covariance, 'MarkerSize',20,'LineWidth',1.5,'Color',EDGE_COLOR,'HandleVisibility','off')

plot(xDIST_VEC(xIND_fitting_lags), xSWOT_fitted_SpecSlopes, '.-', 'MarkerSize',20,'LineWidth',1.5,'Color',CROSS_COLOR); hold on
errorbar(xDIST_VEC(xIND_fitting_lags), xSWOT_fitted_SpecSlopes, xSWOT_fitted_SpecSlopes_covariance, 'MarkerSize',20,'LineWidth',1.5,'Color',CROSS_COLOR,'HandleVisibility','off')

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

LEG = legend(LEG_TEXT,'Location','southeast');
LEG.Position = LEG.Position + [0.0836 0.0331 0 0];
LEG.FontSize = 14;
LEG_TEXT = {};

% set(gcf,'Position',[-1759 392 1466 442])
set(gcf,'Position',[1 355 1041 442])


INPUT = input(['Do you want to save this figure? Enter any number for "yes" or\n push ' ...
    '"Enter" with a blank input for "No".\n']);

if isempty(INPUT)
    disp(['Image not saved'])
else
    % exportgraphics(gcf,...
    %     '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/REV1/F5.pdf',...
    %     'BackgroundColor','none','ContentType','vector')
    % disp(['Image saved'])
    disp('Save the figure after correcting for the structure function method shortcomings (next section).')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Log plot of S.F. of SWOT %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% with adjusted slopes %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Empirical correction (Table S1)
SF_slope_array = [... true; med_inf;   mean_inf;  StDev;
-4.4000   -4.2000   -4.0000   -3.9000   -3.8000   -3.7000   -3.6000   -3.5000   -3.4000   -3.3000   -3.2000   -3.1000   -3.0000   -2.9000   -2.8000   -2.7000   -2.6000   -2.5000   -2.4000   -2.3000   -2.2000   -2.1000   -2.0000   -1.9000   -1.8000   -1.7000   -1.6000   -1.5000   -1.4000   -1.3000   -1.2000   -1.1000   -1.0000;...
-2.9900   -2.9800   -2.9700   -2.9600   -2.9400   -2.9300   -2.9200   -2.8900   -2.8600   -2.8400   -2.7900   -2.7600   -2.7000   -2.6400   -2.5900   -2.5200   -2.4600   -2.3700   -2.3000   -2.2200   -2.1300   -2.0300   -1.9400   -1.8600   -1.7700   -1.6900   -1.5600   -1.5100   -1.4200   -1.3300   -1.2500   -1.1800   -1.1200;...
-2.9900   -2.9800   -2.9700   -2.9600   -2.9400   -2.9300   -2.9100   -2.8900   -2.8600   -2.8400   -2.7900   -2.7500   -2.7000   -2.6500   -2.6000   -2.5200   -2.4600   -2.3800   -2.2900   -2.2100   -2.1300   -2.0200   -1.9400   -1.8700   -1.7700   -1.6800   -1.5700   -1.5100   -1.4200   -1.3400   -1.2500   -1.1900   -1.1300;...
 0.0100    0.0100    0.0100    0.0100    0.0200    0.0200    0.0200    0.0300    0.0300    0.0300    0.0400    0.0500    0.0500    0.0600    0.0700    0.0700    0.0800    0.0700    0.0800    0.0800    0.0900    0.1000    0.1000    0.1200    0.0900    0.0800    0.0900    0.0900    0.0700    0.0700    0.0700    0.0500    0.0500;...
 ];
SF_slope_array = SF_slope_array';
TrueSlope = @(SLOPE,UNC) [interp1(SF_slope_array(:,3),SF_slope_array(:,1),SLOPE) + ...
                          1i*sqrt([UNC].^2 + ...
                                  [SF_slope_array(dsearchn(SF_slope_array(:,3),SLOPE),4)].^2)...
                         ];
% % % % %


close all
figure('Color',[1 1 1])
dBIN = 0.1;
LEG_TEXT = {};
CO = colororder;
m2_to_cm2 = 100*100; % =10^4 if you are to plot cm^2; =1 if m^2
XSCALE = 'lin'; YSCALE = 'lin';
CENTER_COLOR = [77  191  239]/255; % [ 41, 203, 255 ]./255; % Custom color
EDGE_COLOR   = [218   83   25]/255; % [ 244, 179, 0  ]./255; % Custom color
CROSS_COLOR  = [185, 19, 209 ]/255; % [ 244, 179, 0  ]./255; % Custom color

tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');

nexttile;%subplot(121) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % SWOT_ALL
% LEG_TEXT = {LEG_TEXT{:},'SWOT cal/val all data (29.9N-50.2N)'};
LEG_TEXT = {LEG_TEXT{:},['Along-track SWOT cal/val ' newline 'all data (' num2str(round(min([SWOT.lat{end}(:) ; SWOT.lat{end-1}(:)]),1)) 'N-' ...
                                                   num2str(round(max([SWOT.lat{end}(:) ; SWOT.lat{end-1}(:)]),1)) 'N)']};
fill([flip(DIST_VEC); DIST_VEC], ...
     [flip(mean(SF_MAT_swot,2,'omitnan') + std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))); ...
           mean(SF_MAT_swot,2,'omitnan') - std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))], ...
     [1 1 1]*0.7, 'EdgeAlpha',0,'HandleVisibility','off'); hold on
Mean_SF_MAT_swot = mean(SF_MAT_swot,2,'omitnan')*m2_to_cm2;
STD_SF_MAT_swot = std(SF_MAT_swot,[],2,'omitnan')*m2_to_cm2./sqrt(sum(isfinite(SF_MAT_swot),2));
% plot(DIST_VEC, Mean_SF_MAT_swot, 'k.-', 'MarkerSize',10); hold on
errorbar(DIST_VEC, Mean_SF_MAT_swot,STD_SF_MAT_swot, 'k.-', 'MarkerSize',15); hold on


% % % % SWOT 12 most central tracks (6 per swath)
LEG_TEXT = {LEG_TEXT{:},['Along-track SWOT cal/val ' newline '(12 central tracks)']};
fill([flip(DIST_VECcent); DIST_VECcent], ...
     [flip(mean(SFcent_MAT_swot,2,'omitnan') + std(SFcent_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SFcent_MAT_swot(:,:)), 2))); ...
           mean(SFcent_MAT_swot,2,'omitnan') - std(SFcent_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SFcent_MAT_swot(:,:)), 2))], ...
     [1 1 1]*0.7, 'EdgeAlpha',0,'HandleVisibility','off'); hold on
Mean_SFcent_MAT_swot = mean(SFcent_MAT_swot,2,'omitnan')*m2_to_cm2;
STD_SFcent_MAT_swot = std(SFcent_MAT_swot,[],2,'omitnan')*m2_to_cm2./sqrt(sum(isfinite(SFcent_MAT_swot),2));
% plot(DIST_VECcent, Mean_SFcent_MAT_swot, '.-', 'Color',BLUE); hold on
errorbar(DIST_VECcent, Mean_SFcent_MAT_swot, STD_SFcent_MAT_swot, '.-', 'Color',CENTER_COLOR); hold on

% % % % SWOT 12 edge-most tracks (6 per swath)
LEG_TEXT = {LEG_TEXT{:},['Along-track SWOT cal/val ' newline '(12 edge-most tracks)']};
fill([flip(DIST_VECedge); DIST_VECedge], ...
     [flip(mean(SFedge_MAT_swot,2,'omitnan') + std(SFedge_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SFedge_MAT_swot(:,:)), 2))); ...
           mean(SFedge_MAT_swot,2,'omitnan') - std(SFedge_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SFedge_MAT_swot(:,:)), 2))], ...
     [1 1 1]*0.7, 'EdgeAlpha',0,'HandleVisibility','off'); hold on
Mean_SFedge_MAT_swot = mean(SFedge_MAT_swot,2,'omitnan')*m2_to_cm2;
STD_SFedge_MAT_swot = std(SFedge_MAT_swot,[],2,'omitnan')*m2_to_cm2./sqrt(sum(isfinite(SFedge_MAT_swot),2));
% plot(DIST_VECedge, Mean_SFedge_MAT_swot, '.-', 'Color',ORANGE); hold on
errorbar(DIST_VECedge, Mean_SFedge_MAT_swot, STD_SFedge_MAT_swot, '.-', 'Color',EDGE_COLOR); hold on


% Plot cross-track structure function, but with the last 8 elements omitted
% because the last 7 are NaN and the one before is clearly spurious:
LEG_TEXT = {LEG_TEXT{:},['Cross-track']};
Mean_xSF_MAT_swot = mean(xSF_MAT_swot,2,'omitnan')*m2_to_cm2;
STD_xSF_MAT_swot = std(xSF_MAT_swot,[],2,'omitnan')*m2_to_cm2./sqrt(sum(isfinite(xSF_MAT_swot),2));
% loglog(xDIST_VEC(1:[end-8]), Mean_xSF_MAT_swot(1:[end-8],:), 'r.-'); hold on
errorbar(xDIST_VEC(1:[end-8]), Mean_xSF_MAT_swot(1:[end-8],:), STD_xSF_MAT_swot(1:[end-8],:), '.-','Color',CROSS_COLOR); hold on


% % % SWOT_at_S_P_G
LEG_TEXT = {LEG_TEXT{:},['SWOT cal/val at ' newline 'Moorings+Gliders']};
BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_SWOT_moor_12hr,Separation_SWOT_moor_12hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
% plot(Separation_SWOT_moor_12hr, StructFunc_SWOT_moor_12hr, '.', 'MarkerSize',10); hold on
errorbar_stdofmean(BIN_CENTR, Separation_SWOT_moor_12hr, m2_to_cm2*StructFunc_SWOT_moor_12hr, BIN_EDGES, [1 1 1]*0.7); hold on
Mean_SWOT_at_S_P_G = interval_avg(Separation_SWOT_moor_12hr, m2_to_cm2*StructFunc_SWOT_moor_12hr, BIN_EDGES);
% plot(BIN_CENTR,Mean_SWOT_at_S_P_G,'r') % verify that this is the quantity of interest

% %% Plot full CROSS-TRACK structure functions
% 
% close all
% 
% figure('Color',[1 1 1])
% for jj = 1:round(0.1*size(xSF_MAT_swot,2)) % 1:size(SF_MAT_swot,2)
%     % plot(DIST_VEC,SF_MAT_swot(:,jj)*100^2,'-','Color',[1 1 1]*0.7); hold on
%     loglog(xDIST_VEC,xSF_MAT_swot(:,jj)*100^2,'-','Color',[1 1 1]*0.7); hold on
% end
% 
% fill([flip(xDIST_VEC); xDIST_VEC], ...
%      [flip(mean(xSF_MAT_swot,2,'omitnan') + std(xSF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(xSF_MAT_swot(:,:)), 2))); ...
%            mean(xSF_MAT_swot,2,'omitnan') - std(xSF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(xSF_MAT_swot(:,:)), 2))]*100^2, ...
%      [1 1 1]*0.7); hold on
% % plot(DIST_VEC, mean(SF_MAT_swot,2,'omitnan')*100^2, 'k.-'); hold on
% loglog(xDIST_VEC, mean(xSF_MAT_swot,2,'omitnan')*100^2, 'k.-'); hold on
% xlabel('Distance (km)')
% ylabel('SF (cm^2)')
% set(gca,'FontSize',24)


xlabel('Separation (km)')
ylabel('$\overline{\big(\mathrm{SSHA}(x) - \mathrm{SSHA}(x+s)\big)^2}$ \quad (cm$^2$)','Interpreter','latex','FontSize',16)

set(gca,'XLim',DIST_VEC([1 end]))
set(gca,'YLim',[0.1 100]); set(gca,'XLim',[0 700])
set(gca,'yscale',YSCALE,'xscale',XSCALE)
set(gca,'yscale','log','xscale','log')
set(gca,'FontSize',20); grid on



% Fitting parameters
FitBasisFunction = {'a*t.^b','a','b'};
FitParameters0 = [0.1,1.5];
FitBasisDerivatives = {'t.^b','a*log(t).*t.^b'};
ToleranceLevel = 0.0001; % NLLSF parameter, tells change in parameter magnitude, if small then the result was stable
N_nllsf_iterations = 20; % number of iterations before we give up


nexttile;%subplot(122) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear SWOT_fitted_SpecSlopes
if ~exist('SWOT_fitted_SpecSlopes','var')
    % i.e. only run this if we haven't yet calculated the slope as a
    % function of maximum fitting interval.
    
    % Along-track
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

    % Along-track, 12 central tracks
    % IND_fitting_lags = [4:2:10 12:4:20 25:5:50 60:10:100 150]';
    SWOTcent_fitted_SpecSlopes            = nan(size(IND_fitting_lags));
    SWOTcent_fitted_SpecSlopes_covariance = nan(size(IND_fitting_lags));
    for xi = 1:length(IND_fitting_lags)
        XI = IND_fitting_lags(xi);
        [NLLSF_COEF,JACOBIAN,SF_FIT,ConvergenceRecord] = ...
        nonlinear_lsqf(mean(SFcent_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2, DIST_VECcent(1:XI), ...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations);
        NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ [mean(SFcent_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2] - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
        SWOTcent_fitted_SpecSlopes(xi) = -NLLSF_COEF(2) - 1;
        SWOTcent_fitted_SpecSlopes_covariance(xi) = NLLSF_COEF_cov(2,2);
        disp(ConvergenceRecord(end))
    end

    % Along-track, 12 edge-most tracks
    % IND_fitting_lags = [4:2:10 12:4:20 25:5:50 60:10:100 150]';
    SWOTedge_fitted_SpecSlopes            = nan(size(IND_fitting_lags));
    SWOTedge_fitted_SpecSlopes_covariance = nan(size(IND_fitting_lags));
    for xi = 1:length(IND_fitting_lags)
        XI = IND_fitting_lags(xi);
        [NLLSF_COEF,JACOBIAN,SF_FIT,ConvergenceRecord] = ...
        nonlinear_lsqf(mean(SFedge_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2, DIST_VECedge(1:XI), ...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations);
        NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ [mean(SFedge_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2] - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
        SWOTedge_fitted_SpecSlopes(xi) = -NLLSF_COEF(2) - 1;
        SWOTedge_fitted_SpecSlopes_covariance(xi) = NLLSF_COEF_cov(2,2);
        disp(ConvergenceRecord(end))
    end

    % Cross-track
    xIND_fitting_lags = [4:2:10 12:4:20 25:5:50]';
    xSWOT_fitted_SpecSlopes            = nan(size(xIND_fitting_lags));
    xSWOT_fitted_SpecSlopes_covariance = nan(size(xIND_fitting_lags));
    for xi = 1:length(xIND_fitting_lags)
        XI = xIND_fitting_lags(xi);
        [NLLSF_COEF,JACOBIAN,SF_FIT,ConvergenceRecord] = ...
        nonlinear_lsqf(mean(xSF_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2, xDIST_VEC(1:XI), ...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations);
        NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ [mean(xSF_MAT_swot(1:XI,:),2,'omitnan')*m2_to_cm2] - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
        xSWOT_fitted_SpecSlopes(xi) = -NLLSF_COEF(2) - 1;
        xSWOT_fitted_SpecSlopes_covariance(xi) = NLLSF_COEF_cov(2,2);
        disp(ConvergenceRecord(end))
    end

    % Adjust with the empirical correction
    SWOT_fitted_SpecSlopes                = real(TrueSlope(SWOT_fitted_SpecSlopes,SWOT_fitted_SpecSlopes_covariance));
    SWOTcent_fitted_SpecSlopes            = real(TrueSlope(SWOTcent_fitted_SpecSlopes,SWOTcent_fitted_SpecSlopes_covariance));
    SWOTedge_fitted_SpecSlopes            = real(TrueSlope(SWOTedge_fitted_SpecSlopes,SWOTedge_fitted_SpecSlopes_covariance));
    xSWOT_fitted_SpecSlopes               = real(TrueSlope(xSWOT_fitted_SpecSlopes,xSWOT_fitted_SpecSlopes_covariance));
    SWOT_fitted_SpecSlopes_covariance     = imag(TrueSlope(SWOT_fitted_SpecSlopes,SWOT_fitted_SpecSlopes_covariance));
    SWOTcent_fitted_SpecSlopes_covariance = imag(TrueSlope(SWOTcent_fitted_SpecSlopes,SWOTcent_fitted_SpecSlopes_covariance));
    SWOTedge_fitted_SpecSlopes_covariance = imag(TrueSlope(SWOTedge_fitted_SpecSlopes,SWOTedge_fitted_SpecSlopes_covariance));
    xSWOT_fitted_SpecSlopes_covariance    = imag(TrueSlope(xSWOT_fitted_SpecSlopes,xSWOT_fitted_SpecSlopes_covariance));

else
end; disp(' ')

plot(DIST_VEC(IND_fitting_lags), SWOT_fitted_SpecSlopes, 'k.-', 'MarkerSize',20,'LineWidth',1.5); hold on
errorbar(DIST_VEC(IND_fitting_lags), SWOT_fitted_SpecSlopes, SWOT_fitted_SpecSlopes_covariance, 'MarkerSize',20,'LineWidth',1.5,'Color','k','HandleVisibility','off')

plot(DIST_VEC(IND_fitting_lags), SWOTcent_fitted_SpecSlopes, '.-', 'MarkerSize',20,'LineWidth',1.5,'Color',CENTER_COLOR); hold on
errorbar(DIST_VEC(IND_fitting_lags), SWOTcent_fitted_SpecSlopes, SWOTcent_fitted_SpecSlopes_covariance, 'MarkerSize',20,'LineWidth',1.5,'Color',CENTER_COLOR,'HandleVisibility','off')

plot(DIST_VEC(IND_fitting_lags), SWOTedge_fitted_SpecSlopes, '.-', 'MarkerSize',20,'LineWidth',1.5,'Color',EDGE_COLOR); hold on
errorbar(DIST_VEC(IND_fitting_lags), SWOTedge_fitted_SpecSlopes, SWOTedge_fitted_SpecSlopes_covariance, 'MarkerSize',20,'LineWidth',1.5,'Color',EDGE_COLOR,'HandleVisibility','off')

plot(xDIST_VEC(xIND_fitting_lags), xSWOT_fitted_SpecSlopes, '.-', 'MarkerSize',20,'LineWidth',1.5,'Color',CROSS_COLOR); hold on
errorbar(xDIST_VEC(xIND_fitting_lags), xSWOT_fitted_SpecSlopes, xSWOT_fitted_SpecSlopes_covariance, 'MarkerSize',20,'LineWidth',1.5,'Color',CROSS_COLOR,'HandleVisibility','off')

% Do the same (get slope as a function of max. evaluated separation) for
% SWOT at instrument locations:
SWOT_CalVal_fitted_SpecSlopes            = nan(length(BIN_CENTR),1);
SWOT_CalVal_fitted_SpecSlopes_maxlag     = nan(length(BIN_CENTR),1);
SWOT_CalVal_fitted_SpecSlopes_covariance = nan(length(BIN_CENTR),1);
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
SWOT_CalVal_fitted_SpecSlopes            = real(TrueSlope(SWOT_CalVal_fitted_SpecSlopes,SWOT_CalVal_fitted_SpecSlopes_covariance));
SWOT_CalVal_fitted_SpecSlopes_covariance = imag(TrueSlope(SWOT_CalVal_fitted_SpecSlopes,SWOT_CalVal_fitted_SpecSlopes_covariance));
plot(SWOT_CalVal_fitted_SpecSlopes_maxlag, SWOT_CalVal_fitted_SpecSlopes, ...
    '.-', 'Color', [1 1 1]*0.7, 'MarkerSize',20,'LineWidth',1.5)
errorbar(SWOT_CalVal_fitted_SpecSlopes_maxlag, SWOT_CalVal_fitted_SpecSlopes, SWOT_CalVal_fitted_SpecSlopes_covariance, ...
    'Color', [1 1 1]*0.7, 'MarkerSize',20,'LineWidth',1.5)

xlabel({'Max. separation (km) of';'fitted structure function'})
ylabel({'Inferred Spectral Slope \lambda';'(corrected)'})
set(gca,'FontSize',20); grid on

LEG = legend(LEG_TEXT,'Location','southeast');
LEG.Position = LEG.Position + [0.0836 0.0331 0 0];
LEG.FontSize = 14;
LEG_TEXT = {};

% set(gcf,'Position',[-1759 392 1466 442])
set(gcf,'Position',[1 355 1041 461])


INPUT = input(['Do you want to save this figure? Enter any number for "yes" or\n push ' ...
    '"Enter" with a blank input for "No".\n']);

if isempty(INPUT)
    disp(['Image not saved'])
else
    % exportgraphics(gcf,...
    %     '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/F4.pdf',...
    %     'BackgroundColor','none','ContentType','vector')
    exportgraphics(gcf,...
        '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/REV1/F5.pdf',...
        'BackgroundColor','none','ContentType','vector')
    disp(['Image saved'])
end



%% Auxiliary function

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


