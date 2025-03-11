%% Supplemental Material - Data Analysis
%% Examine and compare cross-track and along-track structure functions from SWOT

% SWOT = load('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/SWOTdata_CCS.mat',...
%             'time','ssha_karin_2','ssha_karin_2_qual','name','lat','lon','internal_tide_hret','height_cor_xover','height_cor_xover_qual');
SWOT = load('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/SWOTdata_CCS_33_385.mat',... More limited latitude range than above
            'time','ssha_karin_2','ssha_karin_2_qual','name','lat','lon', ...
            'internal_tide_hret','height_cor_xover','height_cor_xover_qual','cross_track_distance');
warning('Be sure to download the updated version with SWH when it is done processing on the server.')
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
XTRACK = 150;
SLOPE = -5/3; SLOPE = -2;
SEGS = 3;

IND = 1;
ALONGTRACK_DISTANCE = abs([SWOT.lon{IND}(TRACK,:) + 1i*SWOT.lat{IND}(TRACK,:)] - [SWOT.lon{IND}(TRACK,1) + 1i*SWOT.lat{IND}(TRACK,1)]);
CROSSTRACK_DISTANCE = abs([SWOT.lon{IND}(:,XTRACK) + 1i*SWOT.lat{IND}(:,XTRACK)] - [SWOT.lon{IND}(1,XTRACK) + 1i*SWOT.lat{IND}(1,XTRACK)]);

dx = 111*median(diff(ALONGTRACK_DISTANCE));
xdx = 111*median(diff(CROSSTRACK_DISTANCE));

%% Calculate mean ascending and descending tracks in order to remove the geoid per JW's advice

warning(['Also plot the spectra resulting from the removed and unremoved ' ...
         'alongside each other for the supplementary data'])

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

%% ALONG-TRACK Structure functions from full tracks, every cell

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

%% Plot full ALONG-TRACK structure functions

close all

figure('Color',[1 1 1])
for jj = 1:round(0.1*size(SF_MAT_swot,2)) % 1:size(SF_MAT_swot,2)
    % plot(DIST_VEC,SF_MAT_swot(:,jj)*100^2,'-','Color',[1 1 1]*0.7); hold on
    loglog(DIST_VEC,SF_MAT_swot(:,jj)*100^2,'-','Color',[1 1 1]*0.7); hold on
end

fill([flip(DIST_VEC); DIST_VEC], ...
     [flip(mean(SF_MAT_swot,2,'omitnan') + std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))); ...
           mean(SF_MAT_swot,2,'omitnan') - std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))]*100^2, ...
     [1 1 1]*0.7); hold on
% plot(DIST_VEC, mean(SF_MAT_swot,2,'omitnan')*100^2, 'k.-'); hold on
loglog(DIST_VEC, mean(SF_MAT_swot,2,'omitnan')*100^2, 'k.-'); hold on
xlabel('Distance (km)')
ylabel('SF (cm^2)')
set(gca,'FontSize',24)

%% CROSS-TRACK Structure functions from full tracks, every cell

close all

% error(['This calculation takes a while'])

xSF_ARRAY_swot = nan(size(SWOT.ssha_karin_2{1},1) - 1, size(SWOT.ssha_karin_2{1},2), length(SWOT.ssha_karin_2));
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
%%
xSF_MAT_swot = nan(size(xSF_ARRAY_swot,1),size(xSF_ARRAY_swot,2)*size(xSF_ARRAY_swot,3));
for kk = 1:size(xSF_ARRAY_swot,3)
    xSF_MAT_swot(:,[1:size(xSF_ARRAY_swot,2)] + size(xSF_ARRAY_swot,2)*[kk-1]) = xSF_ARRAY_swot(:,:,kk);
end

%% Plot full CROSS-TRACK structure functions

close all

figure('Color',[1 1 1])
for jj = 1:round(0.1*size(xSF_MAT_swot,2)) % 1:size(SF_MAT_swot,2)
    % plot(DIST_VEC,SF_MAT_swot(:,jj)*100^2,'-','Color',[1 1 1]*0.7); hold on
    loglog(xDIST_VEC,xSF_MAT_swot(:,jj)*100^2,'-','Color',[1 1 1]*0.7); hold on
end

fill([flip(xDIST_VEC); xDIST_VEC], ...
     [flip(mean(xSF_MAT_swot,2,'omitnan') + std(xSF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(xSF_MAT_swot(:,:)), 2))); ...
           mean(xSF_MAT_swot,2,'omitnan') - std(xSF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(xSF_MAT_swot(:,:)), 2))]*100^2, ...
     [1 1 1]*0.7); hold on
% plot(DIST_VEC, mean(SF_MAT_swot,2,'omitnan')*100^2, 'k.-'); hold on
loglog(xDIST_VEC, mean(xSF_MAT_swot,2,'omitnan')*100^2, 'k.-'); hold on
xlabel('Distance (km)')
ylabel('SF (cm^2)')
set(gca,'FontSize',24)

%% Compare ALONG-TRACK and CROSS-TRACK SF's side by side

close all

figure('Color',[1 1 1])


fill([flip(DIST_VEC); DIST_VEC], ... STD of mean
     [flip(mean(SF_MAT_swot,2,'omitnan') + std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))); ...
           mean(SF_MAT_swot,2,'omitnan') - std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))]*100^2, ...
     [1 1 1]*0.7,'HandleVisibility','off'); hold on
loglog(DIST_VEC, mean(SF_MAT_swot,2,'omitnan')*100^2, 'k.-'); hold on

fill([flip(xDIST_VEC); xDIST_VEC], ... STD of mean
     [flip(mean(xSF_MAT_swot,2,'omitnan') + std(xSF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(xSF_MAT_swot(:,:)), 2))); ...
           mean(xSF_MAT_swot,2,'omitnan') - std(xSF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(xSF_MAT_swot(:,:)), 2))]*100^2, ...
     [1 0 0]*0.7,'HandleVisibility','off','EdgeColor','r'); hold on
loglog(xDIST_VEC, mean(xSF_MAT_swot,2,'omitnan')*100^2, 'r.-'); hold on
xlabel('Distance (km)')
ylabel('SF (cm^2)')
set(gca,'FontSize',24)
set(gca,'XScale','log','YScale','log')
legend('Along-track','Cross-track')
% %%
% figure
% fill([flip(DIST_VEC); DIST_VEC], ... STD of mean
%      [flip(mean(SF_MAT_swot,2,'omitnan') + std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))); ...
%            mean(SF_MAT_swot,2,'omitnan') - std(SF_MAT_swot,0,2,'omitnan')./sqrt(sum(isfinite(SF_MAT_swot(:,:)), 2))]*100^2, ...
%      [0 1 1]*0.7,'HandleVisibility','off'); hold on
%% Slopes - linear least squares fit in loglog space

% linear least squares fit in loglog space - ALONG TRACK
close all
figure; disp('-----')
loglog(DIST_VEC, mean(SF_MAT_swot,2,'omitnan'),'k.-','LineWidth',3); hold on
YY = mean(SF_MAT_swot,2,'omitnan');
XX = [log10(DIST_VEC(isfinite(YY))), ones(sum(isfinite(YY)),1)]\log10(YY(isfinite(YY))); disp(-1-XX(1))
loglog(DIST_VEC, 10.^[ [log10(DIST_VEC), ones(size(DIST_VEC))]*XX ], '.-')
XX = [log10(DIST_VEC(isfinite(YY) & DIST_VEC<100)), ones(sum(isfinite(YY) & DIST_VEC<100),1)]\...
      log10(YY(isfinite(YY) & DIST_VEC<100)); disp(-1-XX(1))
loglog(DIST_VEC, 10.^[ [log10(DIST_VEC), ones(size(DIST_VEC))]*XX ], '.-')
XX = [log10(DIST_VEC(isfinite(YY) & DIST_VEC<40)), ones(sum(isfinite(YY) & DIST_VEC<40),1)]\...
      log10(YY(isfinite(YY) & DIST_VEC<40)); disp(-1-XX(1))
loglog(DIST_VEC, 10.^[ [log10(DIST_VEC), ones(size(DIST_VEC))]*XX ], '.-')

% linear least squares fit in loglog space - CROSS TRACK
figure; disp(' ')
loglog(xDIST_VEC, mean(xSF_MAT_swot,2,'omitnan'),'k.-','LineWidth',3); hold on
YY = mean(xSF_MAT_swot,2,'omitnan');
XX = [log10(xDIST_VEC(isfinite(YY))), ones(sum(isfinite(YY)),1)]\log10(YY(isfinite(YY))); disp(-1-XX(1))
loglog(xDIST_VEC, 10.^[ [log10(xDIST_VEC), ones(size(xDIST_VEC))]*XX ], '.-')
XX = [log10(xDIST_VEC(isfinite(YY) & xDIST_VEC<100)), ones(sum(isfinite(YY) & xDIST_VEC<100),1)]\...
      log10(YY(isfinite(YY) & xDIST_VEC<100)); disp(-1-XX(1))
loglog(xDIST_VEC, 10.^[ [log10(xDIST_VEC), ones(size(xDIST_VEC))]*XX ], '.-')
XX = [log10(xDIST_VEC(isfinite(YY) & xDIST_VEC<40)), ones(sum(isfinite(YY) & xDIST_VEC<40),1)]\...
      log10(YY(isfinite(YY) & xDIST_VEC<40)); disp(-1-XX(1))
loglog(xDIST_VEC, 10.^[ [log10(xDIST_VEC), ones(size(xDIST_VEC))]*XX ], '.-')

% set(gca,'XScale','lin')
% set(gca,'YScale','lin')


%% Slopes - non-linear least squares fit in regular space

FitBasisFunction = {'a*t.^b','a','b'};
FitParameters0 = [0.1,1];
FitBasisDerivatives = {'t.^b','a*log(t).*t.^b'};
ToleranceLevel = 0.00001; % NLLSF parameter, tells change in parameter magnitude, if small then the result was stable
N_nllsf_iterations = 10; % number of iterations before we give up

% non-linear least squares fit in regular space - ALONG TRACK
close all
figure; disp('-----')
plot(DIST_VEC, mean(SF_MAT_swot,2,'omitnan'),'k.-','LineWidth',3); hold on
[XX,~,YY_fit] = ...
    nonlinear_lsqf(mean(SF_MAT_swot,2,'omitnan'),DIST_VEC, ...
                   FitBasisFunction,FitParameters0,FitBasisDerivatives, ...
                   ToleranceLevel,N_nllsf_iterations);
plot(DIST_VEC(1:length(YY_fit)),YY_fit,'o:','MarkerSize',12,'LineWidth',2); disp(-1-XX(2))
[XX,~,YY_fit] = ...
    nonlinear_lsqf(mean(SF_MAT_swot(DIST_VEC<100,:),2,'omitnan'),DIST_VEC(DIST_VEC<100), ...
                   FitBasisFunction,FitParameters0,FitBasisDerivatives, ...
                   ToleranceLevel,N_nllsf_iterations);
plot(DIST_VEC(1:length(YY_fit)),YY_fit,'o:','MarkerSize',12,'LineWidth',2); disp(-1-XX(2))
[XX,~,YY_fit] = ...
    nonlinear_lsqf(mean(SF_MAT_swot(DIST_VEC<40,:),2,'omitnan'),DIST_VEC(DIST_VEC<40), ...
                   FitBasisFunction,FitParameters0,FitBasisDerivatives, ...
                   ToleranceLevel,N_nllsf_iterations);
plot(DIST_VEC(1:length(YY_fit)),YY_fit,'o:','MarkerSize',12,'LineWidth',2); disp(-1-XX(2))
% set(gca,'XScale','log')
% set(gca,'YScale','log')

% non-linear least squares fit in regular space - ALONG TRACK
% close all
figure; disp('-----')
plot(xDIST_VEC, mean(xSF_MAT_swot,2,'omitnan'),'k.-','LineWidth',3); hold on
[XX,~,YY_fit] = ...
    nonlinear_lsqf(mean(xSF_MAT_swot,2,'omitnan'),xDIST_VEC, ...
                   FitBasisFunction,FitParameters0,FitBasisDerivatives, ...
                   ToleranceLevel,N_nllsf_iterations);
plot(xDIST_VEC(1:length(YY_fit)),YY_fit,'o:','MarkerSize',12,'LineWidth',2); disp(-1-XX(2))
[XX,~,YY_fit] = ...
    nonlinear_lsqf(mean(xSF_MAT_swot(xDIST_VEC<100,:),2,'omitnan'),xDIST_VEC(xDIST_VEC<100), ...
                   FitBasisFunction,FitParameters0,FitBasisDerivatives, ...
                   ToleranceLevel,N_nllsf_iterations);
plot(xDIST_VEC(1:length(YY_fit)),YY_fit,'o:','MarkerSize',12,'LineWidth',2); disp(-1-XX(2))
[XX,~,YY_fit] = ...
    nonlinear_lsqf(mean(xSF_MAT_swot(xDIST_VEC<40,:),2,'omitnan'),xDIST_VEC(xDIST_VEC<40), ...
                   FitBasisFunction,FitParameters0,FitBasisDerivatives, ...
                   ToleranceLevel,N_nllsf_iterations);
plot(xDIST_VEC(1:length(YY_fit)),YY_fit,'o:','MarkerSize',12,'LineWidth',2); disp(-1-XX(2))
% set(gca,'XScale','log')
% set(gca,'YScale','log')