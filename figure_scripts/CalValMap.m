%% FIGURE 1 - CalValMap

% Load SWOT data:
SWOT = load('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/SWOTdata_CCS.mat');
% 871.4267 MB

%% POSTER
%% Large panel - background map

close all

figure
M_P = m_proj('Lambert',...
    'longitudes',[-130 -116],...
    'latitudes',[30 45]);
set(gcf,'color','w')
% COAST = m_gshhs_l('patch',0.8*[1 1 1]);
COAST = m_gshhs_h('patch',[93 165 38]/255);

hold on

% % Large panel - crossover

% Data is at:
% kachelein@eddy:/mnt/flow/swot/KaRIn/SWOT_L2_LR_SSH_2.0_AllToDate
% Files are of this form:
% SWOT_L2_LR_SSH_Basic_015_382_20240521T224612_20240521T233658_PIC0_01.nc

% figure; imagesc(SWOT.ssha_karin_2{150} - SWOT.height_cor_xover{150}); colorbar

one_if_false_nan_if_true = @(IN) (~IN)./(~IN);

II = 173; % max = 180
for ii = [II + [0 1]]
    m_pcolor_centered(SWOT.lon{ii}, SWOT.lat{ii}, ...
        [SWOT.ssha_karin_2{ii} + SWOT.height_cor_xover{ii}] .* ...
        one_if_false_nan_if_true(SWOT.ssha_karin_2_qual{ii} + SWOT.height_cor_xover_qual{ii}) );
    % disp(['Time = ' datestr(mean(SWOT.time{ii})) ' (ii = ' num2str(ii) ')'])
    disp(['Time = ' datestr(SWOT.time{ii}(dsearchn(SWOT.lat{ii}(35,:)', 36)) ) ' (ii = ' num2str(ii) ')'])
    % ^ 36 = approx. xover lat., 35 = index of center of the swath
end
colormap(bwr(64)*0.99)
CB = colorbar;
CB.Label.String = 'SSHA (m)';
CB.Position = [0.5244    0.1620    0.0193    0.1935];
CB.FontSize = 12;
clim([-1 1]*0.2)
CB.Ticks = [-.2 -.1 0 .1 .2];

% % Box indicating in-set:
INSET_LAT = [35 37];
INSET_LON = [-126 -124.5];
collimate = @(IN) IN(:);
m_plot([INSET_LON flip(INSET_LON) , INSET_LON(1)]', ...
       [collimate([INSET_LAT ; INSET_LAT]) ; INSET_LAT(1)], 'k', 'LineWidth', 1)

% m_grid('box','fancy', 'backgroundcolor',[1 1 1]*0.9);
m_grid('box','fancy', 'backgroundcolor',[128 226 255]/255);

set(gcf,'Position',[262    69   962   728])

% exportgraphics(gcf,...
%     '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/component_figures/F1_USWC_map.pdf',...
%     'BackgroundColor','none','ContentType','vector')

%% Small panel - background map

close all

figure
M_P = m_proj('Lambert',...
    'longitudes',INSET_LON,...
    'latitudes',INSET_LAT);
set(gcf,'color','w')
% COAST = m_gshhs_l('patch',0.8*[1 1 1]);
COAST = m_gshhs_l('patch',[93 165 38]/255);

hold on

% % Small panel - crossover
II = 173; % max = 180
for ii = [II + [0 1]]
    m_pcolor_centered(SWOT.lon{ii}, SWOT.lat{ii}, ...
        [SWOT.ssha_karin_2{ii} + SWOT.height_cor_xover{ii}] .* ...
        one_if_false_nan_if_true(SWOT.ssha_karin_2_qual{ii} + SWOT.height_cor_xover_qual{ii}) );
    % disp(['Time = ' datestr(mean(SWOT.time{ii})) ' (ii = ' num2str(ii) ')'])
    disp(['Time = ' datestr(SWOT.time{ii}(dsearchn(SWOT.lat{ii}(35,:)', 36)) ) ' (ii = ' num2str(ii) ')'])
    % ^ 36 = approx. xover lat., 35 = index of center of the swath
end
colormap(bwr(64)*0.999)
clim([-1 1]*0.2)
% % % % % 

m_grid('box','fancy', 'backgroundcolor',[128 226 255]/255, ...
       'xtick', [round(INSET_LON(1)):round(INSET_LON(end))], ...
       'ytick', [round(INSET_LAT(1)):round(INSET_LAT(end))])

% % Small panel - gliders

% % gliders
if ~exist('Glider_loc','var')
    G_DIR = dir('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/GLIDERS/ru*_L2_Processed_2*.nc');
    Glider_loc = {};
    for ii = 1:length(G_DIR)
        Glider_loc{ii} = [ncread([G_DIR(ii).folder '/' G_DIR(ii).name],'LONGITUDE_PROFILE') + 1i*ncread([G_DIR(ii).folder '/' G_DIR(ii).name],'LATITUDE_PROFILE')];

        % % % Old QC correction, no longer needed:
        % GliderQC = ncread([G_DIR(ii).folder '/' G_DIR(ii).name],'QC_FLAG');
        % Glider_loc{ii} = Glider_loc{ii}.*[~GliderQC./~GliderQC];
    end
else
end

GLIDER_COLOR = {[1 1 1]*0.6 , [0 1 1]*0.6};
for ii = 1:length(Glider_loc)
    m_plot(real(Glider_loc{ii}), imag(Glider_loc{ii}), '.-','Color',GLIDER_COLOR{ii})
end

% % Small panel - moorings
if ~exist('Mooring_nom_loc','var')
    M_P_DIR = dir('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/PROF/*.nc');
    Mooring_nom_loc = nan(length(length(M_P_DIR)),2);
    for ii = 1:length(M_P_DIR)
        Mooring_nom_loc(ii,1) = ncread([M_P_DIR(ii).folder '/' M_P_DIR(ii).name],'LONGITUDE');
        Mooring_nom_loc(ii,2) = ncread([M_P_DIR(ii).folder '/' M_P_DIR(ii).name],'LATITUDE');
    end
else
end
for ii = 1:length(M_P_DIR)
    SP_sign = sign(strcmp(M_P_DIR(ii).name(33),'S') - 0.5);
    if SP_sign > 0 % S
        m_plot(Mooring_nom_loc(ii,1), Mooring_nom_loc(ii,2), '.k','MarkerSize',20)
        m_text(Mooring_nom_loc(ii,1) + 0.1, Mooring_nom_loc(ii,2) + 0.02,...
           M_P_DIR(ii).name(33:34), 'FontSize',16)
    else % P
        m_plot(Mooring_nom_loc(ii,1), Mooring_nom_loc(ii,2), '*k','MarkerSize',10,'LineWidth',2)
        m_text(Mooring_nom_loc(ii,1) - 0.2, Mooring_nom_loc(ii,2) - 0.02,...
           M_P_DIR(ii).name(33:34), 'FontSize',16)
    end
end
% m_plot(Mooring_nom_loc(:,1), Mooring_nom_loc(:,2), '.k','MarkerSize',20)

set(gcf,'Position',[255   207   693   520])

% exportgraphics(gcf,...
%     '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/component_figures/F1_xover_inset.pdf',...
%     'BackgroundColor','none','ContentType','vector')

%% PAPER
%%      Large panel - background map

close all

figure
M_P = m_proj('Lambert',...
    'longitudes',[-130 -116],...
    'latitudes',[30 45]);
set(gcf,'color','w')
% COAST = m_gshhs_l('patch',0.8*[1 1 1]);
COAST = m_gshhs_h('patch',[93 165 38]/255);

hold on

% % Large panel - crossover

% Data is at:
% kachelein@eddy:/mnt/flow/swot/KaRIn/SWOT_L2_LR_SSH_2.0_AllToDate
% Files are of this form:
% SWOT_L2_LR_SSH_Basic_015_382_20240521T224612_20240521T233658_PIC0_01.nc

% figure; imagesc(SWOT.ssha_karin_2{150} - SWOT.height_cor_xover{150}); colorbar

one_if_false_nan_if_true = @(IN) (~IN)./(~IN);

II = 173; % max = 180
for ii = [II + [0 1]]
    m_pcolor_centered(SWOT.lon{ii}, SWOT.lat{ii}, ...
        [SWOT.ssha_karin_2{ii} + SWOT.height_cor_xover{ii}] .* ...
        one_if_false_nan_if_true(SWOT.ssha_karin_2_qual{ii} + SWOT.height_cor_xover_qual{ii}) );
    % disp(['Time = ' datestr(mean(SWOT.time{ii})) ' (ii = ' num2str(ii) ')'])
    disp(['Time = ' datestr(SWOT.time{ii}(dsearchn(SWOT.lat{ii}(35,:)', 36)) ) ' (ii = ' num2str(ii) ')'])
    % ^ 36 = approx. xover lat., 35 = index of center of the swath
end
colormap(bwr(64)*0.99)
CB = colorbar;
CB.Label.String = 'SSHA (m)';
CB.Position = [0.5244    0.1620    0.0193    0.1935];
CB.FontSize = 20;
clim([-1 1]*0.2)
CB.Ticks = [-.2 -.1 0 .1 .2];

% % Box indicating in-set:
INSET_LAT = [35 37];
INSET_LON = [-126 -124.5];
collimate = @(IN) IN(:);
m_plot([INSET_LON flip(INSET_LON) , INSET_LON(1)]', ...
       [collimate([INSET_LAT ; INSET_LAT]) ; INSET_LAT(1)], 'k', 'LineWidth', 1)

% m_grid('box','fancy', 'backgroundcolor',[1 1 1]*0.9);
m_grid('box','fancy', 'backgroundcolor',[128 226 255]/255,'fontsize',18);

set(gcf,'Position',[262    69   962   728])

% exportgraphics(gcf,...
%     '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/component_figures/F1_USWC_map.pdf',...
%     'BackgroundColor','none','ContentType','vector')

%%      Small panel - background map

close all

figure
M_P = m_proj('Lambert',...
    'longitudes',INSET_LON,...
    'latitudes',INSET_LAT);
set(gcf,'color','w')
% COAST = m_gshhs_l('patch',0.8*[1 1 1]);
COAST = m_gshhs_l('patch',[93 165 38]/255);

hold on

% % Small panel - crossover
II = 173; % max = 180
for ii = [II + [0 1]]
    m_pcolor_centered(SWOT.lon{ii}, SWOT.lat{ii}, ...
        [SWOT.ssha_karin_2{ii} + SWOT.height_cor_xover{ii}] .* ...
        one_if_false_nan_if_true(SWOT.ssha_karin_2_qual{ii} + SWOT.height_cor_xover_qual{ii}) );
    % disp(['Time = ' datestr(mean(SWOT.time{ii})) ' (ii = ' num2str(ii) ')'])
    disp(['Time = ' datestr(SWOT.time{ii}(dsearchn(SWOT.lat{ii}(35,:)', 36)) ) ' (ii = ' num2str(ii) ')'])
    % ^ 36 = approx. xover lat., 35 = index of center of the swath
end
colormap(bwr(64)*0.999)
clim([-1 1]*0.2)
% % % % % 

m_grid('box','fancy', 'backgroundcolor',[128 226 255]/255, ...
       'xtick', [round(INSET_LON(1)):round(INSET_LON(end))], ...
       'ytick', [round(INSET_LAT(1)):round(INSET_LAT(end))], ...
       'fontsize',18)

% % Small panel - gliders

% % gliders
if ~exist('Glider_loc','var')
    G_DIR = dir('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/GLIDERS/ru*_L2_Processed_2*.nc');
    Glider_loc = {};
    for ii = 1:length(G_DIR)
        Glider_loc{ii} = [ncread([G_DIR(ii).folder '/' G_DIR(ii).name],'LONGITUDE_PROFILE') + 1i*ncread([G_DIR(ii).folder '/' G_DIR(ii).name],'LATITUDE_PROFILE')];

        % % % Old QC correction, no longer needed:
        % GliderQC = ncread([G_DIR(ii).folder '/' G_DIR(ii).name],'QC_FLAG');
        % Glider_loc{ii} = Glider_loc{ii}.*[~GliderQC./~GliderQC];
    end
else
end

GLIDER_COLOR = {[1 1 1]*0.6 , [0 1 1]*0.6};
for ii = 1:length(Glider_loc)
    m_plot(real(Glider_loc{ii}), imag(Glider_loc{ii}), '.-','Color',GLIDER_COLOR{ii})
end

% % Small panel - moorings
if ~exist('Mooring_nom_loc','var')
    M_P_DIR = dir('/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/PROF/*.nc');
    Mooring_nom_loc = nan(length(length(M_P_DIR)),2);
    for ii = 1:length(M_P_DIR)
        Mooring_nom_loc(ii,1) = ncread([M_P_DIR(ii).folder '/' M_P_DIR(ii).name],'LONGITUDE');
        Mooring_nom_loc(ii,2) = ncread([M_P_DIR(ii).folder '/' M_P_DIR(ii).name],'LATITUDE');
    end
else
end
for ii = 1:length(M_P_DIR)
    SP_sign = sign(strcmp(M_P_DIR(ii).name(33),'S') - 0.5);
    if SP_sign > 0 % S
        m_plot(Mooring_nom_loc(ii,1), Mooring_nom_loc(ii,2), '.k','MarkerSize',20)
        m_text(Mooring_nom_loc(ii,1) + 0.1, Mooring_nom_loc(ii,2) + 0.02,...
           M_P_DIR(ii).name(33:34), 'FontSize',16)
    else % P
        m_plot(Mooring_nom_loc(ii,1), Mooring_nom_loc(ii,2), '*k','MarkerSize',10,'LineWidth',2)
        m_text(Mooring_nom_loc(ii,1) - 0.2, Mooring_nom_loc(ii,2) - 0.02,...
           M_P_DIR(ii).name(33:34), 'FontSize',16)
    end
end
% m_plot(Mooring_nom_loc(:,1), Mooring_nom_loc(:,2), '.k','MarkerSize',20)

set(gcf,'Position',[255   207   693   520])


% exportgraphics(gcf,...
%     '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/component_figures/F1_xover_inset.pdf',...
%     'BackgroundColor','none','ContentType','vector')



%% Unused

% m_proj('oblique mercator',...
%     'longitudes',[-128 -118],...
%     'latitudes',[45 30],'direction','vertical','aspect',.5);