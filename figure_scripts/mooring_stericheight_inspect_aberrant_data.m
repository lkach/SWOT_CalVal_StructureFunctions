%% Load the workspace if needed

% error('Only load if this is needed.')

INPUT = input(['Do you want to load ~606 MB into memory? Enter any number for "yes" or push\n' ...
    '"Enter" with a blank input for "No".\n']);

if isempty(INPUT)
    disp(['No data loaded'])
else
    load(['/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/data/' ...
        'SF_plot_fixed_bottom_windowed_SWOTt_data.mat'])
    disp('Loaded data')
end

% Rearrange the spike indices to match the North to South order that the rest
% of this script follows:
SPIKE_IND_NtoS = SPIKE_IND;
SPIKE_IND_NtoS{1}  = SPIKE_IND{7}; % S1
SPIKE_IND_NtoS{2}  = SPIKE_IND{1}; % P1
SPIKE_IND_NtoS{3}  = SPIKE_IND{2}; % P2
SPIKE_IND_NtoS{4}  = SPIKE_IND{8}; % S2
SPIKE_IND_NtoS{5}  = SPIKE_IND{3}; % P3
% No P4
SPIKE_IND_NtoS{6}  = SPIKE_IND{9}; % S3
SPIKE_IND_NtoS{7}  = SPIKE_IND{4}; % P5
SPIKE_IND_NtoS{8}  = SPIKE_IND{5}; % P6
SPIKE_IND_NtoS{9}  = SPIKE_IND{10}; % S4
SPIKE_IND_NtoS{10} = SPIKE_IND{6}; % P7

%%

close all

figure
plot(S1.P.TIME_STERIC_HEIGHT,S1_sh_600,'.-'); hold on
plot(S2.P.TIME_STERIC_HEIGHT,S2_sh_600,'.-')
plot(S3.P.TIME_STERIC_HEIGHT,S3_sh_600,'.-')
plot(S4.P.TIME_STERIC_HEIGHT,S4_sh_600,'.-')

plot(P1.P.TIME_STERIC_HEIGHT,P1_sh_600,'.-')
plot(P2.P.TIME_STERIC_HEIGHT,P2_sh_600,'.-')
plot(P3.P.TIME_STERIC_HEIGHT,P3_sh_600,'.-')
plot(P5.P.TIME_STERIC_HEIGHT,P5_sh_600,'.-')
plot(P6.P.TIME_STERIC_HEIGHT,P6_sh_600,'.-')
plot(P7.P.TIME_STERIC_HEIGHT,P7_sh_600,'.-')


figure
plot(S1.P.TIME_STERIC_HEIGHT,S1_sh_600 + 0*0.1*-1,'.-'); hold on
plot(S2.P.TIME_STERIC_HEIGHT,S2_sh_600 + 3*0.1*-1,'.-')
plot(S3.P.TIME_STERIC_HEIGHT,S3_sh_600 + 6*0.1*-1,'.-')
plot(S4.P.TIME_STERIC_HEIGHT,S4_sh_600 + 9*0.1*-1,'.-')

plot(P1.P.TIME_STERIC_HEIGHT,P1_sh_600 + 1*0.1*-1,'.-')
plot(P2.P.TIME_STERIC_HEIGHT,P2_sh_600 + 2*0.1*-1,'.-')
plot(P3.P.TIME_STERIC_HEIGHT,P3_sh_600 + 4*0.1*-1,'.-')
plot(P5.P.TIME_STERIC_HEIGHT,P5_sh_600 + 7*0.1*-1,'.-')
plot(P6.P.TIME_STERIC_HEIGHT,P6_sh_600 + 8*0.1*-1,'.-')
plot(P7.P.TIME_STERIC_HEIGHT,P7_sh_600 + 10*0.1*-1,'.-')

% %
% close all
figure
HBINS = [-.1:.001:.1];
histogram(S1_sh_600,HBINS,'Normalization','pdf'); hold on
histogram(S2_sh_600,HBINS,'Normalization','pdf'); hold on
histogram(S3_sh_600,HBINS,'Normalization','pdf'); hold on
histogram(S4_sh_600,HBINS,'Normalization','pdf'); hold on

histogram(P1_sh_600,HBINS,'Normalization','pdf'); hold on
histogram(P2_sh_600,HBINS,'Normalization','pdf'); hold on
histogram(P3_sh_600,HBINS,'Normalization','pdf'); hold on
histogram(P5_sh_600,HBINS,'Normalization','pdf'); hold on
histogram(P6_sh_600,HBINS,'Normalization','pdf'); hold on
histogram(P7_sh_600,HBINS,'Normalization','pdf'); hold on

%% Difference between all pairs

close all

MM_list = {'S1','P1','P2','S2','P3',  'S3','P5','P6','S4','P7'};

T0 = datenum('1950-01-01');
T_start = datenum('2023-04-01') - T0;
T_end   = datenum('2023-07-10') - T0;

Pseudo_location = [0 10 20 30 40   60 70 80 90 100];
% ^ along a straight line going south along the descending track

Pair_indices = [];
SH_diff = {};
II = 1;
for i1 = 1:length(MM_list)
    for i2 = i1:length(MM_list)
        Pair_indices = [Pair_indices ; i1, i2];
        P_IND_A = eval([ MM_list{i1} '.P.TIME_STERIC_HEIGHT > T_start & ' MM_list{i1} '.P.TIME_STERIC_HEIGHT < T_end & ' ... <-- Time interval
           '[' MM_list{i1} '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM_list{i1} '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM_list{i1}(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
           '~SPIKE_IND_NtoS{i1}']); % <-- Spike removal
        P_IND_B = eval([ MM_list{i2} '.P.TIME_STERIC_HEIGHT > T_start & ' MM_list{i2} '.P.TIME_STERIC_HEIGHT < T_end & ' ... <-- Time interval
           '[' MM_list{i2} '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM_list{i2} '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM_list{i2}(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
           '~SPIKE_IND_NtoS{i2}']); % <-- Spike removal
        SH_A = eval([MM_list{i1} '_sh_600(P_IND_A)']);
        SH_B = eval([MM_list{i2} '_sh_600(P_IND_B)']);
        T_A = eval([MM_list{i1} '.P.TIME_STERIC_HEIGHT(P_IND_A)']); % units of days
        T_B = eval([MM_list{i2} '.P.TIME_STERIC_HEIGHT(P_IND_B)']);
        if length(SH_A) >= length(SH_B)
            SH_long = SH_A;
            SH_shrt = SH_B;
            T_long = T_A;
            T_shrt = T_B;
        else
            SH_long = SH_B;
            SH_shrt = SH_A;
            T_long = T_B;
            T_shrt = T_A;
        end
        dt_max_allowed = 2/24; % i.e. must be within this many days to be valid
        SH_diff_12 = nan(size(SH_long));
        % for ti=1:length(SH_long) % match the SHORTER data to the times of the LONGER data
        %     diff_squared_ti = SH_long(ti) - ...
        %                       SH_shrt(dsearchn(T_shrt,T_long(ti)));
        %     if [ [T_shrt(dsearchn(T_shrt,T_long(ti))) - T_long(ti)] <= dt_max_allowed ] && ... both time series are close in time
        %        [ T_long(ti) >= T_start ] && [ T_long(ti) <= T_end ] % within SWOT cal/val
        %         SH_diff_12(ti) = diff_squared_ti;
        %     else
        %     end
        % end
        for ti=1:length(SH_shrt) % match the LONGER data to the times of the SHORTER data
            diff_squared_ti = SH_shrt(ti) - ...
                              SH_long(dsearchn(T_long,T_shrt(ti)));
            if [ [T_long(dsearchn(T_long,T_shrt(ti))) - T_shrt(ti)] <= dt_max_allowed ] && ... both time series are close in time
               [ T_shrt(ti) >= T_start ] && [ T_shrt(ti) <= T_end ] % within SWOT cal/val
                SH_diff_12(ti) = diff_squared_ti;
            else
            end
        end
        SH_diff{II} = SH_diff_12;
        II = II + 1;
    end
end

%% Visualize matrices of pertinent quantities

close all

blank_mat_10x10 = zeros(10,10);
[XX,YY] = meshgrid(1:10,1:10);

var_mat_10x10  = nan(10,10);
mean_mat_10x10 = nan(10,10);
meanSqDiff_mat_10x10 = nan(10,10);
Pseudo_distance = abs(Pseudo_location' - Pseudo_location);
for ii = 1:size(Pair_indices,1)
    var_mat_10x10(Pair_indices(ii,1),Pair_indices(ii,2)) = ...
        var(SH_diff{ii},'omitmissing');
    mean_mat_10x10(Pair_indices(ii,1),Pair_indices(ii,2)) = ...
        mean(SH_diff{ii},'omitmissing');
    meanSqDiff_mat_10x10(Pair_indices(ii,1),Pair_indices(ii,2)) = ...
        mean(SH_diff{ii}.^2,'omitmissing');
end

% figure
% % imagesc(blank_mat_10x10); grid on
% imagesc(var_mat_10x10); grid on
% % pcolor_centered(XX,YY,var_mat_10x10); grid on
% xticklabels(MM_list)
% yticklabels(MM_list)


figure
IM = imagesc(var_mat_10x10); grid on
xticklabels(MM_list)
yticklabels(MM_list)
clim([-1 1]*max(abs(IM.CData(:))))
colormap('bwr')
title('Covariance of difference')

figure
IM = imagesc(mean_mat_10x10); grid on
xticklabels(MM_list)
yticklabels(MM_list)
clim([-1 1]*max(abs(IM.CData(:))))
colormap('bwr')
title('Mean difference')

figure
IM = imagesc(meanSqDiff_mat_10x10); %grid on
xticklabels(MM_list)
yticklabels(MM_list)
clim([-1 1]*max(abs(IM.CData(:))))
colormap('bwr')
title('Mean squared difference')

%% Regular old means

% mean(P1.P.STERIC_HEIGHT_ANOMALY,'omitmissing')

close all

SHA_mean = nan(length(MM_list),1);
for ii = 1:length(MM_list)
    SHA_mean(ii) = mean(eval([MM_list{ii} '_sh_600']),'omitmissing');
end

figure
plot([1:length(MM_list)]', SHA_mean, '.-')
xticklabels(MM_list)
ylabel('mean(steric height anomaly to 600 m)')

%% Mean steric height in the cal/val interval

T_start = datenum('2023-04-01') - T0;
T_end   = datenum('2023-07-10') - T0;

% T_start = datenum('2023-05-08') - T0;
% T_end   = datenum('2023-06-21') - T0;

% close all

SHA_mean = nan(length(MM_list),1);
SHA_std = nan(length(MM_list),1);
for ii = 1:length(MM_list)
    P_IND = eval([ MM_list{ii} '.P.TIME_STERIC_HEIGHT > T_start & ' MM_list{ii} '.P.TIME_STERIC_HEIGHT < T_end & ' ... <-- Time interval
        '[' MM_list{ii} '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM_list{ii} '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM_list{ii}(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
        '~SPIKE_IND_NtoS{ii}']); % <-- Spike removal
    SHA_mean(ii) = mean(eval([MM_list{ii} '_sh_600(P_IND)']),'omitmissing');
    SHA_std(ii) = std(eval([MM_list{ii} '_sh_600(P_IND)']),'omitmissing');
end

figure('Color','w')
% plot([1:length(MM_list)]', 100*SHA_mean, 'k.-','LineWidth',1.5,'MarkerSize',20)
errorbar(Pseudo_location, 100*SHA_mean, 100*SHA_std, 'k.-','LineWidth',1.5,'MarkerSize',20); hold on

% Optional for evaluation: linear fit.
plot(Pseudo_location, [ones(size(Pseudo_location')),Pseudo_location']*[[ones(size(Pseudo_location')),Pseudo_location']\[100*SHA_mean]],'r.-')

MM_list_withP4 = {'S1','P1','P2','S2','P3', '(P4)' ,'S3','P5','P6','S4','P7'};;
xticklabels(MM_list_withP4)
title('Mean of 600m steric height anomaly')
ylabel('cm')

set(gca,'FontSize',16, ...
        ...'XLim',[1 10],...
        'YLim',[-10 3])
% set(gcf,'Position',[-1439         622         376         189])
% set(gcf,'Position',[1   525   514   272])
set(gcf,'Position',[1   525   551   272])

INPUT = input(['Do you want to save this figure? Enter any number for "yes" or\n push ' ...
    '"Enter" with a blank input for "No".\n']);

if isempty(INPUT)
    disp(['Image not saved'])
else
    figure(1)
    % exportgraphics(gcf,...
    %     '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/extra_FigG.pdf',...
    %     'BackgroundColor','none','ContentType','vector')
    % disp(['Image saved'])

    % Post-revisions:
    exportgraphics(gcf,...
        '/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/REV1/extra_FigG.pdf',...
        'BackgroundColor','none','ContentType','vector')
    disp(['Image saved'])
end

%% Data in SWOT cal/val window

close all

T1 = datenum('2023-03-28');

figure
plot(T0 + S1.P.TIME_STERIC_HEIGHT(S1.P.TIME_STERIC_HEIGHT >= T_start & S1.P.TIME_STERIC_HEIGHT <= T_end), ...
    S1_sh_600(S1.P.TIME_STERIC_HEIGHT >= T_start & S1.P.TIME_STERIC_HEIGHT <= T_end) + 0*0.1*-1,'.-'); hold on
plot(T0 + S2.P.TIME_STERIC_HEIGHT(S2.P.TIME_STERIC_HEIGHT >= T_start & S2.P.TIME_STERIC_HEIGHT <= T_end), ...
    S2_sh_600(S2.P.TIME_STERIC_HEIGHT >= T_start & S2.P.TIME_STERIC_HEIGHT <= T_end) + 3*0.1*-1,'.-')
plot(T0 + S3.P.TIME_STERIC_HEIGHT(S3.P.TIME_STERIC_HEIGHT >= T_start & S3.P.TIME_STERIC_HEIGHT <= T_end), ...
    S3_sh_600(S3.P.TIME_STERIC_HEIGHT >= T_start & S3.P.TIME_STERIC_HEIGHT <= T_end) + 6*0.1*-1,'.-')
plot(T0 + S4.P.TIME_STERIC_HEIGHT(S4.P.TIME_STERIC_HEIGHT >= T_start & S4.P.TIME_STERIC_HEIGHT <= T_end), ...
    S4_sh_600(S4.P.TIME_STERIC_HEIGHT >= T_start & S4.P.TIME_STERIC_HEIGHT <= T_end) + 9*0.1*-1,'.-')
text(T1,0*0.1*-1,'S1'); text(T1,3*0.1*-1,'S2'); text(T1,6*0.1*-1,'S3'); text(T1,9*0.1*-1,'S4')

plot(T0 + P1.P.TIME_STERIC_HEIGHT(P1.P.TIME_STERIC_HEIGHT >= T_start & P1.P.TIME_STERIC_HEIGHT <= T_end), ...
    P1_sh_600(P1.P.TIME_STERIC_HEIGHT >= T_start & P1.P.TIME_STERIC_HEIGHT <= T_end) + 1*0.1*-1,'.-')
plot(T0 + P2.P.TIME_STERIC_HEIGHT(P2.P.TIME_STERIC_HEIGHT >= T_start & P2.P.TIME_STERIC_HEIGHT <= T_end), ...
    P2_sh_600(P2.P.TIME_STERIC_HEIGHT >= T_start & P2.P.TIME_STERIC_HEIGHT <= T_end) + 2*0.1*-1,'.-')
plot(T0 + P3.P.TIME_STERIC_HEIGHT(P3.P.TIME_STERIC_HEIGHT >= T_start & P3.P.TIME_STERIC_HEIGHT <= T_end), ...
    P3_sh_600(P3.P.TIME_STERIC_HEIGHT >= T_start & P3.P.TIME_STERIC_HEIGHT <= T_end) + 4*0.1*-1,'.-')
plot(T0 + P5.P.TIME_STERIC_HEIGHT(P5.P.TIME_STERIC_HEIGHT >= T_start & P5.P.TIME_STERIC_HEIGHT <= T_end), ...
    P5_sh_600(P5.P.TIME_STERIC_HEIGHT >= T_start & P5.P.TIME_STERIC_HEIGHT <= T_end) + 7*0.1*-1,'.-')
plot(T0 + P6.P.TIME_STERIC_HEIGHT(P6.P.TIME_STERIC_HEIGHT >= T_start & P6.P.TIME_STERIC_HEIGHT <= T_end), ...
    P6_sh_600(P6.P.TIME_STERIC_HEIGHT >= T_start & P6.P.TIME_STERIC_HEIGHT <= T_end) + 8*0.1*-1,'.-')
plot(T0 + P7.P.TIME_STERIC_HEIGHT(P7.P.TIME_STERIC_HEIGHT >= T_start & P7.P.TIME_STERIC_HEIGHT <= T_end), ...
    P7_sh_600(P7.P.TIME_STERIC_HEIGHT >= T_start & P7.P.TIME_STERIC_HEIGHT <= T_end) + 10*0.1*-1,'.-')
text(T1,1*0.1*-1,'P1');text(T1,2*0.1*-1,'P2');text(T1,4*0.1*-1,'P3');text(T1,7*0.1*-1,'P5');
text(T1,8*0.1*-1,'P6');text(T1,10*0.1*-1,'P7');                   text(T1,5*0.1*-1,'(P4)');
yticklabels({})
datetick('x')

figure
PLOTTED_MOORINGS = 1:length(MM_list);[2 6];
for ii = PLOTTED_MOORINGS % 1:length(MM_list)
    P_IND = eval([ MM_list{ii} '.P.TIME_STERIC_HEIGHT > T_start & ' MM_list{ii} '.P.TIME_STERIC_HEIGHT < T_end & ' ... <-- Time interval
        '[' MM_list{ii} '.P.STERIC_HEIGHT_MAX_PROF_DEPTH - ' MM_list{ii} '.P.STERIC_HEIGHT_MIN_PROF_DEPTH] > ' MM_list{ii}(1) 'P_RANGE_LIM & ' ... <-- Profiler depth requirement
        '~SPIKE_IND_NtoS{ii}']); % <-- Spike removal
    plot(T0 + eval([MM_list{ii} '.P.TIME_STERIC_HEIGHT(P_IND)']), eval([MM_list{ii} '_sh_600(P_IND)']), '.-'); hold on
end
datetick('x')
legend(MM_list{PLOTTED_MOORINGS})