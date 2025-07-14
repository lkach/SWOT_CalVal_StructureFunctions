%% Nonlinear Least Squares Fit to obtain the values in T1

warning('First run F5.m to load needed variables')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


if LEVEL == 2
% % % SP.600m.24hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (600m, 36hr windowed)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_36hrwin,Separation_S_P_600m_36hrwin,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_P_600m_36hrwin, m2_to_cm2*StructFunc_S_P_600m_36hrwin, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:[end-1]),BIN_CENTR(2:[end-1]),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:[end-1])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/36hr_win/M to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
% [NLLSF_COEF,~,~,~] = ...
%             nonlinear_lsqf(SF_avg(5:[end-1]),BIN_CENTR(5:[end-1]),...
%             FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
% disp(['600m24hr/M 40 to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1)])
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:5),BIN_CENTR(2:5),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:5)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/36hr_win/M to ' num2str(BIN_CENTR(5)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')


% % % SP.600m.1hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (600m, 1hr)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_1hr,Separation_S_P_600m_1hr,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_P_600m_1hr, m2_to_cm2*StructFunc_S_P_600m_1hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:[end-1]),BIN_CENTR(2:[end-1]),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:[end-1])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/1hr/M to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
% [NLLSF_COEF,~,~,~] = ...
%             nonlinear_lsqf(SF_avg(5:[end-1]),BIN_CENTR(5:[end-1]),...
%             FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
% disp(['600m/1hr/M 40 to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1)])
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:5),BIN_CENTR(2:5),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:5)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/1hr/M to ' num2str(BIN_CENTR(5)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')


% % % SPG.600m.1hr
LEG_TEXT = {LEG_TEXT{:},'Moorings+Gliders (600m, 1hr)'};
BIN_EDGES = [(-dBIN/2/2):[dBIN/2]:[14*dBIN/2/2]]; BIN_EDGES =  [BIN_EDGES, [BIN_EDGES(end)+dBIN]:[dBIN]:1.1]*(0.92)*DD;% initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_G_600m_1hr,Separation_S_P_G_600m_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_P_G_600m_1hr, m2_to_cm2*StructFunc_S_P_G_600m_1hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:[end-1]),BIN_CENTR(2:[end-1]),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:[end-1])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/1hr/MG to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
% [NLLSF_COEF,~,~,~] = ...
%             nonlinear_lsqf(SF_avg(8:[end-1]),BIN_CENTR(8:[end-1]),...
%             FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
% disp(['600m24hr/MG 40 to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1)])
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:8),BIN_CENTR(2:8),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:8)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/1hr/MG to ' num2str(BIN_CENTR(8)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')




% % % SP.1800m.36hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (1800m, 36hr windowed)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_deep_36hrwin, Separation_S_deep_36hrwin,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_deep_36hrwin, m2_to_cm2*StructFunc_S_deep_36hrwin, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:end),BIN_CENTR(2:end),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:end)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['1800m/36hr_win/M to ' num2str(BIN_CENTR(end)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
% 
BIN_CENTR = interval_avg(Separation_S_deep_1hr, Separation_S_deep_1hr,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_deep_1hr, m2_to_cm2*StructFunc_S_deep_1hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:end),BIN_CENTR(2:end),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:end)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['1800m/1hr/M to ' num2str(BIN_CENTR(end)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')
else
end

% % % SWOT_at_S_P_G
LEG_TEXT = {LEG_TEXT{:},'SWOT cal/val at Moorings+Gliders'};
BIN_EDGES = [(-dBIN/2):[dBIN]:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_SWOT_moor_12hr,Separation_SWOT_moor_12hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_SWOT_moor_12hr, m2_to_cm2*StructFunc_SWOT_moor_12hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:[end-1]),BIN_CENTR(2:[end-1]),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:[end-1])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
% The covariance of gamma (SF order) is the same as that of lambda
% (spectral slope) because both are just the other times -1 and shifted -1.
disp(['SWOT to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
% [NLLSF_COEF,~,~,~] = ...
%             nonlinear_lsqf(SF_avg(5:[end-1]),BIN_CENTR(5:[end-1]),...
%             FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
% disp(['SWOT 40 to ' num2str(BIN_CENTR([end-1])) ' km: ' num2str(-NLLSF_COEF(2)-1)])
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg(2:5),BIN_CENTR(2:5),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg(2:5)' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['SWOT to ' num2str(BIN_CENTR(5)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')

[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(Mean_SF_MAT_swot(2:dsearchn(DIST_VEC,90)), DIST_VEC(2:dsearchn(DIST_VEC,90)),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ Mean_SF_MAT_swot(2:dsearchn(DIST_VEC,90)) - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['SWOT (all) to ' num2str(DIST_VEC(dsearchn(DIST_VEC,90))) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(Mean_SF_MAT_swot(2:dsearchn(DIST_VEC,40)), DIST_VEC(2:dsearchn(DIST_VEC,40)),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ Mean_SF_MAT_swot(2:dsearchn(DIST_VEC,40)) - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['SWOT (all) to ' num2str(DIST_VEC(dsearchn(DIST_VEC,40))) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')

% 600m/36hr_win/M to 90.6818 km: -2.5145 +- 0.13978
% 600m/36hr_win/M to 40.9595 km: -2.3272 +- 0.065636
% 
% 600m/1hr/M to 90.8624 km: -2.4122 +- 0.13862
% 600m/1hr/M to 40.9355 km: -2.0894 +- 0.011356
% 
% 600m/1hr/MG to 91.0229 km: -2.4328 +- 0.10239
% 600m/1hr/MG to 40.3921 km: -1.9808 +- 0.11839
% 
% 1800m/36hr_win/M to 90.7835 km: -2.1807 +- 0.3001
% 1800m/1hr/M to 90.7901 km: -2.0682 +- 0.30687
% 
% SWOT to 90.9635 km: -2.3211 +- 0.077906
% SWOT to 40.8199 km: -2.3662 +- 0.039596
% 
% SWOT (all) to 90.1931 km: -1.9573 +- 0.020517
% SWOT (all) to 40.0858 km: -2.2172 +- 0.0050224

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spectral slope of average SFs from above, but with 50 km spike removed:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% /Users/kachelein/Documents/MATLAB/Luke/GaussNewtonParameterUncertaintyTest.m
% https://www8.cs.umu.se/kurser/5DA001/HT07/lectures/lsq-handouts.pdf
% https://www.incertitudes.fr/book.pdf

warning([''])

clc ; close all

FitBasisFunction = {'a*t.^b','a','b'};
FitParameters0 = [0.1,1.5];
FitBasisDerivatives = {'t.^b','a*log(t).*t.^b'};
ToleranceLevel = 0.00001; % NLLSF parameter, tells change in parameter magnitude, if small then the result was stable
N_nllsf_iterations = 10; % number of iterations before we give up




% % % SP.600m.24hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (600m, 36hr windowed)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_36hrwin,Separation_S_P_600m_36hrwin,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_P_600m_36hrwin, m2_to_cm2*StructFunc_S_P_600m_36hrwin, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg([2:5,7:(end-1)]),BIN_CENTR([2:5,7:(end-1)]),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg([2:5,7:(end-1)])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/36hr_win/M to ' num2str(BIN_CENTR(end-1)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')


% % % SP.600m.1hr
LEG_TEXT = {LEG_TEXT{:},'Moorings (600m, 1hr)'};
BIN_EDGES = [(-dBIN/2):dBIN:1.1]*(0.92)*DD; % initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_600m_1hr,Separation_S_P_600m_1hr,BIN_EDGES); BIN_CENTR = BIN_CENTR(isfinite(BIN_CENTR));
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_P_600m_1hr, m2_to_cm2*StructFunc_S_P_600m_1hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg([2:5,7:(end-1)]),BIN_CENTR([2:5,7:(end-1)]),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg([2:5,7:(end-1)])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/1hr/M to ' num2str(BIN_CENTR(end-1)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')


% % % SPG.600m.1hr
LEG_TEXT = {LEG_TEXT{:},'Moorings+Gliders (600m, 1hr)'};
BIN_EDGES = [(-dBIN/2/2):[dBIN/2]:[14*dBIN/2/2]]; BIN_EDGES =  [BIN_EDGES, [BIN_EDGES(end)+dBIN]:[dBIN]:1.1]*(0.92)*DD;% initial bin edges
BIN_CENTR = interval_avg(Separation_S_P_G_600m_1hr,Separation_S_P_G_600m_1hr,BIN_EDGES);
BIN_EDGES = [[ BIN_CENTR(1) - [BIN_CENTR(2) - BIN_CENTR(1)]/2 ], ...
               [BIN_CENTR(2:end) + BIN_CENTR(1:[end-1])]/2, ...
               [BIN_CENTR(end) + [BIN_CENTR(end) - BIN_CENTR(end-1)]/2] ];
SF_avg = interval_avg(Separation_S_P_G_600m_1hr, m2_to_cm2*StructFunc_S_P_G_600m_1hr, BIN_EDGES);
[NLLSF_COEF,JACOBIAN,SF_FIT,~] = ...
            nonlinear_lsqf(SF_avg([2:8,10:(end-1)]),BIN_CENTR([2:8,10:(end-1)]),...
            FitBasisFunction,FitParameters0,FitBasisDerivatives,ToleranceLevel,N_nllsf_iterations,[]);
NLLSF_COEF_cov = sqrt(inv(JACOBIAN'*JACOBIAN))*sqrt(sum([ SF_avg([2:5,7:(end-1)])' - SF_FIT ].^2)/[length(SF_FIT)-length(NLLSF_COEF)]);
disp(['600m/1hr/MG to ' num2str(BIN_CENTR(end-1)) ' km: ' num2str(-NLLSF_COEF(2)-1) ' +- ' num2str(NLLSF_COEF_cov(2,2))])
disp(' ')

% % % Slopes without the 50km bump:
% 
% 600m/36hr_win/M to 90.6818 km: -2.675 +- 0.071553
% 
% 600m/1hr/M to 90.8624 km: -2.5479 +- 0.098005
% 
% 600m/1hr/MG to 91.0229 km: -2.4988 +- 0.1875








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

