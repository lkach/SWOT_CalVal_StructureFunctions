%% Supplemental Material
%% Go through several slopes and compare spread of recovered slopes to true slope

addpath ~; lkaddpath

close all
% error

% SLOPES_TO_TEST = [[-4]:0.5:[-1]]';
% SLOPES_TO_TEST = [SLOPES_TO_TEST; [   -2.2 -2.1 -1.9 -1.8   ]'];
% SLOPES_TO_TEST = unique(SLOPES_TO_TEST);
SLOPES_TO_TEST = [[-4.5]:0.1:[-0.5]]';

% NN = 100000;
NN = 200000;
i0 = round(NN/2);
dx = .100; % km
XX = 0:dx:[dx*[NN - 1]]; % km
    XX = XX';
    XX = XX - XX(i0);
dk = 1/[XX(end) - XX(1)];
k_Ny = 1/(2*dx); % 1/200 cpkm
kk = dk:dk:k_Ny; kk = kk';

% Don't consider the ability for moorings to move at this point:
Mooring_i = dsearchn(XX,[0:10:100]');
X_mooring = XX(Mooring_i);
BIN_EDGES = [[-5]:10:[max(X_mooring)+5]]'; % m
BIN_CENTERS = [BIN_EDGES(1:[end-1])+BIN_EDGES(2:end)]/2;
N_mooring = length(Mooring_i);

N_t = 100;
MC = 100;

SLOPES_INFERRED_FROM_SF = nan(length(SLOPES_TO_TEST),MC);

warning('off')
for slope_i = 1:length(SLOPES_TO_TEST)
    SLOPE_true = SLOPES_TO_TEST(slope_i);
    % % % Spectrum with two slope:
    SS = kk.^SLOPE_true;
    SS = SS/sum(SS);
    for mc = 1:MC
        % % % Make data
        for ti = 1:N_t
            DD = randn_ts(SS); % from a spectrum
            DD = [DD; nan*ones(4,1)];
            % DD = AR_make(); % from an AR process
            StericHeight_t(:,ti) = DD(Mooring_i);
        end
        SF_t = nan(N_mooring,N_t);
        for ti = 1:N_t
            % Use the custom function "nanstructurefunction"
            [SF_i,~,~] = ...
                nanstructurefunction(...
                X_mooring, ...
                StericHeight_t(:,ti), ...
                BIN_EDGES, 2);
            SF_t(:,ti) = SF_i;
        end
        if MC == 1 % plot
        [FIT_PARAM,~,FIT_SF,~,~,~] = ...
            nonlinear_lsqf(mean(SF_t(2:end,:),2), BIN_CENTERS(2:end),...
            {'a*t.^p','a','p'},...
            [max( [mean(SF_t(end,:))/[100^(-SLOPE_true-1)] , 7.3407e-07] ), -SLOPE_true-1],...
            {'t.^p','a*log(t).*t.^p'},...
            0.0001,20,true);
        else % don't plot
        [FIT_PARAM,~,FIT_SF,~,~,~] = ...
            nonlinear_lsqf(mean(SF_t(2:end,:),2), BIN_CENTERS(2:end),...
            {'a*t.^p','a','p'},...
            [max( [mean(SF_t(end,:))/[100^(-SLOPE_true-1)] , 7.3407e-07] ), -SLOPE_true-1],...
            {'t.^p','a*log(t).*t.^p'},...
            0.0001,20);
        end
        SLOPES_INFERRED_FROM_SF(slope_i,mc) = -1 - FIT_PARAM(2);
        % disp([num2str(mc) '/' num2str(MC)])
    end
    disp(['Finished ' num2str(SLOPE_true)])
end
warning('on')

%% Plot results

close all

% figure
% for slope_i = 1:length(SLOPES_TO_TEST)
%     histogram(SLOPES_INFERRED_FROM_SF(slope_i,:), (-5):(0.1):(0)); hold on
% end

figure('Color','w')
HIST2 = histogram2( SLOPES_INFERRED_FROM_SF, repmat(SLOPES_TO_TEST,1,MC), ...
                   (-5.05):(0.1):(0), (-5.05):(0.1):(0),...
                   'HandleVisibility','off','Normalization','pdf');
[X_grid,Y_grid] = meshgrid([HIST2.XBinEdges(1:[end-1]) + HIST2.XBinEdges(2:end)]/2, ...
                           [HIST2.YBinEdges(1:[end-1]) + HIST2.YBinEdges(2:end)]/2);
HIST2_Values = HIST2.Values;
    HIST2_Values = [HIST2_Values.^2]./[HIST2_Values];
pcolor_centered(X_grid,Y_grid,HIST2_Values); shading flat; hold on
plot([min(X_grid(:)) max(X_grid(:))],...
     [min(Y_grid(:)) max(Y_grid(:))],'k')
axis square
set(gca,'XLim',[min(X_grid(:)) max(X_grid(:))],...
        'YLim',[min(Y_grid(:)) max(Y_grid(:))])
colormap(turbo)
CB = colorbar; CB.Label.String = 'PDF';
CB.Position = [0.7303 0.1611 0.0277 0.4000];

xlabel('True slopes')
ylabel('Slopes from S.F.')
set(gca,'FontSize',20)
GCF_P = get(gcf,'Position');
GCF_P(3) = 595;
set(gcf,'Position',GCF_P);

% Also print a table for the publication that summarizes this:
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('TRUE')
for slope_i = 1:length(SLOPES_TO_TEST)
    disp(num2str(SLOPES_TO_TEST(slope_i)))
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('MEDIAN')
for slope_i = 1:length(SLOPES_TO_TEST)
    disp(num2str(median(SLOPES_INFERRED_FROM_SF(slope_i,:),'omitnan'),'%.2f'))
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('MEAN')
for slope_i = 1:length(SLOPES_TO_TEST)
    % disp(num2str(mean(SLOPES_INFERRED_FROM_SF(slope_i,:),'omitnan'),'%.2f'))
    disp(num2str(mean(SLOPES_INFERRED_FROM_SF(slope_i,abs(SLOPES_INFERRED_FROM_SF(slope_i,:))<5),'omitnan'),'%.2f'))
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('STD')
for slope_i = 1:length(SLOPES_TO_TEST)
    % disp(num2str(std(SLOPES_INFERRED_FROM_SF(slope_i,:),'omitnan'),'%.2f'))
    disp(num2str(std(SLOPES_INFERRED_FROM_SF(slope_i,abs(SLOPES_INFERRED_FROM_SF(slope_i,:))<5),'omitnan'),'%.2f'))
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Converged fits')
for slope_i = 1:length(SLOPES_TO_TEST)
    disp(num2str(sum(isfinite(SLOPES_INFERRED_FROM_SF(slope_i,:)))))
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%

error

exportgraphics(gcf,...
['/Users/kachelein/Documents/JPL/papers/my_work/CalVal_StructureFunctions/figures/' ...
 'extra_Slope_TrueVsEstimated.pdf'],...
'BackgroundColor','none','ContentType','vector')

% save('~/MATLAB/11.2024/SLOPES_INFERRED_FROM_SF.mat')
