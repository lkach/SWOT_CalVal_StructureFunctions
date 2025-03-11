%% Evaluate the effect of a chosen time-domain window on wavenumber spectral slope

% S&F24 is the abbreviation for:
% Models of the Sea Surface Height Expression of the Internal-Wave Continuum
% R. M. Samelson and J. T. Farrar (2024)
% https://doi.org/10.1175/JPO-D-23-0178.1

% Numerically integrate the coordinate transfer integral from S&F24,
% equation 14, after taking into account equations 8-12, but using a
% transfer function (square of the FT of a window function) instead of the
% spectrum

% Constants are as follows:
ff = 1.348*10^-5;
% ^ cycles per second, corresponding to 20.6 hours, the inertial frequency
% of 35.5 degrees latitude (close to the center of the mooring array)
ff = ff*2*pi; % radians per second to make work with omega
c1 = 2.3;
% ^ m/s, as given by figure 2 in Chelton et al. (1998) and referenced in
% S&F24. For higher modes, they approximate c_j = c_1/j^2

%% Omega_j(K) = sqrt(c_j^2 K^2 + f^2)

jj = 1;
cj = c1/jj^2;

% That forumla (eq. 9) implies that K is the angular wavenumber, i.e.
% K = 2*pi/wavelength

% This also means that
% % K = sqrt([omega^2 - f^2]/c_j^2)

% omega = 2*pi*f
df = 1/[90*24]; % cph (over 90 days)
df = df*[1/3600]; % cps
f_vec = df:df:(1/3600); f_vec = f_vec';
omega = 2*pi*f_vec;
KK_j = sqrt([omega.^2 - ff^2]/cj^2);

close all
figure
plot(real(KK_j),'.-'); hold on
plot(imag(KK_j),'.-');

% Real K's start at omega = ff
KK_j = real(KK_j);
KK_j = KK_j(KK_j>0);

figure
plot([2*pi./KK_j]/1000,'.-')
title('Wavelength in km')

Omega_j = sqrt(cj^2 * KK_j.^2 + ff^2);

%% Define the transfer function

% Transfer function = square of the FT of the convolution window

% For a boxcar of width of T (note: zeta = ordinary frequency to avoid confusion with f):
% w(t) = Pi(t/T)
% w(zeta) = sin(pi T zeta)/(pi zeta)
% w(omega) = 2 sin(T omega/2)/omega
% W(omega) = 4 sin^2(T omega/2)/omega^2

WINDOW_WIDTH = 36*60*60; % in seconds, here equivalent of 36 hours
TRANFER_FUNCTION = @(OM) [2*sin([WINDOW_WIDTH]*OM/2)./OM].^2;
WW = TRANFER_FUNCTION(Omega_j);
% WW = [2*sin([WINDOW_WIDTH]*Omega_j/2)./Omega_j].^2;


% For a Hann of width of T (note: zeta = ordinary frequency to avoid confusion with f):
% w(t) = cos(pi t / T)^2 / T
% w(zeta) = sin(pi T zeta)/(2 pi zeta T (1 - T^2 zeta^2))
% w(omega) =  sin(T omega/2)/(omega T (1 - (T omega / [2 pi])^2))
% W(omega) = [sin(T omega/2)/(omega T (1 - (T omega / [2 pi])^2))]^2

% WINDOW_WIDTH = 36*60*60; % in seconds, here equivalent of 36 hours
% TRANFER_FUNCTION = @(OM) [sin(WINDOW_WIDTH*OM/2)./(OM*WINDOW_WIDTH.*(1 - (WINDOW_WIDTH*OM/[2*pi]).^2))].^2;
% WW = TRANFER_FUNCTION(Omega_j);

figure
subplot(211)
plot([(-WINDOW_WIDTH/2):(WINDOW_WIDTH/100):(WINDOW_WIDTH/2)], ... For Hann
     cos(pi*[(-WINDOW_WIDTH/2):(WINDOW_WIDTH/100):(WINDOW_WIDTH/2)] / WINDOW_WIDTH).^2 / WINDOW_WIDTH,'.-')
subplot(212)
plot(Omega_j,WW,'.-')

%% Plot the integrand for different values of lowercase k (k_test)

k_test = KK_j(10);
INTEGRAND = ...
    (2/pi)*WW.*(cj^2 .* KK_j)./[sqrt(cj^2 * KK_j.^2 + ff^2).*sqrt(KK_j.^2 - k_test^2)];
INTEGRAND = real(INTEGRAND);

close all
figure
plot(KK_j,INTEGRAND,'.-'); hold on
plot([1 1]*k_test, [0 1]*max(INTEGRAND(INTEGRAND<max(INTEGRAND))), 'k--')
xlabel('angular wavenumber (rad/m)')
ylabel('Integrand of eq. 14')
title(['for $k$ = ' num2str(k_test) ' rad/m'],'Interpreter','latex')

%% Do the same, but with more densely-packed wavenumbers as one approaches k_test from the right

k_test = KK_j(1000);
KK_j_int = unique([KK_j; [KK_j(dsearchn(KK_j,k_test)):...
                          ([KK_j(dsearchn(KK_j,k_test)+1) - KK_j(dsearchn(KK_j,k_test))]/20):...
                          KK_j(dsearchn(KK_j,k_test)+1)]']);

Omega_j_int = sqrt(cj^2 * KK_j_int.^2 + ff^2);
WW = TRANFER_FUNCTION(Omega_j_int);
INTEGRAND = ...
    (2/pi)*WW.*(cj^2 .* KK_j_int)./[sqrt(cj^2 * KK_j_int.^2 + ff^2).*sqrt(KK_j_int.^2 - k_test^2)];
INTEGRAND = real(INTEGRAND);

close all
figure
plot(KK_j_int,INTEGRAND,'.-'); hold on
plot([1 1]*k_test, [0 1]*max(INTEGRAND(INTEGRAND<max(INTEGRAND))), 'k--')
xlabel('angular wavenumber (rad/m)')
ylabel('Integrand of eq. 14')
title(['for $k$ = ' num2str(k_test) ' rad/m'],'Interpreter','latex')

%% Now go through valid lowercase k to get the full coordinate transformed,
%% 1D transfer function we seek:

kk = KK_j;
W1_j = nan(size(kk)); % The 1D transfer function we seek

for ik = 1:[length(kk)-1]
    k_test = kk(ik);
    KK_j_int = unique([KK_j; [KK_j(dsearchn(KK_j,k_test)):...
                              ([KK_j(dsearchn(KK_j,k_test)+1) - KK_j(dsearchn(KK_j,k_test))]/20):...
                              KK_j(dsearchn(KK_j,k_test)+1)]']);
    
    Omega_j_int = sqrt(cj^2 * KK_j_int.^2 + ff^2);
    WW = TRANFER_FUNCTION(Omega_j_int);
    INTEGRAND = ...
        (2/pi)*WW.*(cj^2 .* KK_j_int)./[sqrt(cj^2 * KK_j_int.^2 + ff^2).*sqrt(KK_j_int.^2 - k_test^2)];
    INTEGRAND = real(INTEGRAND);
    % INTEGRAND(INTEGRAND==Inf) = max(INTEGRAND(isfinite(INTEGRAND))); % one possible correction
    INTEGRAND(INTEGRAND==Inf) = 0; % another possible correction
    W1_j(ik) = trapz(KK_j_int,INTEGRAND);
end
%% plot

close all

figure('Color','w')
plot(kk,W1_j,'k.-'); hold on
plot(kk,kk.^-2,'r-')
plot(kk,[10^-10]*kk.^-2,'r-')
xlabel('angular wavenumber (rad/m)')
ylabel(['$W^1_{j=' num2str(jj) '}(k)$'],'Interpreter','latex')
set(gca,'XScale','log','YScale','log')

figure('Color','w')
plot(0.001*2*pi./kk,W1_j,'k.-'); hold on
plot(0.001*2*pi./kk,0.01*[2*pi./kk].^2,'r-')
plot(0.001*2*pi./kk,[10^-34]*[2*pi./kk].^6,'r-')
xlabel('Wavelength (km)')
ylabel(['$W^1_{j=' num2str(jj) '}(k)$'],'Interpreter','latex')
set(gca,'XScale','log','YScale','log')

%%
%%
%% Calculate the estimated transfer function from application to fake data


% % Rectangular Window:
% WINDOW_WIDTH = 36; % in hours
% TT = [-WINDOW_WIDTH/2]:1:[WINDOW_WIDTH/2];
% WINDOW = ones(WINDOW_WIDTH+0,1)/WINDOW_WIDTH;
% TRANFER_FUNCTION = @(OM) [2*sin([WINDOW_WIDTH]*OM/2)./(OM*WINDOW_WIDTH)].^2;

% % % For a Hann of width of T (note: zeta = ordinary frequency to avoid confusion with f):
% % % w(t) = cos(pi t / T)^2 / T
% % % w(zeta) = sin(pi T zeta)/(2 pi zeta T (1 - T^2 zeta^2))
% % % w(omega) =  sin(T omega/2)/(omega T (1 - (T omega / [2 pi])^2))
% % % W(omega) = [sin(T omega/2)/(omega T (1 - (T omega / [2 pi])^2))]^2

% % Hann Window:
WINDOW_WIDTH = 36; % in hours
TT = [-WINDOW_WIDTH/2]:1:[WINDOW_WIDTH/2];
WINDOW = [cos(pi*TT/WINDOW_WIDTH).^2]/WINDOW_WIDTH;
TRANFER_FUNCTION = @(OM) [sin(WINDOW_WIDTH*OM/2)./(OM*WINDOW_WIDTH.*(1 - (WINDOW_WIDTH*OM/[2*pi]).^2))].^2;

% White noise:
TT = [0:1:[(900*24)-1]]';
DD = randn(size(TT));

% Convolution:
DD_conv = conv(DD,WINDOW,'same');

% Plot:
close all
figure('Color','w')
plot(TT,DD,'.-'); hold on
plot(TT,DD_conv,'.-')

% Spectra:
figure('Color','w')
[SS     , f_vec, ~] = nanspectrum(DD     ,1,'hr',5,'.-',true,0); hold on
[SS_conv,   ~  , ~] = nanspectrum(DD_conv,1,'hr',5,'.-',true,0);
% subsampled:
dt_sample = 3;
[SS_ss ,      f_vec_ss, ~] = nanspectrum(DD(     1:dt_sample:end),dt_sample,'hr',5,'.-',true,0); hold on
[SS_conv_ss ,    ~    , ~] = nanspectrum(DD_conv(1:dt_sample:end),dt_sample,'hr',5,'.-',true,0);

% Transfer function:
figure('Color','w')
semilogy(f_vec, SS_conv./SS, 'k-'); hold on
semilogy(f_vec_ss, SS_conv_ss./SS_ss, 'g-');
semilogy(f_vec,TRANFER_FUNCTION(f_vec*2*pi),'r')
