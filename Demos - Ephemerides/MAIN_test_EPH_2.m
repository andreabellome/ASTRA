%% Mars flight-path angle effect

clearDeleteAll;

mu = 132724487690;
AU = 149597870.7;

pl = 4;

% compute with Horizon
horizon = load('horizons_MARS.txt');
horizon = horizon(:,~all(horizon == 0, 1));

indi               = 1;                                  % Mars pericentre
indisum            = 25;
kepHOR(indi,:)     = car2kep(horizon(indi+indisum, 3:8), mu); % keplerian elements
MARS_orbit(indi,:) = horizon(indi+indisum, :);

% extract Mars orbit over one period
while abs(rad2deg(kepHOR(indi,end)) - 360) > 0.5
    indi               = indi + 1;
    kepHOR(indi,:)     = car2kep(horizon(indi+indisum, 3:8), mu);
    MARS_orbit(indi,:) = horizon(indi+indisum, :);
end

posHOR = MARS_orbit(:, 3:5);
velHOR = MARS_orbit(:, 6:8);
vzHOR  = velHOR(:,3);

tmjd2000 = zeros(size(MARS_orbit, 1), 1);
for indi = 1:size(MARS_orbit, 1)
    tmjd2000(indi,1) = jd2mjd2000(MARS_orbit(indi,1));
end

% compute CE
posCE = zeros(size(MARS_orbit, 1), 3);
velCE = zeros(size(MARS_orbit, 1), 3);
kepCE = zeros(size(MARS_orbit, 1), 6);
for indi = 1:size(MARS_orbit,1)
    [posCE(indi,:), velCE(indi,:)] = ephj2000(pl, tmjd2000(indi));
    kepCE(indi,:) = car2kep([posCE(indi,:), velCE(indi,:)], mu);
end
vzCE = velCE(:,3);

%%

% compute differences between CE an HE in vr and vth
vrHOR  = zeros(size(MARS_orbit,1), 1);
vthHOR = zeros(size(MARS_orbit,1), 1);
kHOR   = zeros(size(MARS_orbit,1), 1);
vrCE   = zeros(size(MARS_orbit,1), 1);
vthCE  = zeros(size(MARS_orbit,1), 1);
kCE    = zeros(size(MARS_orbit,1), 1);

for indi = 1:size(MARS_orbit,1)
    
    % HE in perifocal frame
    [~, ~, vrHOR(indi,1), vthHOR(indi,1), kHOR(indi)] = ...
        ECI2perifocal(posHOR(indi,:), velHOR(indi,:), mu);
     
    % CE in perifocal frame
    [~, ~, vrCE(indi,1), vthCE(indi,1), kCE(indi)] = ...
       ECI2perifocal(posCE(indi,:), velCE(indi,:), mu);
end

diffVR  = abs(vrHOR - vrCE);
diffVTH = abs(vthHOR - vthCE);
diffVZ  = abs(vzHOR - vzCE);

THS = linspace(0, 2*pi, length(tmjd2000));

figure;
hold on; grid on;
xlabel('True anomaly - deg'); ylabel('Difference - km/s');
plot(rad2deg(THS), diffVR, 'linewidth', 2);
plot(rad2deg(THS), diffVTH, 'linewidth', 2);
plot(rad2deg(THS), diffVZ, 'linewidth', 2);
legend('Radial', 'Transverse', 'Out of plane');
ylim([0 3])
xlim([0 360])
vline(0,   'k--', 'Pericentre')
vline(180, 'k--', 'Apocentre')

% %%
% 
% figure;
% hold on; grid on;
% xlabel('Epoch - MJD2000'); ylabel('Velocity - km/s');
% plot(tmjd2000, vrHOR, 'linewidth', 2);
% plot(tmjd2000, vthHOR, 'linewidth', 2);
% legend('v_r', 'v_{\theta}');
% 
% figure;
% hold on; grid on;
% xlabel('Epoch - MJD2000'); ylabel('Velocity - km/s');
% plot(tmjd2000, vrCE, 'linewidth', 2);
% plot(tmjd2000, vthCE, 'linewidth', 2);
% legend('v_r', 'v_{\theta}');
