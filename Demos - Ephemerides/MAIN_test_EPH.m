%% Earth June 2005 - September 2013

clearDeleteAll;

mu = 132724487690;
AU = 149597870.7;

pl = 3;

% compute with horizon
horizon = load('horizons_EARTH.txt');

tjd = horizon(:,1);

posHOR = horizon(:, 3:5);
velHOR = horizon(:, 6:8);

tmjd2000 = zeros(length(tjd), 1);
for indi = 1:length(tjd)
    tmjd2000(indi,1) = jd2mjd2000(tjd(indi));
end

% compute with UPLANET and AA
posUPL = zeros(length(tjd), 3);
velUPL = zeros(length(tjd), 3);
posAAE = zeros(length(tjd), 3);
velAAE = zeros(length(tjd), 3);
posCCE = zeros(length(tjd), 3);
velCCE = zeros(length(tjd), 3);
for indi = 1:length(tmjd2000)
    [posUPL(indi,:), velUPL(indi,:)] = EphSS_car(pl, tmjd2000(indi));
    [posAAE(indi,:), velAAE(indi,:)] = EphAA_car(pl, tmjd2000(indi));
    [posCCE(indi,:), velCCE(indi,:)] = ephj2000(pl, tmjd2000(indi));
end

diffPOS_UPL = abs(posHOR - posUPL);
diffVEL_UPL = abs(velHOR - velUPL); 

diffPOS_AAE = abs(posHOR - posAAE);
diffVEL_AAE = abs(velHOR - velAAE); 

diffPOS_CCE = abs(posHOR - posCCE);
diffVEL_CCE = abs(velHOR - velCCE); 


figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\E_diffPOS_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\E_diffVEL_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\E_diffPOS_AAE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\E_diffVEL_AAE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_CCE(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_CCE(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_CCE(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\E_diffPOS_CCE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_CCE(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_CCE(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_CCE(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\E_diffVEL_CCE.jpg']);

%% Venus June 2005 - September 2013

pl = 2;

% compute with horizon
horizon = load('horizons_VENUS.txt');
horizon = horizon(:,~all(horizon == 0, 1));

tjd = horizon(:,1);

posHOR = horizon(:, 3:5);
velHOR = horizon(:, 6:8);

tmjd2000 = zeros(length(tjd), 1);
for indi = 1:length(tjd)
    tmjd2000(indi,1) = jd2mjd2000(tjd(indi));
end

% compute with UPLANET and AA
posUPL = zeros(length(tjd), 3);
velUPL = zeros(length(tjd), 3);
posAAE = zeros(length(tjd), 3);
velAAE = zeros(length(tjd), 3);
posCCE = zeros(length(tjd), 3);
velCCE = zeros(length(tjd), 3);
for indi = 1:length(tmjd2000)
    [posUPL(indi,:), velUPL(indi,:)] = EphSS_car(pl, tmjd2000(indi));
    [posAAE(indi,:), velAAE(indi,:)] = EphAA_car(pl, tmjd2000(indi));
    [posCCE(indi,:), velCCE(indi,:)] = ephj2000(pl, tmjd2000(indi));
end

diffPOS_UPL = abs(posHOR - posUPL);
diffVEL_UPL = abs(velHOR - velUPL); 

diffPOS_AAE = abs(posHOR - posAAE);
diffVEL_AAE = abs(velHOR - velAAE); 

diffPOS_CCE = abs(posHOR - posCCE);
diffVEL_CCE = abs(velHOR - velCCE); 

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\V_diffPOS_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\V_diffVEL_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\V_diffPOS_AAE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\V_diffVEL_AAE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_CCE(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_CCE(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_CCE(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\V_diffPOS_CCE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_CCE(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_CCE(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_CCE(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\V_diffVEL_CCE.jpg']);

%% Mars June 2005 - September 2013

pl = 4;

% compute with horizon
horizon = load('horizons_MARS.txt');
horizon = horizon(:,~all(horizon == 0, 1));

tjd = horizon(:,1);

posHOR = horizon(:, 3:5);
velHOR = horizon(:, 6:8);

tmjd2000 = zeros(length(tjd), 1);
for indi = 1:length(tjd)
    tmjd2000(indi,1) = jd2mjd2000(tjd(indi));
end

% compute with UPLANET and AA
posUPL = zeros(length(tjd), 3);
velUPL = zeros(length(tjd), 3);
posAAE = zeros(length(tjd), 3);
velAAE = zeros(length(tjd), 3);
posCCE = zeros(length(tjd), 3);
velCCE = zeros(length(tjd), 3);
for indi = 1:length(tmjd2000)
    [posUPL(indi,:), velUPL(indi,:)] = EphSS_car(pl, tmjd2000(indi));
    [posAAE(indi,:), velAAE(indi,:)] = EphAA_car(pl, tmjd2000(indi));
    [posCCE(indi,:), velCCE(indi,:)] = ephj2000(pl, tmjd2000(indi));
end

diffPOS_UPL = abs(posHOR - posUPL);
diffVEL_UPL = abs(velHOR - velUPL);

diffPOS_AAE = abs(posHOR - posAAE);
diffVEL_AAE = abs(velHOR - velAAE);

diffPOS_CCE = abs(posHOR - posCCE);
diffVEL_CCE = abs(velHOR - velCCE);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\M_diffPOS_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\M_diffVEL_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\M_diffPOS_AAE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\M_diffVEL_AAE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_CCE(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_CCE(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_CCE(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\M_diffPOS_CCE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_CCE(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_CCE(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_CCE(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\M_diffVEL_CCE.jpg']);

%% Mercury June 2005 - September 2013

pl = 1;

% compute with horizon
horizon = load('horizons_MERCURY.txt');
horizon = horizon(:,~all(horizon == 0, 1));

tjd = horizon(:,1);

posHOR = horizon(:, 3:5);
velHOR = horizon(:, 6:8);

tmjd2000 = zeros(length(tjd), 1);
for indi = 1:length(tjd)
    tmjd2000(indi,1) = jd2mjd2000(tjd(indi));
end

% compute with UPLANET and AA
posUPL = zeros(length(tjd), 3);
velUPL = zeros(length(tjd), 3);
posAAE = zeros(length(tjd), 3);
velAAE = zeros(length(tjd), 3);
for indi = 1:length(tmjd2000)
    [posUPL(indi,:), velUPL(indi,:)] = EphSS_car(pl, tmjd2000(indi));
    [posAAE(indi,:), velAAE(indi,:)] = EphAA_car(pl, tmjd2000(indi));
end

diffPOS_UPL = abs(posHOR - posUPL);
diffVEL_UPL = abs(velHOR - velUPL); 

diffPOS_AAE = abs(posHOR - posAAE);
diffVEL_AAE = abs(velHOR - velAAE); 

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\Y_diffPOS_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\Y_diffVEL_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\Y_diffPOS_AAE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\Y_diffVEL_AAE.jpg']);

%% Jupiter June 2005 - September 2013

pl = 5;

% compute with horizon
horizon = load('horizons_JUPITER.txt');
horizon = horizon(:,~all(horizon == 0, 1));

tjd = horizon(:,1);

posHOR = horizon(:, 3:5);
velHOR = horizon(:, 6:8);

tmjd2000 = zeros(length(tjd), 1);
for indi = 1:length(tjd)
    tmjd2000(indi,1) = jd2mjd2000(tjd(indi));
end

% compute with UPLANET and AA
posUPL = zeros(length(tjd), 3);
velUPL = zeros(length(tjd), 3);
posAAE = zeros(length(tjd), 3);
velAAE = zeros(length(tjd), 3);
for indi = 1:length(tmjd2000)
    [posUPL(indi,:), velUPL(indi,:)] = EphSS_car(pl, tmjd2000(indi));
    [posAAE(indi,:), velAAE(indi,:)] = EphAA_car(pl, tmjd2000(indi));
end

diffPOS_UPL = abs(posHOR - posUPL);
diffVEL_UPL = abs(velHOR - velUPL); 

diffPOS_AAE = abs(posHOR - posAAE);
diffVEL_AAE = abs(velHOR - velAAE); 

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\J_diffPOS_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\J_diffVEL_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\J_diffPOS_AAE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\J_diffVEL_AAE.jpg']);

%% Saturn June 2005 - September 2013

pl = 6;

% compute with horizon
horizon = load('horizons_SATURN.txt');
horizon = horizon(:,~all(horizon == 0, 1));

tjd = horizon(:,1);

posHOR = horizon(:, 3:5);
velHOR = horizon(:, 6:8);

tmjd2000 = zeros(length(tjd), 1);
for indi = 1:length(tjd)
    tmjd2000(indi,1) = jd2mjd2000(tjd(indi));
end

% compute with UPLANET and AA
posUPL = zeros(length(tjd), 3);
velUPL = zeros(length(tjd), 3);
posAAE = zeros(length(tjd), 3);
velAAE = zeros(length(tjd), 3);
for indi = 1:length(tmjd2000)
    [posUPL(indi,:), velUPL(indi,:)] = EphSS_car(pl, tmjd2000(indi));
    [posAAE(indi,:), velAAE(indi,:)] = EphAA_car(pl, tmjd2000(indi));
end

diffPOS_UPL = abs(posHOR - posUPL);
diffVEL_UPL = abs(velHOR - velUPL); 

diffPOS_AAE = abs(posHOR - posAAE);
diffVEL_AAE = abs(velHOR - velAAE); 

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_UPL(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\S_diffPOS_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_UPL(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_UPL(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\S_diffVEL_UPL.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km');
plot(tmjd2000, diffPOS_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffPOS_AAE(:,3), 'linewidth', 2);
legend('x', 'y', 'z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\S_diffPOS_AAE.jpg']);

figure;
hold on; grid on;
xlabel('Epoch - MJD2000'); ylabel('Difference - km/s');
plot(tmjd2000, diffVEL_AAE(:,1), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,2), 'linewidth', 2);
plot(tmjd2000, diffVEL_AAE(:,3), 'linewidth', 2);
legend('v_x', 'v_y', 'v_z');
legend('Location', 'northwest');
% saveas(gcf, [pwd '\images_EPH\S_diffVEL_AAE.jpg']);
