clear all; close all; clc; format long g;

date     = [2020 1 1 0 0 0];   % date 
tmjd2000 = date2mjd2000(date); % from date to mjd2000

date0 = [2000 1 1 12 0 0]; % intial date
date1 = [2050 1 1 12 0 0]; % final date
tdate = 1; % date step - days

t0mjd2000 = date2mjd2000(date0); % initial epoch mjd2000
t1mjd2000 = date2mjd2000(date1); % final epoch mjd2000
tt        = t0mjd2000:tdate:t1mjd2000; % grid of epochs

planets = [1 2 3 4 5 6 7 8]; % planet - Earth

%% compare uplanet with ephj2000

for indpl = 1:length(planets)
    
    pl = planets(indpl);

    rr1  = zeros(length(tt), 3);
    vv1  = zeros(length(tt), 3);
    rr2  = zeros(length(tt), 3);
    vv2  = zeros(length(tt), 3);
    DRR  = zeros(length(tt), 3);
    DR3D = zeros(length(tt), 1);
    DR2D = zeros(length(tt), 1);
    for indi = 1:length(tt)

        % ephemerides using uplanet
        [rr1(indi,:), vv1(indi,:)] = EphSS_car(pl, tt(indi));

        % ephemerides using circular coplanar model
        [rr2(indi,:), vv2(indi,:)] = ephj2000(pl, tt(indi));

        % difference between the two
        DRR(indi,:)  = abs(rr1(indi,:) - rr2(indi,:));
        DR3D(indi,1) = norm(DRR(indi,:));   % consider also z axis 
        DR2D(indi,1) = norm(DRR(indi,1:2)); % only x and y axis

    end

    % save in structure
    Struc1(:,pl).Planet               = pl;
    Struc1(:,pl).StateVector_uplanet  = [rr1 vv1];
    Struc1(:,pl).StateVector_ephj2000 = [rr2 vv2];
    Struc1(:,pl).Difference           = DRR;
    Struc1(:,pl).Difference_3Dnorm    = DR3D;
    Struc1(:,pl).Difference_2Dnorm    = DR2D;

    figure;
    plot(tt, DR2D./1e6);
    xlabel('Epoch - MJD2000'); ylabel('2D Difference - millions of km');
    title(planetIdToName(pl));

end

%% compare uplanet with de405

for indpl = 1:length(planets)
    
    pl = planets(indpl);

    rr1  = zeros(length(tt), 3);
    vv1  = zeros(length(tt), 3);
    rr2  = zeros(length(tt), 3);
    vv2  = zeros(length(tt), 3);
    DRR  = zeros(length(tt), 3);
    DR3D = zeros(length(tt), 1);
    DR2D = zeros(length(tt), 1);
    for indi = 1:length(tt)

        % ephemerides using uplanet
        [rr1(indi,:), vv1(indi,:)] = EphSS_car(pl, tt(indi));

        % ephemerides using JPL 405 model
        [rr2(indi,:), vv2(indi,:)] = planetEphemeris(mjd20002jd(tt(indi)), 'Sun', planetIdToName(pl), '405');

        % difference between the two
        DRR(indi,:)  = abs(rr1(indi,:) - rr2(indi,:));
        DR3D(indi,1) = norm(DRR(indi,:));   % consider also z axis 
        DR2D(indi,1) = norm(DRR(indi,1:2)); % only x and y axis

    end

    % save in structure
    Struc2(:,pl).Planet               = pl;
    Struc2(:,pl).StateVector_uplanet  = [rr1 vv1];
    Struc2(:,pl).StateVector_ephj2000 = [rr2 vv2];
    Struc2(:,pl).Difference           = DRR;
    Struc2(:,pl).Difference_3Dnorm    = DR3D;
    Struc2(:,pl).Difference_2Dnorm    = DR2D;

    figure;
    plot(tt, DR2D./1e6);
    xlabel('Epoch - MJD2000'); ylabel('2D Difference - millions of km');
    title(planetIdToName(pl));

end

