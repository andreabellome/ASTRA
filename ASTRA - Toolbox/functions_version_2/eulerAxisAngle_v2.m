function v1 = eulerAxisAngle_v2(v, n, theta)

v = v(:);
n = n/norm(n);

R11 = cos(theta)+(1-cos(theta)).*n(1)^2;
R12 = (1-cos(theta)).*n(1)*n(2)+sin(theta).*n(3);
R13 = (1-cos(theta)).*n(1)*n(3)-sin(theta).*n(2);
R21 = (1-cos(theta)).*n(1)*n(2)-sin(theta).*n(3);
R22 = cos(theta)+(1-cos(theta)).*n(2)^2;
R23 = (1-cos(theta)).*n(2)*n(3)+sin(theta).*n(1);
R31 = (1-cos(theta)).*n(1)*n(3)+sin(theta).*n(2);
R32 = (1-cos(theta)).*n(2)*n(3)-sin(theta).*n(1);
R33 = cos(theta)+(1-cos(theta)).*n(3)^2;

v1(:,1) = R11.*v(1) + R21.*v(2) + R31.*v(3);
v1(:,2) = R12.*v(1) + R22.*v(2) + R32.*v(3);
v1(:,3) = R13.*v(1) + R23.*v(2) + R33.*v(3);

return