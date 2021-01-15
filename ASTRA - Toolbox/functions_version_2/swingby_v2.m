function v2 = swingby_v2(v1, rp, mu, gam, n_r)

% rp is a vector

mu_rp     = mu./rp;
theta_inf = acos(-mu_rp./(norm(v1)^2 + mu_rp));
delta     = 2.*theta_inf-pi;                  

n_pi = eulerAxisAngle_v2(n_r, v1, gam);
v2   = eulerAxisAngle_v2(v1, n_pi, delta);

return