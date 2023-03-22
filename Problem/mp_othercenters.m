function othercenters = mp_othercenters(prob, xu)
np = 5; % number of peaks

xl_center= prob.get_xlprime(xu);

nx    = size(xu, 1);
w     = ones(np, 1) * 0.02; % width 0.001 for 3d
h     = ones(np, 1) * 50; % height
h(1)  = 70;

v = 0.3 .*  (prob.xl_bu(prob.q+1 : end) - prob.xl_bl(prob.q+1 : end));      % distances between basin centers
v = [zeros(1, prob.q), v];

% alpha parameter, make basin of attractions  round
delta = 80 ./ (prob.xl_bu -prob.xl_bl);

forward  =  xl_center + v;
backward = xl_center - v;

v2 = 0.3 .*  (prob.xl_bu(1 : prob.q) - prob.xl_bl(1 : prob.q));
v2 = [v2, zeros(1, prob.r)];
left  = xl_center + v2;
right =  xl_center - v2;
othercenters = [forward; backward; left; right];

end