function fmp = mp_module_smd6only(prob, xu, xl )
% This function calcluate moving peaks on the lower level 
% for given xu xl
% returns its moving peak f value.
% moving peaks are 5 in total, the center one co-exists with
% orginal lower level function's global optimum
%--------------------------------------------------------

np = 5; % number of peaks
fmp = [];
xl_center= prob.get_xlprime(xu);


nx    = size(xu, 1);
w     = ones(np, 1) * 0.01; % width
h      = ones(np, 1) * 50; % height
h(1)  = 70;

v = 0.3 .*  (prob.xl_bu(prob.q + prob.s+1 : end) - prob.xl_bl(prob.q + prob.s+1 : end));      % distances between basin centers
v = [zeros(1, prob.q + prob.s), v];

% alpha parameter, make basin of attractions  round
delta = 80 ./ (prob.xl_bu -prob.xl_bl);

forward =  xl_center + v;
backward = xl_center - v;

v2 = 0.3 .*  (prob.xl_bu(1 : prob.q + prob.s) - prob.xl_bl(1 : prob.q + prob.s));
v2 = [v2, zeros(1, prob.r)];
left = xl_center + v2;
right =  xl_center - v2;
mp_centers = {xl_center, forward, backward, left, right};

fmp = [];
for j = 1:length(mp_centers)
    fj = h(j) ./ (1 + w(j)  .* sum(((xl - mp_centers{j}) .*  delta).^2,  2));
    fmp = [fmp, fj];
end
fmp = max(fmp, [], 2);
end