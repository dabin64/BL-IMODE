function [functionValue equalityConstrVals inequalityConstrVals]=llTestProblem(llPop, testProblemName, ulMember)

    %This function evaluates the lower level objective values and constraints
    %for a set of lower level members corresponding to a given upper level member.
    global llFunctionEvaluations
    global   llFE 
    noOfMembers = size(llPop,1);
    llFunctionEvaluations = llFunctionEvaluations + noOfMembers;
    llFE                  = llFE  + noOfMembers;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function call here
    fhandle = str2func(testProblemName);
    equalityConstrVals = [];
    inequalityConstrVals = [];
    
    for i=1:noOfMembers
        if size(llPop,1) == size(ulMember,1)
            [functionValue(i,:) equalityConstrValsTemp inequalityConstrValsTemp] = fhandle(ulMember(i,:), llPop(i,:));
        elseif size(ulMember,1) == 1
            [functionValue(i,:) equalityConstrValsTemp inequalityConstrValsTemp] = fhandle(ulMember, llPop(i,:));
        else
            disp('Error in llTestProblem, size of llPop and ulMember mismatch');
        end
        
        if ~isempty(equalityConstrValsTemp)
            equalityConstrVals(i,:) = equalityConstrValsTemp;
        end
        if ~isempty(inequalityConstrValsTemp)
            inequalityConstrVals(i,:) = inequalityConstrValsTemp;
        end
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd1(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - tan(xl2)).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd2(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - log(xl2)).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd3(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                    + q + sum(xl1.^2 - cos(2*pi*xl1)) ...
                    + sum((xu2.^2 - tan(xl2)).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd4(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                        + q + sum(xl1.^2 - cos(2*pi*xl1)) ...
                        + sum((abs(xu2) - log (1+xl2)).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd5(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    term2 = 0;
    for i=1:q-1
        term2 = term2 + (xl1(i+1) - xl1(i).^2).^2 + (xl1(i) - 1).^2;
    end
    
    functionValue = sum((xu1).^2) ...
                        + term2 ...
                        + sum((abs(xu2) - xl2.^2).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd6(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = floor((length(xl) - r)/2 - eps);
    s = ceil((length(xl) - r)/2 + eps);
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q+s);
    xl2 = xl(q+s+1:q+s+r);

    term2 = sum(xl1(1:q).^2);
    for i=q+1:2:q+s-1
        term2 = term2 + (xl1(i+1) - xl1(i)).^2;
    end
    
    functionValue = sum((xu1).^2) ...
                    + term2 ...
                    + sum((xu2 - xl2).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd7(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^3) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - log(xl2)).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd8(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    term2 = 0;
    for i=1:q-1
        term2 = term2 + (xl1(i+1) - xl1(i).^2).^2 + (xl1(i) - 1).^2;
    end
    
    functionValue = sum(abs(xu1)) ...
                        + term2 ...
                        + sum((xu2 - xl2.^3).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd9(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - log(1+xl2)).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals(1) = sum(xl1.^2)+sum(xl2.^2) - floor(sum(xl1.^2)+sum(xl2.^2)+0.5);
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd10(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    b = 2*ones(size(xl1));

    functionValue = sum((xu1).^2) ...
                    + sum((xl1 - b).^2) ...
                    + sum((xu2 - tan(xl2)).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    for i=1:q
        inequalityConstrVals(i) = xl1(i) + xl1(i).^3 - sum(xl1.^3);
    end
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd11(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - log(xl2)).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals(1) = sum((xu2 - log(xl2)).^2) - 1;
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd12(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    b = 2*ones(size(xl1));

    functionValue = sum((xu1).^2) ...
                    + sum((xl1 - b).^2) ...
                    + sum((xu2 - tan(xl2)).^2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    for i=1:q
        inequalityConstrVals(i) = xl1(i) + xl1(i).^3 - sum(xl1.^3);
    end
    inequalityConstrVals(q+1) = sum((xu2 - tan(xl2)).^2) - 1;
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd13(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);
    
    term1 = sum(abs(xu1) + 2*abs(sin(xu1)));
    term2 = 0;
    for i=1:q
        term2 = term2 + sum(xl1(1:i).^2);
    end
    term3 = sum((xu2 - log(xl2)).^2);

    functionValue = term1 + term2 + term3;

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd14(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = floor((length(xl) - r)/2 - eps);
    s = ceil((length(xl) - r)/2 + eps);
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q+s);
    xl2 = xl(q+s+1:q+s+r);

    term1 = sum(floor(xu1));
    term2 = 0;
    for i=1:q
        term2 = term2 + abs(xl1(i))^(i+1);
    end
    for i=q+1:2:q+s-1
        term2 = term2 + (xl1(i+1) - xl1(i)).^2;
    end
    term3 = sum(abs(xu2.^2 - xl2.^2));
    
    functionValue = term1 + term2 + term3;

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = tp1(xu,xl)

    functionValue = (xu(1)- xl(1)).^2 + (xu(2) - xl(2)).^2;

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals = [];
%      inequalityConstrVals(1) = -xl(1);
%      inequalityConstrVals(2) = -xl(2);
%      inequalityConstrVals(3) = xl(1)-10;
%      inequalityConstrVals(4) = xl(2)-10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp2(xu,xl)

    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);
    functionValue = (y1-x1+20).^2 + (y2-x2+20).^2;

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = 10-x1+2*y1;
    inequalityConstrVals(2) = 10-x2+2*y2;
%     inequalityConstrVals(3) = -20+y1;
%     inequalityConstrVals(4) = -20+y2;
%     inequalityConstrVals(5) = -10-y1;
%     inequalityConstrVals(6) = -10-y2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp3(xu,xl)
    
    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);
    functionValue = 2*x1.^2+y1.^2-5*y2;

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = -3-x1.^2+2*x1-x2.^2+2*y1-y2;
    inequalityConstrVals(2) = 4-x2-3*y1+4*y2;
%     inequalityConstrVals(3) = -y1;
%     inequalityConstrVals(4) = -y2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp4(xu,xl)
    
    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);
    y3 = xl(3);
    functionValue = x1+2*x2+y1+y2+2*y3;

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = y2+y3-y1-1;
    inequalityConstrVals(2) = 2*x1-y1+2*y2-0.5*y3-1;
    inequalityConstrVals(3) = 2*x2+2*y1-y2-0.5*y3-1;
%     inequalityConstrVals(4) = -y1;
%     inequalityConstrVals(5) = -y2;
%     inequalityConstrVals(6) = -y3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp5(xu,xl)
    
    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);
    
    H = [1 3; 3 10];
    b = [-1 2; 3 -3];
    x = [x1 x2]';
    y = [y1 y2]';
    
    functionValue = 0.5*y'*H*y + (b*x)'*y;

    functionValue = functionValue';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = -0.333*y1 + y2 - 2;
    inequalityConstrVals(2) = y1 - 0.333*y2 -2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp6(xu,xl)
    
    x1 = xu(1);
    y1 = xl(1);
    y2 = xl(2);
    
    functionValue = (2*y1-4).^2 + (2*y2-1).^2 + x1*y1;

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = 4*x1+5*y1+4*y2-12;
    inequalityConstrVals(2) = 4*y2-4*x1-5*y1+4;
    inequalityConstrVals(3) = 4*x1-4*y1+5*y2-4;
    inequalityConstrVals(4) = 4*y1-4*x1+5*y2-4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp7(xu,xl)
    
    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);

    functionValue = (x1+y1).*(x2+y2)./(1+x1.*y1+x2.*y2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = y1-x1;
    inequalityConstrVals(2) = y2-x2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp8(xu,xl)
    
    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);

    functionValue = (y1-x1+20).^2+(y2-x2+20).^2;

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = 2*y1-x1+10;
    inequalityConstrVals(2) = 2*y2-x2+10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp9(xu,xl)
    
    llDim = length(xl);
    x = xu;
    y = xl;

    exponent = 1+1/4000*sum(y.^2)-prod(cos(y./sqrt(1:llDim)));
    exponent = exponent*sum(x.^2);
    functionValue = exp(exponent);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = tp10(xu,xl)
    
    llDim = length(xl);
    x = xu;
    y = xl;
    
    exponent = 1+1/4000*sum((y.*x).^2)-prod(cos((y.*x)./sqrt(1:llDim)));
    functionValue = exp(exponent);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = msmd1(xu, xl)

    %SMD1 modified to a multi-objective problem.

    ulMember = xu;
    pop = xl;
    r = floor(size(ulMember,2)/2);
    p = size(ulMember,2) - r;
    q = size(pop,2) - r;

    noOfMembers = size(pop,1);
    
    xu1 = ones(noOfMembers,1)*ulMember(1:p);
    xu2 = ones(noOfMembers,1)*ulMember(p+1:p+r);

    xl1 = pop(:,1:q);
    xl2 = pop(:,q+1:q+r);

    b = zeros(size(xl1));
    d = zeros(size(xu1));

    %Different from sdm1 (third line has tan of xl2 instead of xl2)
    functionValue = sum((xu1 - d).^2,2) ...
                    + sum((xl1 - b).^2,2) ...
                    + sum((xu2 - (xl2)).^2,2);

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = externalProblem(xu,xl)

    [functionValue equalityConstrVals inequalityConstrVals] = llExternalProblem(xu, xl);

function [functionValue equalityConstrVals inequalityConstrVals] = GoldMining(xu,xl)
    
    % beta_ = 0.5;
    tau = xu(1);
    alpha_ = xu(2);
    q = xl(1);

    f1 = -((100-q)*q - (q^2+q) - tau*q);
    f2 = q;

% --WS--
    % functionValue = alpha_ * f1 + (1-alpha_) * f2;

    z_utopia1 = -1300;
    z_utopia2 = -1;
%     assert(z_utopia1 < f1);
%     assert(z_utopia2 < f2);

% --TCH--
%     functionValue = max(alpha_ * (f1-z_utopia1),(1-alpha_) * (f2-z_utopia2));

% --PBI--
    d1 = abs((f1-z_utopia1)*alpha_ + (f2-z_utopia2)*(1-alpha_))/sqrt(alpha_^2+(1-alpha_)^2);
    d2 = sqrt((f1-z_utopia1-d1*alpha_)^2 + (f2-z_utopia2-d1*(1-alpha_))^2);
    functionValue = d1 + 5*d2;

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = f1;

function [functionValue equalityConstrVals inequalityConstrVals] = DecisionMaking(xu,xl)
    
    % beta_ = 0.5;
    x1 = xu(1); x2 = xu(2);
    alpha_ = xu(3);
    y1 = xl(1); y2 = xl(2); y3 = xl(3);

    f1 = -(4.47*x1 + 5.46*x2 - 6.23*y1 - 4.78*y2 + 7.34*y3);
    f2 = -(5.34*x1 + 3.74*y2 + 9.45*y1 + 6.37*y2 + 5.45*y3);

% --WS--
    % functionValue = alpha_ * f1 + (1-alpha_) * f2;

    z_utopia1 = -3000;
    z_utopia2 = -3000;
%     assert(z_utopia1 < f1);
%     assert(z_utopia2 < f2);

% --TCH--
%     functionValue = max(alpha_ * (f1-z_utopia1),(1-alpha_) * (f2-z_utopia2));

% --PBI--
    d1 = abs((f1-z_utopia1)*alpha_ + (f2-z_utopia2)*(1-alpha_))/sqrt(alpha_^2+(1-alpha_)^2);
    d2 = sqrt((f1-z_utopia1-d1*alpha_)^2 + (f2-z_utopia2-d1*(1-alpha_))^2);
    functionValue = d1 + 5*d2;

    functionValue = functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];

    inequalityConstrVals(1) = 3.67*x1 - 7.84*x2 - 6.78*y1 - 5.87*y2 + 1.26*y3 - 105;
    inequalityConstrVals(2) = 4.34*x1 + 9.26*x2 + 8.33*y1 - 2.11*y2 - 2.67*y3 - 830;
    inequalityConstrVals(3) = 4.51*x1 - 3.56*x2 + 1.46*y1 + 2.42*y2 + 7.43*y3 - 565;
    inequalityConstrVals(4) = -(4.47*x1 + 5.46*x2 - 6.23*y1 - 4.78*y2 + 7.34*y3);
    inequalityConstrVals(5) = -(5.34*x1 + 3.74*x2 + 9.45*y1 + 6.37*y2 + 5.45*y3);
    function [functionValue equalityConstrVals inequalityConstrVals] = smd1mp(xu, xl)
    prob1         = smd1mp1(1,2,1);
    [fl, flc] = prob1.evaluate_l(xu, xl);
    functionValue = fl;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [functionValue equalityConstrVals inequalityConstrVals] = smd2mp(xu, xl)
    prob1         = smd2mp1(1,2,1);
    [fl, flc] = prob1.evaluate_l(xu, xl);
    functionValue = fl;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [functionValue equalityConstrVals inequalityConstrVals] = smd3mp(xu, xl)
    prob1         = smd3mp1(1,2,1);
    [fl, flc] = prob1.evaluate_l(xu, xl);
    functionValue = fl;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [functionValue equalityConstrVals inequalityConstrVals] = smd4mp(xu, xl)
    prob1         = smd4mp1(1,2,1);
    [fl, flc] = prob1.evaluate_l(xu, xl);
    functionValue = fl;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [functionValue equalityConstrVals inequalityConstrVals] = smd5mp(xu, xl)
    prob1         = smd5mp1(1,2,1);
    [fl, flc] = prob1.evaluate_l(xu, xl);
    functionValue = fl;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [functionValue equalityConstrVals inequalityConstrVals] = smd6mp(xu, xl)
    prob1         = smd6mp1(1,2,1);
    [fl, flc] = prob1.evaluate_l(xu, xl);
    functionValue = fl;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [functionValue equalityConstrVals inequalityConstrVals] = smd7mp(xu, xl)
    prob1         = smd7mp1(1,2,1);
    [fl, flc] = prob1.evaluate_l(xu, xl);
    functionValue = fl;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [functionValue equalityConstrVals inequalityConstrVals] = smd8mp(xu, xl)
    prob1         = smd8mp1(1,2,1);
    [fl, flc] = prob1.evaluate_l(xu, xl);
    functionValue = fl;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
