function [functionValue equalityConstrVals inequalityConstrVals]=ulTestProblem(ulPop, llPop, testProblemName)

    %This function evaluates the upper level objective values and constraints
    %for a set of upper level members and their corresponding lower level members.
    global ulFunctionEvaluations;
    global UUFE;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function call here
    fhandle = str2func(testProblemName);
    noOfMembers = size(ulPop,1);
    ulFunctionEvaluations = ulFunctionEvaluations+noOfMembers;
    UUFE                  = UUFE+noOfMembers;

    equalityConstrVals = [];
    inequalityConstrVals = [];
    
    for i=1:noOfMembers
        [functionValue(i,:) equalityConstrValsTemp inequalityConstrValsTemp] = fhandle(ulPop(i,:), llPop(i,:));
        
        if ~isempty(equalityConstrValsTemp)
            equalityConstrVals(i,:) = equalityConstrValsTemp;
        end
        if ~isempty(inequalityConstrValsTemp)
            inequalityConstrVals(i,:) = inequalityConstrValsTemp;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd1(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);
    
    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2).^2) + sum((xu2 - tan(xl2)).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd2(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);   
    
    functionValue = sum((xu1).^2) ...
                    - sum((xl1).^2) ...
                    + sum((xu2).^2) - sum((xu2 - log(xl2)).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd3(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    c = zeros(size(xu2));    
    
    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - c).^2) + sum((xu2.^2 - tan(xl2)).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd4(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);
    
    functionValue = sum((xu1).^2) ...
                - sum((xl1).^2) ...
                + sum((xu2).^2) - sum((abs(xu2) - log(1+xl2)).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd5(xu, xl)

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
    
    %Same as smd5
    functionValue = sum((xu1).^2) ...
                - term2 ...
                + sum((xu2).^2) - sum((abs(xu2) - xl2.^2).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd6(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = floor((length(xl) - r)/2 - eps);
    s = ceil((length(xl) - r)/2 + eps);

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q+s);
    xl2 = xl(q+s+1:q+s+r);   
    
    functionValue = sum((xu1).^2) ...
                    - sum(xl1(1:q).^2) + sum(xl1(q+1:q+s).^2) ...
                    + sum((xu2).^2) - sum((xu2 - xl2).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd7(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);  
    
    m = [1:p];
    functionValue = 1+1/400*sum((xu1).^2) - prod(cos(xu1./sqrt(m))) ...
                    - sum((xl1).^2) ...
                    + sum((xu2).^2) - sum((xu2 - log(xl2)).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd8(xu, xl)

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
    
    functionValue = 20+exp(1)-20*exp(-0.2*sqrt(1/p*sum((xu1).^2))) - exp(1/p*sum(cos(2*pi*xu1)))  ...
                - term2 ...
                + sum((xu2).^2) - sum((xu2 - xl2.^3).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd9(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);  
    
    functionValue = sum((xu1).^2) ...
                    - sum((xl1).^2) ...
                    + sum((xu2).^2) - sum((xu2 - log(1+xl2)).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals(1) = sum(xu1.^2)+sum(xu2.^2) - floor(sum(xu1.^2)+sum(xu2.^2)+0.5);
    inequalityConstrVals(2) = sum(xl1.^2)+sum(xl2.^2) - floor(sum(xl1.^2)+sum(xl2.^2)+0.5);
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd10(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    a = 2*ones(size(xu1));
    c = 2*ones(size(xu2));    

    
    functionValue = sum((xu1 - a).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - c).^2) - sum((xu2 - tan(xl2)).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    for i=1:p
        inequalityConstrVals(i) = xu1(i) + xu1(i).^3 - sum(xu1.^3) - sum(xu2.^3);
    end

    for i=1:r
        inequalityConstrVals(p+i) = xu2(i) + xu2(i).^3 - sum(xu2.^3) - sum(xu1.^3);
    end

    inequalityConstrVals = -inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd11(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);
    
    functionValue = sum((xu1).^2) ...
                    - sum((xl1).^2) ...
                    + sum((xu2).^2) - sum((xu2 - log(xl2)).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = xu2 - 1/sqrt(r) - log(xl2);
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd12(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    a = 2*ones(size(xu1));
    c = 2*ones(size(xu2));    

    functionValue = sum((xu1 - a).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - c).^2) + sum(tan(abs(xl2))) - sum((xu2 - tan(xl2)).^2);

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    for i=1:p
        inequalityConstrVals(i) = xu1(i) + xu1(i).^3 - sum(xu1.^3) - sum(xu2.^3);
    end

    for i=1:r
        inequalityConstrVals(p+i) = xu2(i) + xu2(i).^3 - sum(xu2.^3) - sum(xu1.^3);
    end

    inequalityConstrVals(p+r+1:p+2*r) = xu2 - tan(xl2);
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd13(xu, xl)

    %Optimum at xu1=1, xu2=0, xl1=0, xl2=1, f=1+2sin(1), F=0
    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);
    
    term1 = (xu1(1) - 1).^2; %One term has been added here
    for i=1:p-1
        term1 = term1 + (xu1(i)-1)^2 + (xu1(i+1) - xu1(i).^2).^2;
    end
    term2 = 0;
    for i=1:q
        term2 = term2 - sum(xl1(1:i).^2);
    end
    term3 = -sum((xu2 - log(xl2)).^2);
    for i=1:r
        term3 = term3 + sum(xu2(1:i).^2);
    end
    
    functionValue = term1 + term2 + term3;

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd14(xu, xl)

    %Optimum at xu1=1, xu2=0, xl1=0, xl2=0, f=0, F=0
    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = floor((length(xl) - r)/2 - eps);
    s = ceil((length(xl) - r)/2 + eps);

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q+s);
    xl2 = xl(q+s+1:q+s+r);   
    
    term1 = (xu1(1) - 1).^2; %One term has been added here
    for i=1:p-1
        term1 = term1 + (xu1(i)-1)^2 + (xu1(i+1) - xu1(i).^2).^2;
    end
    term2 = 0;
    for i=1:q
        term2 = term2 - abs(xl1(i))^(i+1);
    end
    term2 = term2 + sum(xl1(q+1:q+s).^2);
    
    term3 = (sum(abs(xl2)));
    for i=1:r
        term3 = term3 + i*xu2(i).^2;
    end
    
    functionValue = term1 + term2 + term3;

    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = tp1(xu, xl)

    %At optima xu1=20,xu2=5,xl1=10,xl2=5,fu=-225
    functionValue = (xu(1)-30).^2 + (xu(2)-20).^2 - 20*xl(1) + 20*xl(2);
        
    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = 30 - xu(1) - 2*xu(2);
    inequalityConstrVals(2) = xu(1) + xu(2) - 25;
    inequalityConstrVals(3) = xu(2) - 15;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp2(xu, xl)

    %At optima x1=0,x2=30,y1=-10,xl2=10,fu=0,fl=-100
    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);
    functionValue = 2*x1+2*x2-3*y1-3*y2-60;
        
    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = x1+x2+y1-2*y2-40;
%     inequalityConstrVals(2) = x1-50;
%     inequalityConstrVals(3) = -x1;
%     inequalityConstrVals(4) = x2-50;
%     inequalityConstrVals(5) = -x2;
    %Lower level constraints included at upper level
    inequalityConstrVals(2) = 10-x1+2*y1;
    inequalityConstrVals(3) = 10-x2+2*y2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp3(xu, xl)

    %At optima x1=0,x2=2,y1=1.8750,xl2=0.9062,fu=18.6787,fl=1.0156
    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);
    functionValue = -x1.^2-3*x2.^2-4*y1+y2.^2;
        
    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = x1.^2+2*x2-4;
    %Lower level constraints included at upper level
    inequalityConstrVals(2) = -3-x1.^2+2*x1-x2.^2+2*y1-y2;
    inequalityConstrVals(3) = 4-x2-3*y1+4*y2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp4(xu, xl)

    %At optima x1=0.29,x2=0.70,y1=0,y2=0.27,y3=0.27,fu=29.2,fl=-3.2
    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);
    y3 = xl(3);
    functionValue = -8*x1-4*x2+4*y1-40*y2-4*y3;
        
    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    %inequalityConstrVals = [];
    %Lower level constraints included at upper level
    inequalityConstrVals(1) = y2+y3-y1-1;
    inequalityConstrVals(2) = 2*x1-y1+2*y2-0.5*y3-1;
    inequalityConstrVals(3) = 2*x2+2*y1-y2-0.5*y3-1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp5(xu, xl)

    %At optima x1=2,x2=0,y1=2,y2=0,fu=3.6,fl=2
    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);
        
    r = 0.1;
    x = [x1 x2]';
    y = [y1 y2]';
    
    functionValue = r*x'*x - 3*y1 - 4*y2 + 0.5*y'*y;
        
    functionValue = functionValue';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    %inequalityConstrVals = [];
    %Lower level constraints included at upper level
    inequalityConstrVals(1) = -0.333*y1 + y2 - 2;
    inequalityConstrVals(2) = y1 - 0.333*y2 -2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 function [functionValue equalityConstrVals inequalityConstrVals] = tp6(xu, xl)

    %At optima x1=1.888,y1=0.888,y2=0,fu=1.2098,fl=-7.61
    x1 = xu(1);
    y1 = xl(1);
    y2 = xl(2);

    functionValue = (x1-1).^2 + 2*y1 - 2*x1;
        
    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    %inequalityConstrVals = [];
    %Lower level constraints included at upper level
    inequalityConstrVals(1) = 4*x1+5*y1+4*y2-12;
    inequalityConstrVals(2) = 4*y2-4*x1-5*y1+4;
    inequalityConstrVals(3) = 4*x1-4*y1+5*y2-4;
    inequalityConstrVals(4) = 4*y1-4*x1+5*y2-4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp7(xu, xl)

    %At optima x1=5.24,x2=5.63,y1=5.24,y2=0,y3=0.27,fu=2.0714,fl=-2.0714
    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);

    functionValue = -(x1+y1).*(x2+y2)./(1+x1.*y1+x2.*y2);
        
    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = x1.^2+x2.^2 - 100;
    inequalityConstrVals(2) = x1-x2; %New constraint added at upper level to correct the bilevel problem
    %Lower level constraints included at upper level
    inequalityConstrVals(3) = y1-x1;
    inequalityConstrVals(4) = y2-x2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = tp8(xu, xl)

    %At optima x1=0,x2=30,y1=-10,y2=10,fu=0,fl=-100
    x1 = xu(1);
    x2 = xu(2);
    y1 = xl(1);
    y2 = xl(2);
    
    functionValue = abs(2*x1+2*x2-3*y1-3*y2-60);
        
    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals(1) = x1+x2+y1-2*y2-40;
    %Lower level constraints included at upper level
    inequalityConstrVals(2) = 2*y1-x1+10;
    inequalityConstrVals(3) = 2*y2-x2+10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = tp9(xu, xl)

    %At optima x1=0,x2=30,y1=-10,y2=10,fu=0,fl=-100
    x = xu(:); %x(i) = xu(i);
    y = xl(:); %y(i) = xl(i);
    
    functionValue = sum(abs(x-1))+sum(abs(y));
        
    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = tp10(xu, xl)

    %At optima x1=0,x2=30,y1=-10,y2=10,fu=0,fl=-100
    x = xu(:); %x(i) = xu(i);
    y = xl(:); %y(i) = xl(i);
    
    functionValue = sum(abs(x-1))+sum(abs(y));
        
    functionValue = functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    equalityConstrVals = [];
    inequalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = msmd1(xu, xl)

    %SMD1 modified to a multi-objective problem.

    ulPop = xu;
    llPop = xl;

    r = floor(size(ulPop,2)/2);
    p = size(ulPop,2) - r;
    q = size(llPop,2) - r;

    xu1 = ulPop(:,1:p);
    xu2 = ulPop(:,p+1:p+r);

    xl1 = llPop(:,1:q);
    xl2 = llPop(:,q+1:q+r);
    
    functionValue(:,1) = sum(xu1.^2,2) ...
                    + sum((xl1).^2,2) ...
                    + sum((xu2).^2,2) + sum((xu2 - (xl2)).^2,2);
    
    functionValue(:,2) = 1-abs(xu1(:,1)) + sum(xu1(:,2:end).^2,2) ...
                    + sum((xl1).^2,2) ...
                    + sum((xu2).^2,2) + sum((xu2 - (xl2)).^2,2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = externalProblem(xu,xl)

    [functionValue equalityConstrVals inequalityConstrVals] = ulExternalProblem(xu, xl);

function [functionValue equalityConstrVals inequalityConstrVals] = GoldMining(xu, xl)

    beta_ = 0.5;
    tau = xu(1);
    % alpha_ = xu(2);
    q = xl(1);

    F1 = -tau*q;
    F2 = q;
    functionValue = beta_ * F1 + (1-beta_) * F2;
    functionValue = functionValue;
    
    equalityConstrVals = [];
    inequalityConstrVals = [];

function [functionValue equalityConstrVals inequalityConstrVals] = DecisionMaking(xu, xl)

    beta_ = 0.5;
    x1 = xu(1); x2 = xu(2);
    % alpha_ = xu(3);
    y1 = xl(1); y2 = xl(2); y3 = xl(3);

    F1 = -(3.38*x1 + 7.78*x2 + 8.54*y1 - 2.35*y2 + 4.97*y3);
    F2 = -(6.64*x1 + 4.26*x2 + 4.67*y1 + 4.59*y2 + 3.73*y3); 

    functionValue = beta_ * F1 + (1-beta_) * F2;
    functionValue = functionValue;
    
    equalityConstrVals = [];
    inequalityConstrVals(1) = 4.55*x1 + 7.35*x2 + 9.65*y1 - 6.23*y2 + 4.24*y3 - 987;
    inequalityConstrVals(2) = -5.33*x1 - 1.35*x2 + 2.67*y1 - 4.22*y2 + 1.75*y3 - 135;
    inequalityConstrVals(3) = -2.11*x1 + 2.67*x2 + 4.34*y1 + 9.26*y2 + 8.33*y3 - 830;
    inequalityConstrVals(4) = 2.42*x1 + 7.43*x2 + 4.51*y1 - 3.56*y2 + 1.46*y3 - 565;
function [functionValue equalityConstrVals inequalityConstrVals] = smd1mp(xu, xl)
    prob         = smd1mp1(1, 2, 1);
    [fu, fuc]    = prob.evaluate_u(xu, xl);
    functionValue = fu;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [functionValue equalityConstrVals inequalityConstrVals] = smd2mp(xu, xl)
    prob         = smd2mp1(1, 2, 1);
    [fu, fuc]    = prob.evaluate_u(xu, xl);
    functionValue = fu;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function [functionValue equalityConstrVals inequalityConstrVals] = smd3mp(xu, xl)
    prob         = smd3mp1(1, 2, 1);
    [fu, fuc]    = prob.evaluate_u(xu, xl);
    functionValue = fu;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [functionValue equalityConstrVals inequalityConstrVals] = smd4mp(xu, xl)
    prob         = smd4mp1(1, 2, 1);
    [fu, fuc]    = prob.evaluate_u(xu, xl);
    functionValue = fu;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [functionValue equalityConstrVals inequalityConstrVals] = smd5mp(xu, xl)
    prob         = smd5mp1(1, 2, 1);
    [fu, fuc]    = prob.evaluate_u(xu, xl);
    functionValue = fu;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [functionValue equalityConstrVals inequalityConstrVals] = smd6mp(xu, xl)
    prob         = smd6mp1(1, 2, 1);
    [fu, fuc]    = prob.evaluate_u(xu, xl);
    functionValue = fu;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [functionValue equalityConstrVals inequalityConstrVals] = smd7mp(xu, xl)
    prob         = smd7mp1(1, 2, 1);
    [fu, fuc]    = prob.evaluate_u(xu, xl);
    functionValue = fu;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [functionValue equalityConstrVals inequalityConstrVals] = smd8mp(xu, xl)
    prob         = smd8mp1(1, 2, 1);
    [fu, fuc]    = prob.evaluate_u(xu, xl);
    functionValue = fu;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


