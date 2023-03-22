function [bestLX,bestLF,bestLC] = lowerLevelSearch(xu,level_Archive,BI)
global   llFE
%% Parameter setting
Archive = [];
record = [];
minN     = 4;
aRate    = 2.6;
MCR = zeros(20*BI.l_dim,1) + 0.2;
MF  = zeros(20*BI.l_dim,1) + 0.2;
k   = 1;
MOP = ones(1,5)/5;
maxIter = ceil(BI.LmaxFEs / BI.l_N);
ImprIter = ceil(BI.LmaxImprFEs / BI.l_N);
llFE     = 0;
best.stop = 0;
sqp.stop = 0;
tt       = 0;
%% Generate random population
[B_Init,Population] = Level_Initialization(xu,level_Archive,BI);
for iter = 1 : maxIter
    if(best.stop == 0 )
    % Reduce the population size
    N          = ceil((minN-BI.l_N)*llFE/BI.LmaxFEs) + BI.l_N; 
    [~,rank]   = sort(FitnessSingle(Population));
    Population = Population(rank(1:N));
    Archive    = Archive(randperm(end,min(end,ceil(aRate*N))));
    % Eigencoordinate System
    [~,I]  = sort(Population.objs,'ascend');
    num    = min(N,BI.l_dim);
    TopDec = Population(I(1:num)).decs;
    B      = eig(cov(TopDec));
      % Eigencoordinate System
    if isempty(B_Init)
    [~,I]  = sort(Population.objs,'ascend');
    num    = min(N,BI.l_dim);
    TopDec = Population(I(1:num)).decs;
    B1      = eig(cov(TopDec));
    else
      B1    =  B_Init;
    end
      R1= deal(zeros(size(B,2)));
    R1(logical(eye(size(B,2)))) = rand(1,size(B,2));  
        R2= deal(zeros(size(B1,2)));
    R2(logical(eye(size(B1,2)))) = rand(1,size(B1,2));
    
    % Generate parents, CR, F, and operator for each offspring
    Xp1 = Population(ceil(rand(1,N).*max(1,0.25*N))).decs;
    Xp2 = Population(ceil(rand(1,N).*max(2,0.5*N))).decs;
    Xr1 = Population(randi(end,1,N)).decs;
    Xr3 = Population(randi(end,1,N)).decs;
    P   = [Population,Archive];
    Xr2 = P(randi(end,1,N)).decs;
    CR  = randn(N,1).*sqrt(0.1) + MCR(randi(end,N,1));
    CR  = sort(CR);
    CR  = repmat(max(0,min(1,CR)),1,BI.l_dim);
    F   = min(1,trnd(1,N,1).*sqrt(0.1) + MF(randi(end,N,1)));
    while any(F<=0)
        F(F<=0) = min(1,trnd(1,sum(F<=0),1).*sqrt(0.1) + MF(randi(end,sum(F<=0),1)));
    end
    F  = repmat(F,1,BI.l_dim);
    OP = arrayfun(@(S)find(rand<=cumsum(MOP),1),1:N);
    OP = arrayfun(@(S)find(OP==S),1:length(MOP),'UniformOutput',false);
    % Generate offspring
    PopDec = Population.decs;
    OffDec = PopDec;
    OffDec(OP{1},:) = PopDec(OP{1},:) + F(OP{1},:).*(Xp1(OP{1},:)-PopDec(OP{1},:)+Xr1(OP{1},:)-Xr2(OP{1},:));
    OffDec(OP{2},:) = PopDec(OP{2},:) + F(OP{2},:)*B*R1*B'.*(Xp1(OP{2},:)-PopDec(OP{2},:)+Xr1(OP{2},:)-Xr3(OP{2},:));
    OffDec(OP{3},:) = F(OP{3},:)*B*R1*B'.*(Xr1(OP{3},:)+Xp2(OP{3},:)-Xr3(OP{3},:));
    OffDec(OP{4},:) = PopDec(OP{4},:) + F(OP{4},:)*B1*R2*B1'.*(Xp1(OP{4},:)-PopDec(OP{4},:)+Xr1(OP{4},:)-Xr3(OP{4},:));
    OffDec(OP{5},:) = F(OP{5},:)*B1*R2*B1'.*(Xr1(OP{5},:)+Xp2(OP{5},:)-Xr3(OP{5},:));
    if rand < 0.4
        Site = rand(size(CR)) > CR;
        OffDec(Site) = PopDec(Site);
    else
        p1 = randi(BI.l_dim,N,1);
        p2 = arrayfun(@(S)find([rand(1,BI.l_dim),2]>CR(S,1),1),1:N);
        for i = 1 : N
            Site = [1:p1(i)-1,p1(i)+p2(i):BI.l_dim];
            OffDec(i,Site) = PopDec(i,Site);
        end
    end
     %  offspring evaluation
     for i = 1 : N
        [F(i,:),C(i,:),PopDec(i,:)]= LL_evaluate(xu,OffDec(i,:),BI);  
     end
    Offspring = SOLUTION(BI,PopDec,F,C);
    % Update the population and archive
    delta   = FitnessSingle(Population) - FitnessSingle(Offspring);
    replace = delta > 0;
    Archive = [Archive,Population(replace)];
    Archive = Archive(randperm(end,min(end,ceil(aRate*N))));
    Population(replace) = Offspring(replace);
    % Update CR, F, and probabilities of operators
    if any(replace)
        w      = delta(replace)./sum(delta(replace));
        MCR(k) = (w'*CR(replace,1).^2)./(w'*CR(replace,1));
        MF(k)  = (w'*F(replace,1).^2)./(w'*F(replace,1));
        k      = mod(k,length(MCR)) + 1;
    else
        MCR(k) = 0.5;
        MF(k)  = 0.5;
    end
    delta = max(0,delta./abs(FitnessSingle(Population)));
    if any(cellfun(@isempty,OP))
        MOP = ones(1,5)/5;
    else
        MOP = cellfun(@(S)mean(delta(S)),OP);
        MOP = max(0.1,min(0.9,MOP./sum(MOP)));
    end
    if isempty(Population.best)
    else
        tt=tt+1;
        record(tt,1) =  Population.best.obj;
    if (tt>ImprIter && abs(record(tt,1)-record(tt-ImprIter+1,1)) <  1e-1)&&(sqp.stop == 0)
        x0   = Population.best.dec;
        [level_obj,lever_y, ~] = SQP(xu,x0,BI);
        if(level_obj<Population.best.obj)
            Population.best.obj =level_obj;
            Population.best.dec = lever_y;
            Population.best.con = 0;
            best.stop = 1;
        end
          sqp.stop = 1;
    end
     if (tt>ImprIter && abs(record(tt,1)-record(tt-ImprIter+1,1)) <  1e-5)
           best.stop = 1;
     end
    end
    else
    end
end
if isempty(Population.best)
  [~,Best] = min(Population.objs);
bestLX = Population(Best).dec;
bestLF = Population(Best).obj;
bestLC = Population(Best).con;
else
bestLX = Population.best.dec;
bestLF = Population.best.obj;
bestLC = Population.best.con;
end
end
function [F,C,xl] = LL_evaluate(xu,xl,BI)
for j = 1 : size(xl,2)
    xl(:,j) = max(min(xl(:,j),repmat(BI.l_ub(1,j),size(xl,1),1)),repmat(BI.l_lb(1,j),size(xl,1),1));
end
[F,~,C] = llTestProblem(xl,BI.fn,xu);
 C = sum(max(0,C));
end
function [B,Population] = Level_Initialization(xu,level_Archive,BI)
 PopDec1 = unifrnd(repmat(BI.l_lb,BI.l_N,1),repmat(BI.l_ub,BI.l_N,1));
 if(size(level_Archive,1)>=5)
     matrix1=xu;
     matrix2 =level_Archive(:,1:BI.u_dim);
     d = computeDistance(matrix1, matrix2);
     [~, closestParent] = sort( d,'ascend');
%       PopDec2 = unifrnd(repmat(BI.l_lb, 0.6*BI.l_N,1),repmat(BI.l_ub,0.6*BI.l_N,1));
     PopDec2 =level_Archive(closestParent(1:0.6*BI.l_N),BI.u_dim+1:end);
     PopDec_B =level_Archive(1:5,BI.u_dim+1:end);
%       if(size(level_Archive,1)>=10)
%            PopDec_B =level_Archive(1:10,BI.u_dim+1:end);
%       else
%      PopDec_B =level_Archive(1:end,BI.u_dim+1:end);
%        end
      % Eigencoordinate System
      B      = eig(cov(PopDec_B));
 else
      B      = [];
%   PopDec2 = level_Archive(randi(end,1,0.6*BI.l_N),1:BI.l_dim);
   PopDec2 = unifrnd(repmat(BI.l_lb, 0.6*BI.l_N,1),repmat(BI.l_ub,0.6*BI.l_N,1));
 end
  PopDect  =[PopDec2;PopDec1];
  PopDec  =PopDect(1:BI.l_N,:);
for i = 1 : BI.l_N
  [F(i,:),C(i,:),PopDec(i,:)]= LL_evaluate(xu,PopDec(i,:),BI);  
end
Population = SOLUTION(BI,PopDec,F,C);
end
function F = SQP_LLevaluate_F(xu,xl,BI)
global LC
for j = 1 : size(xl,2)
    xl(:,j) = max(min(xl(:,j),repmat(BI.l_ub(1,j),size(xl,1),1)),repmat(BI.l_lb(1,j),size(xl,1),1));
end
[F,~,C] = llTestProblem(xl,BI.fn,xu);
C = sum(max(0,C));
LC = C;
end

%% 下层函数得到最优化时对应的下层变量
function [level_obj,lever_y, level_funcCount] = SQP(xu,x0,BI)
global SQP_FEs
 options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
 options = optimoptions(options,'Display', 'off');
% options = optimoptions('fmincon','iter','Algorithm','sqp');
%  options = optimset('Display','off');
options.MaxFunctionEvaluations = SQP_FEs;
options.OptimalityTolerance =1e-10;
problem.options = options;
problem.solver = 'fmincon';
problem.objective = @(x)SQP_LLevaluate_F(xu,x,BI);
problem.nonlcon = @unitdisk2;
problem.x0 = x0;
problem.lb = BI.l_lb;
problem.ub = BI.l_ub;
[x,fval,exitflag,output] = fmincon(problem);
level_funcCount                = output.funcCount;
lever_y                        = x;
level_obj                      = fval;
end
function [c,ceq] = unitdisk2(xl)
 global LC
% global SQP_xu 
% global SQP_BI 
% for j = 1 : size(xl,2)
%     xl(:,j) = max(min(xl(:,j),repmat(SQP_BI.l_ub(1,j),size(xl,1),1)),repmat(SQP_BI.l_lb(1,j),size(xl,1),1));
% end
% [F,~,C] = llTestProblem(xl,SQP_BI.fn,SQP_xu);
c = LC';
ceq = [];
end
%% 计算两个上层解决方案之间的距离
function d = computeDistance(matrix1, matrix2)
 
    %Computes pairwise distance between rows of matrix1 and matrix2
    sz1 = size(matrix1, 1);
    sz2 = size(matrix2, 1);
    
    for i = 1:sz1
        for j = 1:sz2
            d(i,j) = sqrt(sum((matrix1(i,:)-matrix2(j,:)).^2));
        end
    end
end