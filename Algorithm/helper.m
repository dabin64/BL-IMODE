function [helper_UX] = helper(BI,helper_iter)
%% Parameter setting
Archive = [];
minN     = 4;
aRate    = 2.6;
MCR = zeros(20*(BI.l_dim+BI.u_dim),1) + 0.2;
MF  = zeros(20*(BI.l_dim+BI.u_dim),1) + 0.2;
k   = 1;
MOP = ones(1,3)/3;
%% Generate random population
[Population] = helper_Initialization(BI);
for iter = 1 : helper_iter
    % Reduce the population size
    N          =  BI.u_N; 
    [~,rank]   = sort(FitnessSingle(Population));
    Population = Population(rank(1:N));
    Archive    = Archive(randperm(end,min(end,ceil(aRate*N))));
    % Generate parents, CR, F, and operator for each offspring
    Xp1 = Population(ceil(rand(1,N).*max(1,0.25*N))).decs;
    Xp2 = Population(ceil(rand(1,N).*max(2,0.5*N))).decs;
    Xr1 = Population(randi(end,1,N)).decs;
    Xr3 = Population(randi(end,1,N)).decs;
    P   = [Population,Archive];
    Xr2 = P(randi(end,1,N)).decs;
    CR  = randn(N,1).*sqrt(0.1) + MCR(randi(end,N,1));
    CR  = sort(CR);
    CR  = repmat(max(0,min(1,CR)),1,(BI.l_dim+BI.u_dim));
    F   = min(1,trnd(1,N,1).*sqrt(0.1) + MF(randi(end,N,1)));
    while any(F<=0)
        F(F<=0) = min(1,trnd(1,sum(F<=0),1).*sqrt(0.1) + MF(randi(end,sum(F<=0),1)));
    end
    F  = repmat(F,1,(BI.l_dim+BI.u_dim));
    OP = arrayfun(@(S)find(rand<=cumsum(MOP),1),1:N);
    OP = arrayfun(@(S)find(OP==S),1:length(MOP),'UniformOutput',false);
    % Generate offspring
    PopDec = Population.decs;
    OffDec = PopDec;
    OffDec(OP{1},:) = PopDec(OP{1},:) + F(OP{1},:).*(Xp1(OP{1},:)-PopDec(OP{1},:)+Xr1(OP{1},:)-Xr2(OP{1},:));
    OffDec(OP{2},:) = PopDec(OP{2},:) + F(OP{2},:).*(Xp1(OP{2},:)-PopDec(OP{2},:)+Xr1(OP{2},:)-Xr3(OP{2},:));
    OffDec(OP{3},:) = F(OP{3},:).*(Xr1(OP{3},:)+Xp2(OP{3},:)-Xr3(OP{3},:));
    if rand < 0.4
        Site = rand(size(CR)) > CR;
        OffDec(Site) = PopDec(Site);
    else
        p1 = randi((BI.l_dim+BI.u_dim),N,1);
        p2 = arrayfun(@(S)find([rand(1,(BI.l_dim+BI.u_dim)),2]>CR(S,1),1),1:N);
        for i = 1 : N
            Site = [1:p1(i)-1,p1(i)+p2(i):(BI.l_dim+BI.u_dim)];
            OffDec(i,Site) = PopDec(i,Site);
        end
    end
     %  offspring evaluation
     for i = 1 : N
        [F(i,:),C(i,:),PopDec(i,:)]= UL_evaluate(OffDec(i,1:BI.u_dim),OffDec(i,BI.u_dim+1:end),BI);  
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
        MOP = ones(1,3)/3;
    else
        MOP = cellfun(@(S)mean(delta(S)),OP);
        MOP = max(0.1,min(0.9,MOP./sum(MOP)));
    end
end 
for i =1:N
helper_LX(i,:) = Population(i).decs;
end
helper_UX  = helper_LX(:,1:BI.u_dim);
end
function [Population] = helper_Initialization(BI)
 PopDec = [unifrnd(repmat(BI.u_lb,BI.u_N,1),repmat(BI.u_ub,BI.u_N,1)),unifrnd(repmat(BI.l_lb,BI.u_N,1),repmat(BI.l_ub,BI.u_N,1))];
for i = 1 : BI.u_N
  [F(i,:),C(i,:),PopDec(i,:)]= UL_evaluate(PopDec(i,1:BI.u_dim),PopDec(i,BI.u_dim+1:end),BI);  
end
Population = SOLUTION(BI,PopDec,F,C);
end
function [F,C,POPdec] = UL_evaluate(UPop,LPOP,BI)
for j = 1 : size(UPop,2)
    UPop(:,j) = max(min(UPop(:,j),repmat(BI.u_ub(1,j),size(UPop,1),1)),repmat(BI.u_lb(1,j),size(UPop,1),1));
end
for j = 1 : size(LPOP,2)
    LPOP(:,j) = max(min(LPOP(:,j),repmat(BI.l_ub(1,j),size(LPOP,1),1)),repmat(BI.l_lb(1,j),size(LPOP,1),1));
end
POPdec=[UPop,LPOP];
[F,~,C] = ulTestProblem(UPop, LPOP, BI.fn);
C = sum(max(0,C));
end