function ins = BLIMODE(BI)
global ulFunctionEvaluations;
global llFunctionEvaluations;
global UUFE;
global ACC_ER;
ulFunctionEvaluations = 0;
llFunctionEvaluations = 0;
%% Parameter setting
Archive = [];
elite = [];
record = [];
Population_model =[];
minN     = 4;
aRate    = 2.6;
MCR = zeros(20*BI.u_dim,1) + 0.2;
MF  = zeros(20*BI.u_dim,1) + 0.2;
k   = 1;
MOP = ones(1,3)/3;
maxIter = ceil(BI.UmaxFEs / BI.u_N);
ImprIter = ceil(BI.UmaxImprFEs / BI.u_N);
best.stop = 0;
% [helper_UX] = helper(BI,1);
%% Generate random population
[Population,level_Archive] = Upper_Initialization(BI);
UUFE = 0;
tt=0;
for iter = 1 : maxIter
    if(best.stop == 0 )
    % Reduce the population size
    N          = ceil((minN-BI.u_N)*UUFE/BI.UmaxFEs) + BI.u_N; 
    [~,rank]   = sort(FitnessSingle(Population));
    Population = Population(rank(1:N));
    Archive    = Archive(randperm(end,min(end,ceil(aRate*N))));
    % Eigencoordinate System
    [~,I]  = sort(Population.objs,'ascend');
    num    = min(N,BI.u_dim);
    TopDec = Population(I(1:num)).decs;
    B      = eig(cov(TopDec));
    R1= deal(zeros(size(B,2)));
    R1(logical(eye(size(B,2)))) = rand(1,size(B,2));
    % Generate parents, CR, F, and operator for each offspring
    Xp1 = Population(ceil(rand(1,N).*max(1,0.25*N))).decs;
    Xp2 = Population(ceil(rand(1,N).*max(2,0.5*N))).decs;
    Xr1 = Population(randi(end,1,N)).decs;
    Xr3 = Population(randi(end,1,N)).decs;
    P   = [Population,Archive];
    Xr2 = P(randi(end,1,N)).decs;
    CR  = randn(N,1).*sqrt(0.1) + MCR(randi(end,N,1));
    CR  = sort(CR);
    CR  = repmat(max(0,min(1,CR)),1,BI.u_dim);
    F   = min(1,trnd(1,N,1).*sqrt(0.1) + MF(randi(end,N,1)));
    while any(F<=0)
        F(F<=0) = min(1,trnd(1,sum(F<=0),1).*sqrt(0.1) + MF(randi(end,sum(F<=0),1)));
    end
    F  = repmat(F,1,BI.u_dim);
    OP = arrayfun(@(S)find(rand<=cumsum(MOP),1),1:N);
    OP = arrayfun(@(S)find(OP==S),1:length(MOP),'UniformOutput',false);
    % Generate offspring
    PopDec = Population.decs;
    OffDec = PopDec;
    OffDec(OP{1},:) = PopDec(OP{1},:) + F(OP{1},:).*(Xp1(OP{1},:)-PopDec(OP{1},:)+Xr1(OP{1},:)-Xr2(OP{1},:));
    OffDec(OP{2},:) = PopDec(OP{2},:) + F(OP{2},:)*B*R1*B'.*(Xp1(OP{2},:)-PopDec(OP{2},:)+Xr1(OP{2},:)-Xr3(OP{2},:));
    OffDec(OP{3},:) = F(OP{3},:)*B*R1*B'.*(Xr1(OP{3},:)+Xp2(OP{3},:)-Xr3(OP{3},:));
    if rand < 0.4
        Site = rand(size(CR)) > CR;
        OffDec(Site) = PopDec(Site);
    else
        p1 = randi(BI.u_dim,N,1);
        p2 = arrayfun(@(S)find([rand(1,BI.u_dim),2]>CR(S,1),1),1:N);
        for i = 1 : N
            Site = [1:p1(i)-1,p1(i)+p2(i):BI.u_dim];
            OffDec(i,Site) = PopDec(i,Site);
        end
    end
%     OffDec
     %  offspring evaluation
    [Offspring] = Offspring_solution(OffDec,BI,level_Archive);
    [Population,Offspring,elite] = Cal_solution(Population,Offspring,BI,elite,level_Archive);
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
%     %%更新下层精英库
%     for i =1:N
%         level_Archive_LX(i,:) = [Population(i).dec,Population(i).lx];
%     end
%     level_Archive = [level_Archive;level_Archive_LX];
%     level_Archive = level_Archive(randperm(end,min(end,BI.levelArchive_N)),:);
    Population_model          = [Population,Population_model];
    Builtmodel_data           = [Population_model.decs,Population_model.lxs,Population_model.objs];
    Builtmodel_data_uni        = unique(Builtmodel_data,'rows');
     [~,rank]                  = sort(Builtmodel_data_uni(:,end));
    build_model_data           = Builtmodel_data_uni(rank(1:end),1:end-1);
    level_Archive              = build_model_data;
      %%停止条件
    if isempty(Population.best)
    else
        tt=tt+1;
        record(tt,1) =  Population.best.obj;
    if (tt>ImprIter && abs(record(tt,1)-record(tt-ImprIter+1,1)) <  ACC_ER)
         best.stop = 1;
    end
    end
    if isempty(Population.best)
    else
        reachTheTarget_a = abs(Population.best.obj-BI.u_fopt)< ACC_ER;
        reachTheTarget_b = abs(Population.best.lf-BI.l_fopt)< ACC_ER ;
        if(reachTheTarget_a)&&(reachTheTarget_b)
            best.stop = 1;
        end
    end
    else
    end
end
if isempty(Population.best)
    [~,Best] = min(Population.objs);
    ins.UF = Population(Best).obj;
    ins.LF = Population(Best).lf;
    ins.UX = Population(Best).dec;
    ins.LX = Population(Best).lx;
    ins.UFEs = ulFunctionEvaluations;
    ins.LFEs = llFunctionEvaluations;
else
    ins.UF = Population.best.obj;
    ins.LF = Population.best.lf;
    ins.UX = Population.best.dec;
    ins.LX = Population.best.lx;
    ins.UFEs = ulFunctionEvaluations;
    ins.LFEs = llFunctionEvaluations;
end
end
%%产生子代
function[Offspring] = Offspring_solution(PopDec,BI,level_Archive)
for i = 1 : size(PopDec,1)
    [bestLX(i,:),bestLF(i,:),bestLC(i,:)] = lowerLevelSearch(PopDec(i,:),level_Archive,BI);
     [F(i,:),C(i,:),UX(i,:),LX(i,:)]= UL_evaluate(PopDec(i,:),bestLX(i,:),BI);  
end
C=[C,bestLC];
Offspring = SOLUTION(BI,UX,F,C,LX,bestLF);
end
%%父代子代于最优解比较
function[Population,Offspring,elite] = Cal_solution(Population,Offspring,BI,elite,level_Archive)
if (isempty(Offspring.best))||(isempty(Population.best))
else
    if isempty(elite)
        if isempty(Offspring.best)&&(isempty(Population.best)==0)
            elite       = Population.best;
        end
        if (isempty(Offspring.best)==0)&&isempty(Population.best)
            elite       = Offspring.best;
        end
        if (isempty(Offspring.best)==0)&&(isempty(Population.best)==0)
            if(Offspring.best.obj<Population.best.obj)
                elite       = Offspring.best;
            else
                elite       = Population.best;
            end
        end
        elite     = refine(elite,BI,level_Archive) ;
    else
        if upperLevelComparator(Offspring,elite)
            RF_idx = find(Offspring.objs<elite.obj);
            for i=1:size(RF_idx,1)
                a = refine(Offspring(RF_idx(i)),BI,level_Archive);
                Offspring(RF_idx(i)).obj =a.obj;
            end
            if upperLevelComparator(Offspring,elite)
                elite = Offspring.best;
            end
        end
        if upperLevelComparator(Population,elite)
            RF_idx = find(Population.objs<elite.obj);
            for i=1:size(RF_idx,1)
                a = refine(Population(RF_idx(i)),BI,level_Archive);
                Population(RF_idx(i)).obj =a.obj;
            end
            if upperLevelComparator(Population,elite)
                elite = Population.best;
            end
        end
        if rand>0.5
            elite = refine(elite,BI,level_Archive);
        end
    end
end
end

%%初始化种群
function[Population,level_Archive] = Upper_Initialization(BI)
PopDec1 = unifrnd(repmat(BI.u_lb,BI.u_N,1),repmat(BI.u_ub,BI.u_N,1));
% PopDec2 = help_dec;
% level_Archive = unifrnd(repmat(BI.l_lb, BI.levelArchive_N,1),repmat(BI.l_ub, BI.levelArchive_N,1));
level_Archive = [];
for i = 1 : BI.u_N
    [bestLX(i,:),bestLF(i,:),bestLC(i,:)] = lowerLevelSearch(PopDec1(i,:),level_Archive,BI);
     [F(i,:),C(i,:),UX(i,:),LX(i,:)]= UL_evaluate(PopDec1(i,:),bestLX(i,:),BI);  
end
C=[C,bestLC];
Population = SOLUTION(BI,UX,F,C,LX,bestLF);
for i =1:BI.u_N
level_Archive_LX(i,:) = [Population(i).dec,Population(i).lx];
end
level_Archive = [level_Archive;level_Archive_LX];
level_Archive = level_Archive(randperm(end,min(end,BI.levelArchive_N)),:);
end

%%
function [F,C,UX,LX] = UL_evaluate(UPop,LPOP,BI)
for j = 1 : size(UPop,2)
    UPop(:,j) = max(min(UPop(:,j),repmat(BI.u_ub(1,j),size(UPop,1),1)),repmat(BI.u_lb(1,j),size(UPop,1),1));
end
for j = 1 : size(LPOP,2)
    LPOP(:,j) = max(min(LPOP(:,j),repmat(BI.l_ub(1,j),size(LPOP,1),1)),repmat(BI.l_lb(1,j),size(LPOP,1),1));
end
UX=UPop;
LX=LPOP;
[F,~,C] = ulTestProblem(UPop, LPOP, BI.fn);
C = sum(max(0,C));
end
function [isNoWorseThan] = upperLevelComparator(P,Q)
if isempty(P)
    isNoWorseThan = flase;
else
    isNoWorseThan = P.best.obj <= Q.best.obj;
end
end
function Q = refine(P,BI,level_Archive)
Q = P;
[R.LX,Q.lf,R.LC] = lowerLevelSearch(Q.dec,level_Archive,BI);
 if lowerLevelComparator(Q,P)&&(R.LC<=0)
    [Q.obj,~,Q.dec,Q.lx]= UL_evaluate(Q.dec,R.LX,BI);
else
    Q = P;
 end
end
function isNoWorseThan = lowerLevelComparator(Q,P)
    isNoWorseThan = Q.lf <= P.lf;
end