

pwd
% parpool('local',4);
%  parpool('threads');
maxRuns = 30;
BI_list = [];
% BI_list = [BI_list;getBLOPinfo('SMDmp',1:8,5)];
% BI_list = [BI_list;getBLOPinfo('TP',9:10)];
BI_list = [BI_list;getBLOPinfo('SMD',1:12,5)];
% BI_list = [BI_list;getBLOPinfo('SMD',1:12,20)];

% BI_list = [BI_list;getBLOPinfo('GoldMining')];
% BI_list = [BI_list;getBLOPinfo('DecisionMaking')];


algName = 'BLIMODE';
all_output_data =[];
all_output_data1 =[];
t=0;
j=0;
for BI = BI_list'
    if BI.dim == 3
        % for gold mining problem
        BI.UmaxFEs = 350;
        BI.UmaxImprFEs = 150;
        BI.LmaxFEs = 300;
        BI.LmaxImprFEs = 10;
    elseif  BI.dim == 6
        BI.UmaxFEs = 650;
        BI.UmaxImprFEs = 250;
        BI.LmaxFEs = 350;
        BI.LmaxImprFEs = 15;
    elseif BI.dim == 5 || strcmp(BI.fn(1:2),'tp')
        % for generic benchmark test problem or decision making problem 
        BI.UmaxFEs = 2500;
        BI.UmaxImprFEs = 250;
        BI.LmaxFEs = 350;
        BI.LmaxImprFEs = 25;
    elseif BI.dim == 10
        BI.UmaxFEs = 3500;
        BI.UmaxImprFEs = 500;
        BI.LmaxFEs = 350;
        BI.LmaxImprFEs = 30;
    elseif BI.dim == 20
        BI.UmaxFEs = 5000;
        BI.UmaxImprFEs = 750;
        BI.LmaxFEs = 500;
        BI.LmaxImprFEs = 50;
    else
        error('unknown dimensionality');
    end
   
    for runNo = 1:maxRuns
        t=t+1;
		tic;
        BI.u_N =5;
        BI.l_N =5;
        BI.levelArchive_N = 10;
        ins = BLIMODE(BI);
		ins.runPath = pwd;
		ins.runNo = runNo;
		ins.BI = BI;
		ins.alg = algName;
%         fprintf('%s %s #%d [%g,%g] %d %d %d %d  \n', ins.alg, ins.BI.fn, ins.runNo, ins.UF, ins.LF,ins.UFEs,ins.LFEs,ins.UX,ins.LX);
        fprintf('%s %s #%d [%g,%g] %d %d\n', ins.alg, ins.BI.fn, ins.runNo, ins.UF, ins.LF,ins.UFEs,ins.LFEs);
        all_output_data(t,1)=abs(ins.UF-BI.u_fopt);
        all_output_data(t,2)=abs(ins.LF-BI.l_fopt);
        all_output_data(t,3)=ins.UFEs;
        all_output_data(t,4)=ins.LFEs;
    end	
    j=j+1;
       for y=1:4
       all_output_data1(j,2*y-1)=prctile(all_output_data((j-1)*30+1:30*j,y),50);
       all_output_data1(j,2*y)=prctile(all_output_data((j-1)*30+1:30*j,y),75)-prctile(all_output_data((j-1)*30+1:30*j,y),25);
       end
end
xlswrite('BL-IMODE-5-IQR.xls',all_output_data1);
save('BL-IMODE-5.mat','all_output_data');



