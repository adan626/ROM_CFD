%% File for PP optimization
clear all;clc;
rng default
options = optimoptions('gamultiobj','PopulationSize',60,...
          'ParetoFraction',0.35,'PlotFcn',@gaplotpareto);
lb = [0.3 0.45];
ub = [0.6 0.75];
tic
[solution,ObjectiveValue] = gamultiobj(@run_prediction,2,...
                          [],[],[],[],lb,ub,options);

paretotoc=toc

%%
figure
scatter(ObjectiveValue(:,1),-ObjectiveValue(:,2),'^b','LineWidth',2);
hold on
xlabel('Mean Residence Time,s','FontSize',20)
ylabel('Variance','FontSize',20)
title('Pareto Optimization','FontSize',20)
set(gca, 'YDir', 'reverse');

vals_sorted=sortrows(ObjectiveValue,1);
yy=smooth(vals_sorted(:,1),-vals_sorted(:,2));
%yy=smooth(vals_sorted(:,1),-vals_sorted(:,2),0.8,'loess')

plot(xx,yy,'--b','Linewidth',1.5)
legend({'Optimized Value','Moving Average'},'FontSize',12,'Location','northeast')

%% PP optimization weighted function
clear all;clc;

options=optimset('PlotFcns',@optimplotfval,'Display','Iter','MaxIter',50);
y0=[0.5 0.6];
y0=[0.45 0.702];
lb = [0.3 0.45];
ub = [0.7 0.75];

% [yopt_,fval,flag] = fminsearch(@run_prediction,y0,options);
% options=optimset('PlotFcns',@optimplotfval,'Display','Iter','Algorithm','sqp',...
%     'TolX',1e-10,'TolFun',1e-10);
% options=optimoptions('fmincon','StepTolerance',1e-10,'Display','Iter');
% [yopt,fval,flag,output]=fmincon(@run_prediction,y0,[],[],[],[],lb,ub,[],options);


tic
rng default % For reproducibility
opts = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','PlotFcns',@optimplotfval,'Maxiter',20);
problem = createOptimProblem('fmincon','objective',...
    @run_prediction,'x0',y0,'lb',lb,'ub',ub,'options',opts);

ms = MultiStart('FunctionTolerance',2e-4,'UseParallel',true)
gs = GlobalSearch(ms)
gs.MaxTime = 1800;
[x,f] = run(gs,problem)
gs=toc


%%
function out = run_prediction(y)

    Q=y(1)*1e-7*10;
    rpm=y(2)*100;
    z = readmatrix('Data_Input.xlsx','Sheet','v_axial','Range','A:A');
    z = z(2:end);
    vz = [0.0029, 0.0030];
    %vx pred
    X_valid = [z repmat(Q,[181 1]).*1e7 repmat(rpm,[181 1])];

    load vxmdl2.mat
    predvx = predict(mdl, X_valid);
    
    % D
    v=mean(predvx);
    predD=(-6.7387e-06) - (33.7226.*Q) - (4.0753e-07.*rpm)+(0.025795.*v) + ...
        (Q-4.786e-7).*((v-0.0021)*(-1.7119e+04)) + (Q - 4.786e-07).*((Q - 4.786e-07).*1.2983e+08);
    vx=predvx;
    D=predD;
    
    yopt=[D; 1.4492];
    [pred] = CD_Model(yopt,vx,vz,Q,0,0,'opt');
    MRT=pred(1);
    var=pred(2);
    skew=pred(3);
    Pe=pred(4);
    out =[MRT -var];
    % 
    % out = abs(MRT-90)/90 +abs(var-500)/500;
end
