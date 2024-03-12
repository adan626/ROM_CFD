%optimization file -2D opt
clear all; clc; close all;

vx = readmatrix('Data_Input.xlsx','Sheet','v_axial');
vx = vx(:,2:end);
vz = readmatrix('Data_Input.xlsx','Sheet','v_theta');
vz = vz(:,2:end);
exp_cond = readmatrix('Data_Input.xlsx','Sheet','Experiments');
Q = exp_cond(:,2);
rpm = exp_cond(:,3);
RTD = readmatrix('Data_Input.xlsx','Sheet','RTD');


%Training, validation sets
all_n=17;
runs = 1:1:17;


%%
% options=optimset('PlotFcns',@optimplotfval,'Display','Iter','MaxIter',50);
% y0=zeros(1,6);
% y0=[6.52E-06 7.34E-06 1.07E-05 6.60E-06 6.52E-06 6.68E-06 6.12E-06 5.34E-06 7.20E-06 5.5982e-06...
% 8.02E-06 7.86E-06 9.86E-06 8.16E-06 8.63E-06 5.598e-06 9.65E-06 1.4492];
% %y0=[0.000365426049552414	0.000463567519779019	0.000344137482312466	0.000407804592492688	0.000436735488003996	3.68896365859823];
% tic;
% func = @(y)obj(y,all_n,runs,vx,vz,Q,RTD);
% [yopt,fval,flag] = fminsearch(func,y0,options)

% yopt = fmincon(func,y0,[],[],[],[],[1e-6 0.1],[1e-4 10],[],options)
% opt_toc=toc;
% % 
% for i=1:all_n
%     func = @(y)obj(y,1,i,vx,vz,Q,RTD);
%     [yopt,fval,flag] = fminsearch(func,y0(i),options);
%     D_col(i)=yopt;
% end



%% predictions
all_n=3;
runs=[6 9 11];
y0=[6.52E-06 7.34E-06 1.07E-05 6.60E-06 6.52E-06 6.68E-06 6.12E-06 5.34E-06 7.20E-06 5.5982e-06...
8.02E-06 7.86E-06 9.86E-06 8.16E-06 8.63E-06 5.598e-06 9.65E-06 1.4492];
MRT_model=zeros(all_n,1);
MRT_exp=zeros(all_n,1);
vari_model=zeros(all_n,1);
vari_exp=zeros(all_n,1);
skew_model=zeros(all_n,1);
skew_exp=zeros(all_n,1);
E_model(1:10000,1:all_n)=0;
time(1:10000,1:all_n)=0;
for i=1:all_n
    j=runs(i);
    time_exp=RTD(:,2*j-1);
    time_exp=time_exp(~isnan(time_exp));
    Et_exp=RTD(:,2*j);
    Et_exp=Et_exp(~isnan(Et_exp));
    yopt_=[y0(j) y0(end)];
    [pred] = CD_Model_2D(yopt_,vx(:,j),vz(:,j),Q(j),time_exp,Et_exp,'run');
    RSD_model(1:length(pred),i)=pred(:,1);
    time(1:length(pred),i)=pred(:,2);
    
    % MRT_model(i)=pred(1);
    % MRT_exp(i)=pred(2);
    % vari_model(i)=pred(3);
    % vari_exp(i)=pred(4);
    % skew_model(i)=pred(5);
    % skew_exp(i)=pred(6);
end
%%

figure
plot(time(1:3000,1),RSD_model(1:3000,1),'-b','Linewidth',2); hold on
plot(time(1:2852,2),RSD_model(1:2852,2),'-r','Linewidth',2);hold on
plot(time(1:2000,3),RSD_model(1:2000,3),'-g','Linewidth',2);hold on
legend({'Run 6','Run 9','Run 11'},'FontSize',12)

%legend({'Run 6-CFD','Run 6-ROM','Run9-CFD','Run9-ROM','Run11-CFD','Run11-ROM'},'FontSize',12)
xlim([0 250])
xlabel('Time, s','FontSize',20)
ylabel('RSD','FontSize',20)
%title('RTD Comparison for Validation Runs','FontSize',20)
set(gcf,'color','white')

z=1:1:10;
figure
plot(z,vz(:,6),'-b','Linewidth',2); hold on
plot(z,vz(:,9),'-r','Linewidth',2);hold on
plot(z,vz(:,11),'-g','Linewidth',2);hold on
legend({'Run 6','Run 9','Run 11'},'FontSize',12)
xlabel('Discretized bin (j)','FontSize',20)
ylabel('Velocity','FontSize',20)
ylim([0.01 0.04])
set(gcf,'color','white')


% 
% ParityPlot(MRT_model,MRT_exp)
% title('MRT','fontsize',15) 
% 
% ParityPlot(vari_model,vari_exp)
% title('Variance','fontsize',15) 
% 
% 
% ParityPlot(skew_model,skew_exp)
% title('Skewness','fontsize',15)
%}

%%
function [err_tot] = obj(y0,n,runs,vx,vz,Q,RTD)
err_tot=0;
for i=1:n
    j=runs(i);
    time_exp=RTD(:,2*j-1);
    time_exp=time_exp(~isnan(time_exp));
    Et_exp=RTD(:,2*j);
    Et_exp=Et_exp(~isnan(Et_exp));
    yopt=[y0(i) y0(end)];
    [err] = CD_Model_2D(yopt,vx(:,j),vz(:,j),Q(j),time_exp,Et_exp,'fit');
    err_tot=err_tot+err;

end

end




