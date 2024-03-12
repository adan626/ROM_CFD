%optimization file
clear all; clc; close all;

vx = readmatrix('Data_Input.xlsx','Sheet','v_axial');
vx = vx(:,2:end);
% vz = readmatrix('Data_Input.xlsx','Sheet','v_theta');
% vz = vz(:,2:end);
vz = [0.0029, 0.0030];
exp_cond = readmatrix('Data_Input.xlsx','Sheet','Experiments');
Q = exp_cond(:,2);
rpm = exp_cond(:,3);
RTD = readmatrix('Data_Input.xlsx','Sheet','RTD');


%Training, validation sets
train_n=14;
test_n=3;
train_runs=[1 2 3 4 5 7 8 10 12 13 14 15 16 17];
test_runs=[6 9 11];
all_n=17;
runs = 1:1:17;


%%
options=optimset('PlotFcns',@optimplotfval,'Display','Iter','MaxIter',50 ...
    ,'MaxFunEvals',2d3,'Algorithm','SQP');
%y0=1.4492;
y0=[6.52E-06 7.34E-06 1.07E-05 6.60E-06 6.52E-06 6.68E-06 6.12E-06 5.34E-06 7.20E-06 5.5982e-06...
8.02E-06 7.86E-06 9.86E-06 8.16E-06 8.63E-06 5.598e-06 9.65E-06 1.4492];
y0(1:17)=5e-6;
y0(18)=1.5;

tic;
func = @(y)obj(y,1,1,vx,vz,Q,RTD);
[yopt_,fval,flag] = fminsearch(func,y0,options);
opt_toc=toc;
% 
% for i=1:all_n
%     func = @(y)obj(y,1,i,vx,vz,Q,RTD);
%     [yopt,fval,flag] = fminsearch(func,y0(i),options);
%     D_col(i)=yopt;
% end



%% predictions

y0=[6.52E-06 7.34E-06 1.07E-05 6.60E-06 6.52E-06 6.68E-06 6.12E-06 5.34E-06 7.20E-06 5.5982e-06...
8.02E-06 7.86E-06 9.86E-06 8.16E-06 8.63E-06 5.598e-06 9.65E-06 1.4492];
MRT_model=zeros(train_n,1);
MRT_exp=zeros(train_n,1);
vari_model=zeros(train_n,1);
vari_exp=zeros(train_n,1);
skew_model=zeros(train_n,1);
skew_exp=zeros(train_n,1);
for i=1:train_n
    j=train_runs(i);
    time_exp=RTD(:,2*j-1);
    time_exp=time_exp(~isnan(time_exp));
    Et_exp=RTD(:,2*j);
    Et_exp=Et_exp(~isnan(Et_exp));
    yopt_=[y0(j) y0(end)];
    [pred] = CD_Model(yopt_,vx(:,j),vz,Q(j),time_exp,Et_exp,'run');
    MRT_model(i)=pred(1);
    MRT_exp(i)=pred(2);
    vari_model(i)=pred(3);
    vari_exp(i)=pred(4);
    skew_model(i)=pred(5);
    skew_exp(i)=pred(6);
end


ParityPlot(MRT_model,MRT_exp)
title('MRT','fontsize',15) 

ParityPlot(vari_model,vari_exp)
title('Variance','fontsize',15) 


ParityPlot(skew_model,skew_exp)
title('Skewness','fontsize',15)


%%
function [err_tot] = obj(y0,n,runs,vx,vz,Q,RTD)
err_tot=0;
for i=1:n
    j=runs(i);
    time_exp=RTD(:,2*j-1);
    time_exp=time_exp(~isnan(time_exp));
    Et_exp=RTD(:,2*j);
    Et_exp=Et_exp(~isnan(Et_exp));
    yopt=[y0(j) y0(end)];
    [err] = CD_Model(yopt,vx(:,j),vz,Q(j),time_exp,Et_exp,'fit');
    err_tot=err_tot+err;

end

end






