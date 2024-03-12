clear all; clc;
exp_cond = readmatrix('Data_Input.xlsx','Sheet','Experiments');
Q = exp_cond(:,2);
rpm = exp_cond(:,3);
RTD = readmatrix('Data_Input.xlsx','Sheet','RTD');
dt=0.05; 
for i=1:17
    time_exp=RTD(:,2*i-1);
    time_exp=time_exp(~isnan(time_exp));
    Et_exp=RTD(:,2*i);
    Et_exp=Et_exp(~isnan(Et_exp));
    time_mod=0:dt:time_exp(end);
    Et_mod=spline(time_exp,Et_exp,time_mod);
    MRT_exp(i)=trapz(time_exp,time_exp.*Et_exp);
    %var_exp(i)=var(Et_exp);
    %var_exp(i)=trapz(time_exp,(time_exp-MRT_exp(i)).^2.*Et_exp);
    var_exp(i)=trapz(time_mod,(time_mod-MRT_exp(i)).^2.*Et_mod);
    denom(i)=trapz(time_mod,Et_mod);
    norm_var(i)=var_exp(i)/(MRT_exp(i)^2);
    %skew_exp(i)=skewness(Et_exp);
    skew_exp(i)=1/(sqrt(var_exp(i))^(3/2))*trapz(time_exp,(time_exp-MRT_exp(i)).^3.*Et_exp);
    norm_skew(i)=skew_exp(i)/(MRT_exp(i)^3);
    % if Et_exp(end)>1e-6
    %     i
    %     plot(time_exp,Et_exp);hold on;
    % end
end

%%
% mcv=(var_exp./denom)';
% norm_var=rescale(norm_var)';
var_exp_=round(rescale(var_exp)',2);
tbl = table(Q,rpm,var_exp_);
figure
heatmap(tbl,'Q','rpm','ColorVariable','var_exp_','FontSize',14,Colormap=sky(4))
ylabel('Screw Speed, rpm')
xlabel('Throughput, Q')
title('Variance')

MRT_exp_=round(rescale(MRT_exp)',2);
tbl = table(Q,rpm,MRT_exp_);
figure
heatmap(tbl,'Q','rpm','ColorVariable','MRT_exp_','FontSize',14,Colormap=sky(4))
ylabel('Screw Speed, rpm')
xlabel('Throughput, Q')
title('MRT')

%%

T=array2table(exp_cond,'VariableNames',{'Run','Q','rpm'});
T.MRT=MRT_exp';
T.var=var_exp';
pivotT=pivot(T,Columns='Q',Rows='rpm',DataVariable='MRT',Method='mean');
pivotT2=pivot(T,Columns='Q',Rows='rpm',DataVariable='var',Method='mean');
MRT_arr=table2array(pivotT);
var_arr=table2array(pivotT2);
%%
figure
plot(MRT_arr(:,2),var_arr(:,2),'-o',MRT_arr(:,3),var_arr(:,3),'-o',...
    MRT_arr(:,4),var_arr(:,4),'-o',MRT_arr(:,5),var_arr(:,5),'-o','Linewidth',2)
leg = legend({'3e-7','4e-7','5e-7','6e-7',''},'FontSize',12,'Location','northwest');
htitle=get(leg,'Title');
set(htitle,'String','Throughput');
xlabel('Mean Residence Time (MRT), s','FontSize',17)
ylabel('Variance','FontSize',17)

hold on;
xmin=70;
xmax=110;
ymin=150;
ymax=280;
r = rectangle('Position',[xmin ymin xmax-xmin ymax-ymin],'FaceColor', [1, 0, 0, 0.3], ...
                'EdgeColor', [1, 0, 0, 0.3]);

annotation('textbox',[.25 0.2 0 0],'string','Steady Variance Zone','FitBoxToText','on','EdgeColor','none','FontSize',12) 

annotation('textbox',[.67 0.55 0 0],'string','From left to right:','FitBoxToText','on','EdgeColor','none','FontSize',11) 
annotation('textbox',[.67 0.50 0 0],'string','High to low RPM','FitBoxToText','on','EdgeColor','none','FontSize',11) 

%%
a=table2array(pivotT);
figure
plot([3 4 5 6],a(:,2:end),'-o')

legend('45','55','65','75')

xx=[45 55 65 75];
yy=[3 4 5 6];


%%
%{
figure;
plot(RTD(:,1),RTD(:,2),RTD(:,3),RTD(:,4),RTD(:,5),RTD(:,6),RTD(:,7),RTD(:,8))
legend('45','55','65','75')

figure;
plot(RTD(:,9),RTD(:,10),RTD(:,11),RTD(:,12),RTD(:,13),RTD(:,14),RTD(:,15),RTD(:,16))
legend('3','4','5','6')

figure;
plot(RTD(:,9),RTD(:,10),RTD(:,17),RTD(:,18),RTD(:,19),RTD(:,20),RTD(:,21),RTD(:,22))
legend('45','55','65','75')

figure;
plot(RTD(:,21),RTD(:,22),RTD(:,33),RTD(:,34),RTD(:,7),RTD(:,8),RTD(:,27),RTD(:,28))
legend('3','4','5','6')

figure
plot(RTD(:,33),RTD(:,34))

figure
plot(RTD(:,9),RTD(:,10),RTD(:,31),RTD(:,32))
%}