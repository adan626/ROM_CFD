clear all;%clc;
RTD = readmatrix('Data_Input.xlsx','Sheet','RTD');
exp_cond = readmatrix('Data_Input.xlsx','Sheet','Experiments');
Q_arr = exp_cond(:,2);
rpm_arr = exp_cond(:,3);
z = readmatrix('Data_Input.xlsx','Sheet','v_axial','Range','A:A');
z = z(2:end);%vx = vx(:,2:end);


vz = [0.0029, 0.0030];

runs=[6,9,11];
n=3;
% runs=1:1:17;
% n=17;
E_model(1:10000,1:n)=0;
time(1:10000,1:n)=0;
% Vx
load vxmdl2.mat
tic
for i=1:n
    j=runs(i);
    Q=Q_arr(j);
    rpm=rpm_arr(j);
    Q_=repmat(Q,[181 1]);
    rpm_=repmat(rpm, [181 1]);
    %vx pred
    X_valid = [z Q_.*1e7 rpm_];
    predvx = predict(mdl, X_valid);

    % D
    v=mean(predvx);
    predD=(-6.7387e-06) - (33.7226.*Q) - (4.0753e-07.*rpm)+(0.025795.*v) + ...
        (Q-4.786e-7).*((v-0.0021)*(-1.7119e+04)) + (Q - 4.786e-07).*((Q - 4.786e-07).*1.2983e+08);
    vx=predvx;
    D=predD;

    y0=[D; 1.4492];

%% 
    time_exp=RTD(:,2*j-1);
    time_exp=time_exp(~isnan(time_exp));
    Et_exp=RTD(:,2*j);
    Et_exp=Et_exp(~isnan(Et_exp));
    yopt_=y0;%[y0(i) y0(end)];
    [pred] = CD_Model(yopt_,vx,vz,Q,time_exp,Et_exp,'run');
    E_model(1:length(pred),i)=pred(:,1);
    time(1:length(pred),i)=pred(:,2);
    % plot(time_exp,Et_exp,'r',time,E_model,'--r','Linewidth',2)
    % hold on;

    MRT_model(i)=pred(1);
    MRT_exp(i)=pred(2);
    vari_model(i)=pred(3);
    vari_exp(i)=pred(4);
    skew_model(i)=pred(5);
    skew_exp(i)=pred(6);

end
toc

%% 

MRT_m_avg=mean(MRT_model);
MRT_e_avg=mean(MRT_exp);
vari_m_avg=mean(vari_model);
vari_e_avg=mean(vari_exp);
%%
figure
plot(RTD(:,11),RTD(:,12),'b',time(1:7000,1),E_model(1:7000,1),'--b','Linewidth',2); hold on
plot(RTD(:,17),RTD(:,18),'r',time(1:3000,2),E_model(1:3000,2),'--r','Linewidth',2);hold on
plot(RTD(:,21),RTD(:,22),'g',time(1:3000,3),E_model(1:3000,3),'--g','Linewidth',2)
legend({'Run 6-CFD','Run 6-ROM','Run9-CFD','Run9-ROM','Run11-CFD','Run11-ROM'},'FontSize',12)
xlim([0 250])
xlabel('Time, s','FontSize',20)
ylabel('E(t)','FontSize',20)
%title('RTD Comparison for Validation Runs','FontSize',20)
set(gcf,'color','white')

%%
% train_runs=[1 2 3 4 5 7 10 12 13 14 15 16 17];
% test_runs=[6 9 11];
% ParityPlot2(MRT_model(train_runs)',MRT_exp(train_runs)',MRT_model(test_runs)',MRT_exp(test_runs)')
% title('MRT Parity Plot','fontsize',15) 
% 
% ParityPlot2(vari_model(train_runs)',vari_exp(train_runs)',vari_model(test_runs)',vari_exp(test_runs)')
% title('Variance Parity Plot','fontsize',15) 


% ParityPlot2(skew_model(train_runs),skew_exp(train_runs),skew_model(test_runs),skew_exp(test_runs))
% title('Skewness','fontsize',15)



