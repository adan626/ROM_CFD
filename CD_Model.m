function [output] = CD_Model(yopt,vx1,vz,Q,time_exp,Et_exp,runtype)
%define variables
xmin=0;
zmin=0;
xmax=0.27;    %length of extruder in m % 0.27 
zmax=0.025;   % width of the extruder in m 
L=xmax;  %extruder length 
W=zmax;  %extruder width 

factor=181;

%1,1% of the extruder is the pulse injection  (1.1%) 
dx=(1.1/(factor*1.1))*xmax;   % axial length of bin
dz=zmax/1;  % to create 10 bins in the radial direction 
dt=0.05;       % time step for integration     

t=0;
if (runtype=='opt')
    tmax=350;
else 
    tmax=time_exp(end);
end

%discretize the domain (linear)
x = xmin-dx : dx: 1.1*xmax;  %axial 
z = zmin:dz:zmax;   %radial 

N=201;
vx1(length(vx1)+1:N)=vx1(end)*3;%*1.5;

vx=vx1*yopt(end);%yopt(end);

%D=1*0.20*vx./80;  %diffusion constants   %optimize this value 
D=ones(length(vx),1)*yopt(1);

CFL=max(vx)*dt/dx+max(vz)*dt/dz;  % CFL condition should be less than 1. 

% initial conditions
u0=zeros(length(x),length(z));
min_num=find(x==0); % find the position where x=0. initial pulse 
max_num=find(x>0.269 & x<0.271); % find the position where x=0.27 (exit)
%u0(min_num,:)=10/length(z); %initialize with position. 
u0(min_num,1)=100/2; %initialize with position. 
u0(min_num,2)=100/2; %initialize with position. 
maxu0=sum(sum(u0));
u=u0;

%compute average residence time 
frac=1-0.6681;
vol=L*2*pi*(zmax/2)^2;
tau=frac*vol/Q;
 
% check if pulse is added to 1.1% of the extruder 
check_percentage=x(min_num+1)/x(max_num)*100;

% loop through time 
nsteps=round(tmax/dt);
a=1; 

mean_vx1=mean(vx1);

%vx=vx';

for n=1:nsteps
    time(n)=t; 
    unp1=u;    

    Cend_num(n)=sum(u(max_num,:));  %RTD C versus t 

    % integration + numerical diff
    %forward in time, centered in space for diffusion
    %and backward in space for advection

             unp1(2:length(x)-1,2:length(z))=u(2:length(x)-1,2:length(z))+D(1:length(x)-2).*u(1:length(x)-2,2:length(z)).*dt./(dx.^2)+D(3:length(x)).*u(3:length(x),2:length(z)).*dt./(dx.^2)-2.*D(2:length(x)-1).*u(2:length(x)-1,2:length(z)).*dt./(dx.^2)...
                 -(dt/dx).*[u(2:length(x)-1,2:length(z)).*vx(2:length(x)-1)-u(1:length(x)-2,2:length(z)).*vx(1:length(x)-2)]...
                 -vz(2:length(z)).*(dt/dz).*u(2:length(x)-1,2:length(z))+vz(1:length(z)-1).*(dt/dz).*u(2:length(x)-1,1:length(z)-1);
   

             unp1(2:length(x)-1,1)=u(2:length(x)-1,1)-vz(1).*(dt/dz).*u(2:length(x)-1,1)+vz(end)*(dt/dz).*u(2:length(x)-1,end)...
                 +D(1:length(x)-2).*u(1:length(x)-2,1).*dt/(dx.^2)+D(3:length(x)).*u(3:length(x),1).*dt/(dx.^2)-2.*D(2:length(x)-1).*u(2:length(x)-1,1).*dt/(dx.^2)...
                 -(dt/dx).*[u(2:length(x)-1,1).*vx(2:length(x)-1)-u(1:length(x)-2,1).*vx(1:length(x)-2)];
  
    %boundary conditions in the ith direction 
    unp1(end,2:length(z))=u(end,2:length(z))+vx(end-1)*(dt/dx)*u(end-1,2:length(z))-D(end)*(dt/dx.^2)*u(end,2:length(z))+D(end-1)*(dt/dx.^2)*u(end-1,2:length(z))...
    +vz(1:length(z)-1)*(dt/dz)*u(end,1:length(z)-1)-vz(2:length(z))*(dt/dz)*u(end,2:length(z)); 
 
    unp1(1,2:length(z))=u(1,2:length(z))-vx(1)*(dt/dx)*u(1,2:length(z))+D(2)*(dt/(dx.^2))*u(2,2:length(z))-D(1)*(dt/dx.^2)*u(1,2:length(z))...
        -vz(2:length(z))*(dt/dz)*u(1,2:length(z))+vz(1:length(z)-1)*(dt/dz)*u(1,1:length(z)-1); 
 
    %boundary conditions at 1,1 position 
    unp1(1,1)=u(1,1)-vx(1)*(dt/dx)*u(1,1)+D(2)*(dt/(dx.^2))*u(2,1)-D(1)*(dt/dx.^2)*u(1,1)... 
           -vz(1)*(dt/dz)*u(1,1)+vz(end)*(dt/dz)*u(1,end); 

    %boundary conditions at end,1 position 
    unp1(end,1)=u(end,1)+vx(end-1)*(dt/dx)*u(end-1,1)+D(end-1)*(dt/(dx.^2))*u(end-1,1)-D(end)*(dt/dx.^2)*u(end,1)...  
           -vz(1)*(dt/dz)*u(end,1)+vz(end)*(dt/dz)*u(end,end); 
     u=unp1; 
     t=t+dt;
     check=sum(sum(u)); % must add up to initial pulse conc. 
  
end 

%compute E(t) curve
Cend_num_bound=Cend_num(end)*1.1;
Cend_num(Cend_num<Cend_num_bound)=0;
den=trapz(time,Cend_num);
E=Cend_num./den;


%compute metrics for model
num=trapz(time,time.*Cend_num);
den=trapz(time,Cend_num);
MRT_model=num/den;
vari_model=trapz(time,((time-MRT_model).^2).*E);
skew_model=1/(sqrt(vari_model)^1.5)*trapz(time,((time-MRT_model).^3).*E);


%compute metrics for exp
if (runtype=='fit')|(runtype=='run')

%resample experimental RTD graph
time_mod=0:dt:tmax;
Et_mod=spline(time_exp,Et_exp,time_mod);
MRT_exp=trapz(time_exp,time_exp.*Et_exp);
error_MRT=abs(MRT_model-MRT_exp)/max(MRT_model,MRT_exp);

vari_exp=trapz(time_mod,((time_mod-MRT_exp).^2).*Et_mod);
error_vari=abs(vari_model-vari_exp)/max(vari_model,vari_exp);

skew_exp=1/(sqrt(vari_exp)^1.5)*trapz(time_mod,((time_mod-MRT_exp).^3).*Et_mod);
error_skew=abs(skew_model-skew_exp)/max(skew_model,skew_exp);

%compute error for corr factor
[M,I]=max(Et_mod);
max_time_exp=time_mod(I);

[M_model,I_model]=max(E);
max_time_model = time(I_model);

err_time=abs(max_time_model-max_time_exp)/max(max_time_exp,max_time_model);

% 
% figure
% plot(time,E,'--',time_exp,Et_exp,'Linewidth',1)
% legend original resampled 
 
RTD_error=norm(Et_mod(1:length(E))-E,1);
end

Pe=L*mean(vx)/mean(D);   %compute Peclet number

if (runtype == 'fit')
    output = err_time+error_MRT+error_vari;%+error_skew;%(err_time+error_vari+error_skew)*100;%error_vari+error_MRT;%RTD_error;
elseif (runtype =='run')
    %output = [MRT_model,MRT_exp,vari_model,vari_exp,skew_model,skew_exp];
    output = [E' time'];
elseif (runtype == 'opt')
    output = [MRT_model,vari_model,skew_model,Pe];
end


end