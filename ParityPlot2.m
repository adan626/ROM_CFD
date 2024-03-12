function [ ] = ParityPlot2(x1,y,x1valid,yvalid)

figure
maxi=1.5*max(max(x1,y));
mini=0.5*min(max(x1,y));
xx=[mini maxi];
yy1=xx*0.85;
yy2=xx*1.15;
yy=[yy1;(yy2-yy1)]';

xlim([mini maxi])
ylim([mini maxi])
plot(xlim,ylim,'--k','LineWidth',1.25)   % reference line 

hold on
ha=area(xx,yy);
set(ha(1), 'FaceColor', 'none') % this makes the bottom area invisible
set(ha(2),'FaceColor','k')
set(ha, 'LineStyle', 'none','FaceAlpha',0.15)
axis tight
hold on

X = [ones(size(x1)) x1];
[b,bint] = regress(y,X) ;
xval = min(x1)-0.05:0.01:max(x1)+0.05;
yhat = b(1)+b(2)*xval;
ylow = bint(1,1)+bint(2,1)*xval;
yupp = bint(1,2)+bint(2,2)*xval;
plot(y,x1,'ks', 'LineWidth', 3, 'MarkerSize', 9);  % observed vals 
plot(yvalid,x1valid,'rs', 'LineWidth', 3, 'MarkerSize', 9);  % observed vals 

p5.Color(4) = 0.5;
hold on;
%p6=plot(xval,yhat,'k','linewidth',3);
%p6.Color(4) = 0.5;
%axis([0.04 0.3 0.03 .35])

fontSize = 20;
hold on
 

alpha(0.3);
xlim([mini maxi])
ylim([mini maxi])
leg=legend('Reference line','','15% error interval','Train','Valid');
set(leg,'location','northwest')
xlabel('Actual', 'FontSize', fontSize);
ylabel('Prediction', 'FontSize', fontSize);
set(gcf,'color','white')
bias=sum(y-x1)/length(y);
biasv=sum(yvalid-x1valid)/length(yvalid);

tbl = table(y , x1);
mdl = fitlm(tbl,'linear');

%---compute R2 and RMSE
e=abs(y-x1);
ybar=mean(y);
SST=sumsqr(y-ybar);
SSR=sumsqr(e);
R_sq=1-SSR/SST;
Rootmean=rmse(x1,y);
nrmse1=sqrt(1/length(y)*sum(abs((x1valid-yvalid)./x1valid).^2));

mdl = fitlm(x1,y);
R_sq2=mdl.Rsquared.ordinary;

evalid=abs(yvalid-x1valid);
ybar_v=mean(yvalid);
SST=sumsqr(yvalid-ybar_v);
SSR=sumsqr(evalid);
R_sq_v=1-SSR/SST;
Rootmeanv=rmse(x1valid,yvalid);
nrmsev=sqrt(1/length(y)*sum(abs((x1-y)./x1).^2));

 
mdlv = fitlm(x1valid,yvalid);
R_sq2_v=mdlv.Rsquared.ordinary;
str=[    'Train: N = ',sprintf('%d',mdl.NumObservations),...
', R^2 = ',sprintf('%.3f',R_sq2),... 
 ', NRMSE = ',sprintf('%.2f',nrmse1),];

strv=[    'Valid: N = ',sprintf('%d',mdlv.NumObservations),...
', R^2 = ',sprintf('%.3f',R_sq2_v),... 
 ', NRMSE = ',sprintf('%.2f',nrmsev),];
% 
% strv=['Train: R^2 = ',sprintf('%.3f',R_sq2_v)];
% str=['Valid: R^2 = ',sprintf('%.3f',R_sq2)];
annotation('textbox',[.35 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black') 
annotation('textbox',[.35 0.83 0 0],'string',strv,'FitBoxToText','on','EdgeColor','black') 

% str=[ ' Train R^2 = ',sprintf('%.3f',R_sq2)];
%  annotation('textbox',[.45 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black')  

end