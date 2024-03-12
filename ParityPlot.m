function [ ] = ParityPlot(x1,y)

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
plot(x1,y,'ks', 'LineWidth', 3, 'MarkerSize', 9);  % observed vals 
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
leg=legend('Reference line','','15% error interval','Observed vals');
set(leg,'location','northwest')
xlabel('Model', 'FontSize', fontSize);
ylabel('Experiment', 'FontSize', fontSize);
set(gcf,'color','white')
bias=sum(y-x1)/length(y);
tbl = table(y , x1);
mdl = fitlm(tbl,'linear');

%---compute R2 and RMSE
e=abs(y-x1);
ybar=mean(y);
SST=sumsqr(y-ybar);
SSR=sumsqr(e);
R_sq=1-SSR/SST;
Rootmean=rmse(x1,y);

mdl = fitlm(x1,y);
R_sq2=mdl.Rsquared.ordinary;
 

str=[    'N = ',sprintf('%d',mdl.NumObservations),...
%', Bias = ',sprintf('%.3f',bias),...    
', R^2 = ',sprintf('%.3f',R_sq2),... 
 ', RMSE = ',sprintf('%.2f',Rootmean),];

annotation('textbox',[.35 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black') 
% str=[ ' Train R^2 = ',sprintf('%.3f',R_sq2)];
%  annotation('textbox',[.45 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black')  

end