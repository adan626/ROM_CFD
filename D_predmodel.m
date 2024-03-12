clear all; clc
data = readmatrix('D_pred.csv');
Q = data(:,1);
rpm = data(:,2);
D = data(:,4);
v = data(:,3);

%%
[sf, g] = fit([Q, rpm],D,'poly33')
figure
plot(sf,[Q, rpm],D)

%%
X = [Q rpm v];
y=D;
coefficients = mvregress(X,y)

% Predictions
y_pred = X * coefficients;

% Calculate the Mean Squared Error
mse = mean((y - y_pred).^2);
fprintf('Mean Squared Error: %f\n', mse);

% Calculate R-squared
SST = sum((y - mean(y)).^2); % Total sum of squares
SSE = sum((y - y_pred).^2);  % Sum of squared errors
R2 = 1 - SSE / SST;         % Coefficient of determination (R-squared)

fprintf('R-squared (R2): %f\n', R2);

% plot([Q rpm],D)
%%
figure
scatter3(Q,rpm,D,'filled')

%%
Q = data(:,1);
rpm = data(:,2);
D = data(:,4);
v = data(:,3);

% Q=Q([6 9 11]);
% rpm=rpm([6 9 11]);
% v=v([6 9 11]);

% Q=Q([1 2 3 4 5 7 8 10 12 13 14 15 16 17]);
% rpm=rpm([1 2 3 4 5 7 8 10 12 13 14 15 16 17]);
% v=v([1 2 3 4 5 7 8 10 12 13 14 15 16 17]);
%predD=(-6.328e-6) + (Q - 4.7857e-7) .* ((Q - 4.7857e-7).*88353568.5)...
  %  -27.89.*Q - 3.515e-7.*rpm + 0.0227686342974021.*v;


predD=(-6.7387e-06) - (33.7226.*Q) - (4.0753e-07.*rpm)+(0.025795.*v) + ...
    (Q-4.786e-7).*((v-0.0021)*(-1.7119e+04)) + (Q - 4.786e-07).*((Q - 4.786e-07).*1.2983e+08);
y=D;%([6 9 11]);
y_pred=predD;
% Calculate the Mean Squared Error
mse = mean( ((y - y_pred).^2)./(y.^2));
rmse=sqrt(mse);
fprintf('Mean Squared Error: %f\n', rmse);
%%
% Calculate R-squared


SST = sum((y - mean(y)).^2); % Total sum of squares
SSE = sum((y - y_pred).^2);  % Sum of squared errors
R2 = 1 - SSE / SST;         % Coefficient of determination (R-squared)

fprintf('R-squared (R2): %f\n', R2);

nrmse = sqrt(mean((y - y_pred).^2))/(max(y) - min(y));
disp(['Mean Squared Error: ' num2str(nrmse)]);

trainy_pred = y_pred([1 2 3 4 5 7 8 10 12 13 14 15 16 17]);
trainy=y([1 2 3 4 5 7 8 10 12 13 14 15 16 17]);
testy_pred = y_pred([6 9 11]);
testy=y([6 9 11]);

ParityPlot2(trainy_pred,trainy,testy_pred,testy)
title('Diffusion Coefficient','fontsize',20) 

%%
y=testy;
y_pred=testy_pred;
nrmse = sqrt(mean((y - y_pred).^2))/(max(y) - min(y));
disp(['Mean Squared Error: ' num2str(nrmse)]);

