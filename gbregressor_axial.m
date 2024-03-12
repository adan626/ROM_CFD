%update file path

clear all; clc;

%filepath = '\axial_velocity.csv';
df = readtable('axial_velocity.csv');




%validation data
valid1 = df((df.Q == 3) & (df.RPM == 55), :);
valid2 = df((df.Q == 3) & (df.RPM == 75), :);
valid3 = df((df.Q == 4) & (df.RPM == 45), :);

valid_data = [valid1; valid2; valid3];

%training data
exclude_conditions = (df.Q == 3 & df.RPM == 55)|(df.Q == 3 & df.RPM == 75)|(df.Q == 4 & df.RPM == 45);
train_data = df(~exclude_conditions, :);

%separation of fearures and target variable for training and validation
X_train = train_data(:, {'z','Q', 'RPM'});
y_train = train_data.axial_velocity;

X_valid = valid_data(:, {'z','Q', 'RPM'});
y_valid = valid_data.axial_velocity;

%print shape of training and validation data
disp(['X_train: ' num2str(size(X_train)) ', y_train: ' num2str(size(y_train)) ',' ...
    'X_valid: ' num2str(size(X_valid)) ', y_valid: ' num2str(size(y_valid))]);
%model definition
tic
mdl = fitrensemble(X_train, y_train, 'Method', 'LSBoost', 'NumLearningCycles', 200, 'LearnRate', 0.05, ...
    'Learners', 'tree');
d=toc;
%predictions
train_pred = predict(mdl, X_train);
y_pred = predict(mdl, X_valid);

%print Normalized Root Mean Squared Error(NRMSE) for training and validation
nrmse = sqrt(mean((y_train - train_pred).^2))/(max(y_train) - min(y_train));
disp(['Mean Squared Error for Training: ' num2str(nrmse)]);

nrmse = sqrt(mean((y_valid - y_pred).^2))/(max(y_valid) - min(y_valid));
disp(['Mean Squared Error for Validation: ' num2str(nrmse)]);

%print R-squared for training and validation
r2_train = 1 - sum((y_train - train_pred).^2) / sum((y_train - mean(y_train)).^2);
disp(['R-squared for Training: ' num2str(r2_train)]);

r2_valid = 1 - sum((y_valid - y_pred).^2) / sum((y_valid - mean(y_valid)).^2);
disp(['R-squared for Validation: ' num2str(r2_valid)]);
%%

valid = df((df.Q == 4) & (df.RPM == 45), :);
X_valid = valid(:, {'z','Q', 'RPM'});

y_pred = predict(mdl, X_valid);
mean2=mean(y_pred)
mean1=mean(valid.axial_velocity)



figure
plot(valid.z,valid.axial_velocity,valid.z,y_pred)