
clear
seeds = (1:25)+804;
trials = length(seeds)
p = 104;
poi = zeros(trials,p); 
real_saved = zeros(trials,4);
real_saved_neut = zeros(trials,4);
poi_neut = zeros(trials,p); 
lasso_coef = zeros(trials, p);
misclas = zeros(trials,1);
misclas_neut = zeros(trials,1);
misclas_lasso = zeros(trials,1);
 
for  s=1:trials
    rng(seeds(s)*100) 
    
    run('simdataright.m')
    run('InputSimGBM.m')
    
    real_saved(s,:) = real;
    poi(s,:) = MargGam;
    misclas(s) = numError; 
    
    rng(seeds(s)*100)  
    run('simdataneutral.m')
    run('InputSimGBM.m')
    
    real_saved_neut(s,:) = real;
    poi_neut(s,:) = MargGam;
    misclas_neut(s) = numError; 

     [B fitInfo] = lasso(X,Y,'CV',10);
     lambda1SE = fitInfo.Lambda1SE;
     idxLambda1SE = fitInfo.Index1SE;
 
     coef = B(:,idxLambda1SE);
     lasso_coef(s,:) = coef; 
     coef0 = fitInfo.Intercept(idxLambda1SE);
 
     yhat = Xf*coef + coef0;
     predYf = yhat > 0.5;
     misclas_lasso(s) = 25- sum(predYf == Yf);
end

% csvwrite('poiJan9.csv', poi)
% csvwrite('poi_neutJan9.csv', poi_neut)
% csvwrite('realJan9.csv', real_saved)
% csvwrite('real_neutJan9.csv', real_saved_neut)
% csvwrite('misclasJan9.csv', misclas)
% csvwrite('misclas_neutJan9.csv', misclas_neut) 
% 

    % for i = 1:3
    %     find(poi(i,:) >= 0.1) 
    %     find(poi_neut(i,:) >= 0.1)
    % end

% x = 1:1:104;
% y = MargGam;
% 
% X = [x; x];
% Y = [repelem(0,104); y.'];
% figure; hold on;
% plot(x,y,'*'),  xlim([1 104]),ylim([0 1.05])
% title('Var Selection');
% line(X,Y); 
% hold off;



[misclas misclas_neut misclas_lasso ]

tpr = sum(sum(poi(:,1:4) > 0.5)) / 4*trials
fpr = sum(sum(poi(:,5:104) > 0.5)) / 100*trials

tpr_neut = sum(sum(poi_neut(:,1:4) > 0.5)) / 4*trials
fpr_neut = sum(sum(poi_neut(:,5:104) > 0.5)) / 100*trials

tpr_lasso = sum(sum(abs(lasso_coef(:,1:4)) > 1e-6)) / 4*trials
fpr_lasso = sum(sum(abs(lasso_coef(:,5:104)) > 1e-6)) /100*trials

%% matthew's correlation coeff

tp = sum(sum(poi(:,1:4) > 0.5))
fp = sum(sum(poi(:,5:104) > 0.5))
tn = 100*trials - fp
fn = 4*trials - tp

mcc = (tp * tn - fp * fn) / (sqrt((tp+fp)*(tp + fn)*(tn+fp)*(tn + fn)))

tp = sum(sum(poi_neut(:,1:4) > 0.5))
fp = sum(sum(poi_neut(:,5:104) > 0.5))
tn = 100*trials - fp
fn = 4*trials - tp

mcc_neut = (tp * tn - fp * fn) / (sqrt((tp+fp)*(tp + fn)*(tn+fp)*(tn + fn)))

tp = sum(sum(abs(lasso_coef(:,1:4)) > 1e-6))
fp = sum(sum(abs(lasso_coef(:,5:104)) > 1e-60))
tn = 100*trials - fp
fn = 4*trials - tp

mcc_lasso = (tp * tn - fp * fn) / (sqrt((tp+fp)*(tp + fn)*(tn+fp)*(tn + fn)))

acc = (25*trials - sum(misclas)) / (25*trials)
acc_neut = (25*trials - sum(misclas_neut)) / (25*trials)
acc_lasso = (25*trials - sum(misclas_lasso)) / (25*trials)

%% look at ROI of LASSO

lambda = [10:-0.05:0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00001, 0];
[B fitInfo] = lasso(X,Y,'CV',10, 'Lambda', lambda);

n_lambda = size(lambda, 2);
tpr_roc_lasso = zeros(n_lambda, 1);
fpr_roc_lasso = zeros(n_lambda, 1);

for beta_ind = 1:n_lambda 
  var_sel = abs(B(:, beta_ind)) > 1e-6;
  tpr_roc_lasso(beta_ind) = sum(var_sel(1:4)) / 4;
  fpr_roc_lasso(beta_ind) = sum(var_sel(5:104)) / 100;
end


      
% Fit lasso with 10-fold CV
% Default setting for alpha is 1 (i.e. lasso)
lasso_opts = glmnetSet;
lambda = [10:-0.05:0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00001, 0];
lasso_cv = cvglmnet(X, Y,  'gaussian', 10, [], 'response', lasso_opts, 0);
 
% Results for all parameters tried
beta_sel = lasso_cv.glmnet_object.beta;
n_lambda = size(lambda, 2);
tpr_roc_lasso = zeros(n_lambda, 1);
fpr_roc_lasso = zeros(n_lambda, 1);
for beta_ind = 1:n_lambda 
  var_sel = beta_sel(:, beta_ind) ~= 0;
  [tpr_roc_lasso(beta_ind), fpr_roc_lasso(beta_ind)] = tpr_fpr_var(var_sel, gamma_true);
end
 
csvwrite(strcat('tpr_fpr_roc_lasso_model', num2str(model), '_iter', ...
              num2str(cur_iter), '.csv'), [tpr_roc_lasso, fpr_roc_lasso]);
