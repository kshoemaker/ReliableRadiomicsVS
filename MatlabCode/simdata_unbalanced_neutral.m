
% clear 
% datasetseed = 1023;
% rng(datasetseed);

% setup params
n = 100; n1 = 60; n2 = n-n1;
p1 = 4; p2 = 100; % num of "real" variables and noise ones
Y = [repelem(1,n1) repelem(2,n2)].';
mu1 = [-0.8,0.8];
%beta = repelem(0.8,2);
sigma = repmat(.2,p1,p1)+diag(repmat(0.8,1,p1));

sigma_noise = 0.3*ones(p2,p2) + 0.4*eye(p2);


% make random normal(0,1) data for p1 features on n1 (n2) subjects
% x1 = randn(n1,p1);
% x2 = randn(n2,p1);
% S = chol(sigma);

% just use the mvnrnd command now that I have it! 
x1 = mvnrnd(repelem(mu1(1),p1), sigma, n1);
x2 = mvnrnd(repelem(mu1(2),p1), sigma, n2);

%dat_1 = x1*S + ones(n1,p1)*mu1(1);% + beta(1).*SNPsig(1:n1,1)*ones(1,p1) +beta(2).*SNPsig(1:n1,2)*ones(1,p1);  
%dat_2 = x2*S + ones(n2,p1)*mu1(2);% +%beta(1).*SNPsig(n1+1:n1+n2,1)*ones(1,p1) +beta(2).*SNPsig(n1+1:n1+n2,2)*ones(1,p1);  

%dat = [dat_1;dat_2];
dat = [x1;x2];

mu_noise = repelem(0, p2);
%samples = randn(n, p2);
%A = chol(sigma_noise);
%dat_noise = samples*A;

dat_noise =   mvnrnd(mu_noise, sigma_noise, n);

Data = [dat, dat_noise];
diag(cov(Data));

% new_order = randperm(p1+p2);
new_order = 1:(p1+p2);
X = Data(:,new_order);
% X_mean = mean(X);
% X = bsxfun(@minus, X, X_mean);
real = find(new_order<(p1+1));
false = find(new_order>(p1+1));
N = zeros(1,p1+p2)+0.15;
% N(real) = 0.25;
%N(false(1:4)) = 0.5;


out = randperm(n);
train = out(1:75);
test = out(76:n);

%X = X-repmat(mean(X),size(X,1),1);
Xf = X(test,:);
Xf_mean = mean(Xf);
Xf = bsxfun(@minus, Xf, Xf_mean);
X = X(train,:);
X_mean = mean(X);
X = bsxfun(@minus, X, X_mean);

Yf = Y(test)-1;
Y = Y(train)-1;

Y1 = find(Y==0); Y2 = find(Y==1); Y1f = find(Yf==0); Y2f = find(Yf==1);
n1 = sum(Y==0); n2 = sum(Y==1);
G=max(Y)+1;

