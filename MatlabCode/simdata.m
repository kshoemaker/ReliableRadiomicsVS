function [X, Xf, Y, Yf, N] = ... 
    simdata(n, balance, p1, p2, mu1, informed)

% n is sample size
% balance is proportion of sample in group 1
% p1 is the number of discriminating features
% p2 is the number of noise features
% mu1 is the vector of true group mean values, 1 value per group
% informed - informed N or not?, 1 for informed, 0 for neutral

% setup params
% n = 100; 
n1 = floor(balance*n); 
n2 = n-n1;
% p1 = 4; p2 = 100; % num of "real" variables and noise ones

Y = [repelem(1,n1) repelem(2,n2)].';
% mu1 = [-.8,.8];

sigma = repmat(.2,p1,p1)+diag(repmat(0.8,1,p1));

sigma_noise = 0.3*ones(p2,p2) + 0.4*eye(p2);

% create true features
x1 = mvnrnd(repelem(mu1(1),p1), sigma, n1);
x2 = mvnrnd(repelem(mu1(2),p1), sigma, n2);


dat = [x1;x2];

% create noise features
mu_noise = repelem(0, p2);

dat_noise =   mvnrnd(mu_noise, sigma_noise, n);


Data = [dat, dat_noise];
% diag(cov(Data));

% new_order = randperm(p1+p2);
new_order = 1:(p1+p2);
X = Data;
% X_mean = mean(X);
% X = bsxfun(@minus, X, X_mean);
real = 1:p1;
false = (p1+1):(p1+p2);

if informed == 1
    N = zeros(1,p1+p2)+0.15;
    N(real) = 0.35;
elseif informed == 0
    N = zeros(1,p1+p2)+0.15;
end 

out = randperm(n);
% 75/25 test train split
train_size = floor(n*0.75);
% test_size = n - train_size;
train = out(1:train_size);
test = out((train_size+1):n);

%X = X-repmat(mean(X),size(X,1),1);
Xf = X(test,:);
Xf_mean = mean(Xf);
Xf = bsxfun(@minus, Xf, Xf_mean); % centering at 0
X = X(train,:);
X_mean = mean(X);
X = bsxfun(@minus, X, X_mean);

Yf = Y(test)-1; % because we set y as 1 and 2
Y = Y(train)-1;


