%  %%%%% First step of the model code %%%%%
 clear
 addpath('../CaseStudyData')
 importCaseStudy('HNModelData_ALL_Oct20HPV.mat')
% 
  
   N = zeros(1,162) +median(N); 
    N = N/3;
  % data = [X; Xf]; 
% 
% % centering the case study data, together! 
%  data = normalize(data);
% % note, 9/29/2020: the  data is centered and scaled in R
% 
%   rng(2389);
%   order = randperm(101); 
%   X = data(order(1:101), :);
%   Xf = data(order(76:101), :);
% 
% dataY = [Y; Yf];
% Y = dataY(order(1:75));
% Yf = dataY(order(76:101));


%%%%% Setting up the Parameters and Hyperparameters %%%%%

p = size(X, 2); % number of radiomic features

%MCMC setting
gam_prior = []; % selection latent variablel

n_iter =  100000; % desired number of MCMC iterations
r1 = 2; % number of starting variables (gamma)
mu_j1 = zeros(p, 1); % class specific mu_jk
mu_j2 = zeros(p, 1); % class specific mu_jk

% Hyperparameters setting
% probability of Add/Delete vs Swap
pPar = 0.5; % as recommended in the motif code

% Prior and hyperpriors on gamma
 alpha_0 = -2.75;
 alpha_1 = 3;

% KS: alpha_1 is set later,  where used. Should move alpha_0 there too. 

% These aren't used since we fixed alpha_1
% w = 0.05;   % mean of alpha_1 normal distrub % changed from 2 to 0,  6.14.18
% t = .05;    % removed them 

% Prior on sigma0j
a = 3;  
b = 0.1;

% Prior on sigmaj1 and sigmaj2
ak = 3; 
bk = 0.1;

% Prior on Sigma1 and Sigma2
c = 0.6;   % doesn't get passed to the function,  used for computing bj 
Q = c*eye(p); %% standard setting for the IW, times c
d_k = 3;      %% standard setting for the IW

% that implies for the Inv-Gamma Prior on sigmaj's
aj = d_k/2; 
bj = 2/c; 

% Prior on mu0's
h1 = 1; % scale parameter


%% Run 3 chains 
%%%%%% Calls the "main program" and runs the model %%%%%%

% Run the mainprog - Chain #1
[mu_1_mat, mu_2_mat, log_prob, ADgam, Sgam, ...
   acceptgam, acceptmu1, acceptmu2, GammaA]= ... 
mainprogGBM(X, Y, gam_prior, n_iter, r1, mu_j1, mu_j2, pPar, alpha_0, alpha_1, ...
    a, b, ak, bk, Q, d_k, aj, bj, N, h1);

% set the burn-in 
bi = 20000;

% save the non-burnin sections
GammaBI1 = GammaA((bi+1):end);
log_prob1 = log_prob((bi+1):end);
mu_1_mat1 = mu_1_mat(:, (bi+1):end);
mu_2_mat1 = mu_2_mat(:, (bi+1):end);

%Beta01_PostMean_1 = Beta01_PostMean; Beta02_PostMean_1 = Beta02_PostMean;
%MargGam_1 = MargGam; MargDel_11 = MargDel1; MargDel_12 = MargDel2;

clear GammaA  log_prob mu_1_mat mu_2_mat;
%save filename_1.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the mainprog - Chain #2
[mu_1_mat, mu_2_mat, log_prob, ADgam, Sgam, ...
   acceptgam, acceptmu1, acceptmu2, GammaA]= ... 
mainprogGBM(X, Y, gam_prior, n_iter, r1, mu_j1, mu_j2, pPar, alpha_0, alpha_1, a, b, ak, bk, Q, d_k, aj, bj, N, h1);


GammaBI2 = GammaA((bi+1):end);
log_prob2 = log_prob((bi+1):end);
mu_1_mat2 = mu_1_mat(:, (bi+1):end);
mu_2_mat2 = mu_2_mat(:, (bi+1):end);

clear GammaA  log_prob mu_1_mat mu_2_mat;
%save filename_2.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the mainprog - Chain #3
[mu_1_mat, mu_2_mat, log_prob, ADgam, Sgam, ...
   acceptgam, acceptmu1, acceptmu2, GammaA]= ... 
mainprogGBM(X, Y, gam_prior, n_iter, r1, mu_j1, mu_j2, pPar, alpha_0, alpha_1, a, b, ak, bk, Q, d_k, aj, bj, N, h1);


GammaBI3 = GammaA((bi+1):end);
log_prob3 = log_prob((bi+1):end);
mu_1_mat3 = mu_1_mat(:, (bi+1):end);
mu_2_mat3 = mu_2_mat(:, (bi+1):end);

clear GammaA  log_prob mu_1_mat mu_2_mat;
%save filename_3.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pool together 3 MCMC Chains

GammaA = [GammaBI1 ; GammaBI2; GammaBI3];
log_prob = [log_prob1,  log_prob2,  log_prob3];
mu_1_mat = [mu_1_mat1 mu_1_mat2 mu_1_mat3];
mu_2_mat = [mu_2_mat1 mu_2_mat2 mu_2_mat3];


disp(' ')
disp('------- Calculating marginal probabilities for gamma (ROIs)')
disp(' ')
sss = size(GammaA); 
GammaBI=GammaA; % (bi:sss(1)); 
aa = size(GammaBI);
div = 1; 
bb = round(aa(1)/div);
aaa = 1; 
freTot = zeros(p, 1);

for j=1:div
    bbb = bb*j;
    if j==div 
        bbb=size(GammaBI);
    end
    FreqUni = [];
    for i=aaa:bbb
        FreqUni = [FreqUni str2num(GammaBI{i})];
    end
    fre = tabulate(FreqUni); fsize = size(fre);
    m = p - fsize(1); Z0 = zeros(m, 1);
    fre = [fre(:,2);Z0]; freTot = fre+freTot;
    aaa = bbb+1;
end
disp(' ')
disp('------- Selected ROIs ...')
disp(' ')
MargGam = freTot./aa(1); 	
feature_thresh = 0.5;
find(MargGam>feature_thresh)


    ROI_sel = find(MargGam>feature_thresh); %SNP1_sel = find(MargDel1>0.5); SNP2_sel = find(MargDel2>0.5);
    % biB = 100000;
    Beta01_PostMean = mean(mu_1_mat(ROI_sel, :), 2);
    Beta02_PostMean = mean(mu_2_mat(ROI_sel, :), 2);

    test = isempty(ROI_sel);
    if test == 0
        Beta01_gam = []; Beta02_gam = [];
        sss = size(GammaA)-1; 
        GammaBI=GammaA; %(bi:sss(1)); 
        aa = size(GammaBI);
        div = 1; 
        bb = round(aa(1)/div); 
        aaa = 1; 
        freTot = zeros(p, 1);
        bbb=size(GammaBI);
        Beta01_gamBI = mu_1_mat(ROI_sel, :);
        Beta02_gamBI = mu_2_mat(ROI_sel, :);
        for i=aaa:bbb
            if length(str2num(GammaBI{i}))==length(ROI_sel)
                if str2num(GammaBI{i})==ROI_sel'
                    Beta01_gam = [Beta01_gam; Beta01_gamBI(:, i)'];
                    Beta02_gam = [Beta02_gam; Beta02_gamBI(:, i)'];
                end
            end    
        end
        disp(' ')
        disp('------- Posterior means (group 1) ...')
        disp(' ')
        mean(Beta01_gam)
        disp(' ')
        disp('------- Posterior means (group 2) ...')
        disp(' ')
        mean(Beta02_gam)
    end

    selected = size(ROI_sel);

[obvs, vars]  = size(X);
     plotx = 1:1:vars;
     ploty = MargGam;
%     
%     
     plotX = [plotx; plotx];
     plotY = [repelem(0, vars); ploty.'];
     figure; hold on;
     plot(plotx, ploty, '*'),   xlim([1 vars]), ylim([0 1.05])
     %title(mytitle);
     line(plotX, plotY,  'color', 'r'); 
       % line(plotX(:, real), plotY(:, real),  'color', 'b');
%     
     hold off;
 
%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%    Prediction 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(' ')
    disp('------- Prediction  -------')
    disp(' ')
    threshold = 0.5; %  feature_thresh;
    numError = zeros(1, length(threshold)); 
    numVar = zeros(1, length(threshold));

    biB = bi;
    Yf = double(Yf);
    Y1 = find(Y==0); 
    Y2 = find(Y==1); 
    Y1f = find(Yf==0); 
    Y2f = find(Yf==1);
    n1 = sum(Y==0); 
    n2 = sum(Y==1);
    n = n1 + n2;
    G=max(Y)+1;

    % a_k^prime
    ak_p1 = ak + n1/2; %% error found 7/28, for n1 =!= n2
    ak_p2 = ak + n2/2;
    
    % a_k^new
    akf1 = ak_p1 + 1/2;
    akf2 = ak_p2 + 1/2;    

    numErrorV = [];
    numVarV = []; 


        gammaf = MargGam>threshold;
        gammaList = find(gammaf);
        Xgam = X(:, logical(gammaf));
        XgamC = X(:, logical(1-gammaf));
        Xfgam = Xf(:, logical(gammaf));
        XfgamC = Xf(:, logical(1-gammaf));
        Nk =[];
        X1gam = Xgam(Y1, :);
        X2gam = Xgam(Y2, :); 
        ROI_sel = gammaList;
        % Beta01_PostMean = mean(mu_1_mat(ROI_sel, biB:end), 2);
        % Beta02_PostMean = mean(mu_2_mat(ROI_sel, biB:end), 2);
        % changing so that we only use the model average from the models
        % that select the 4
        % changing back so we have means.
        
        % only when the selectected are selected
        beta01f =  mean(Beta01_gam); % Beta01_PostMean;
        beta02f =  mean(Beta02_gam);
        
        % all values 
        beta01f =  Beta01_PostMean;
        beta02f =  Beta02_PostMean;
        [nf,  pf] = size(Xf);
        ClassP = [];
        ClassPM = [];
        likelis = [];
        PpostM = [];
        PostProb = [];
        [nselected , nothing] = size(gammaList);
        % nselected = selected(1);


     for i=1:nf        
            AA = 1; 
            BB = 1; 
            b1_p = 1; 
            b2_p = 1;
            b1f = 1; 
            b2f = 1;

            test = isempty(ROI_sel);
            if test == 0
                for j=1:nselected
                    Xj1 = X1gam(:, j);
                    Xj1f = Xfgam(i, j);

                    M = (Xj1-repmat(beta01f(j)', n1, 1));  
                    %b_k^prime
                    b1_p = b1_p * (bk +  ((Xj1-repmat(beta01f(j)', n1, 1))'*(Xj1-repmat(beta01f(j)', n1, 1))/2));                                
                    
                    %b_k^new
                    b1f = b1f * (bk + ((Xj1f-beta01f(j)')'*(Xj1f-beta01f(j)')+M'*M)/2);


                    Xj2 = X2gam(:, j);
                    Xj2f = Xfgam(i, j);

                    M = (Xj2-repmat(beta02f(j)', n2, 1)); 
                    %b_k^prime
                    b2_p = b2_p * (bk + ((Xj2-repmat(beta02f(j)', n2, 1))'*(Xj2-repmat(beta02f(j)', n2, 1))/2));
                    
                    %b_k^new
                    b2f = b2f * (bk + ((Xj2f-beta02f(j)')'*(Xj2f-beta02f(j)')+M'*M)/2);   
                    

                end
            end 
            % b1_p = bk + b1_p/2;
            % b2_p = bk + b2_p/2;
            % b1f = bk + b1f/2;
            % b2f = bk + b2f/2;
            AA =  exp(ak_p1 * log(b1_p) - akf1 * log(b1f)) ;% (b1_p^ak_p1)/(b1f^akf1);  
            BB =  exp(ak_p2 * log(b2_p) - akf2 * log(b2f)) ;    % (b2_p^ak_p2)/(b2f^akf2);

            likf1 = (-(1/2)*log((2*pi)) + log(n1/n) + log(AA) + gammaln(akf1) - gammaln(ak_p1))
            likf2 = (-(1/2)*log((2*pi)) + log(n2/n) + log(BB) + gammaln(akf2) - gammaln(ak_p2))
            %likf1 = 1/(2*pi)^(n1/2)*(n1/n)*AA*gamma(akf1)/gamma(ak_p1); 
            %likf2 = 1/(2*pi)^(n2/2)*(n2/n)*BB*gamma(akf2)/gamma(ak_p2);

            Ppost = [likf1 likf2];
            likelis = [likelis ; Ppost] ;% storing the computed likelihoods of class 1 and 2
            ExpPpost = [exp(likf1) exp(likf2)];
            PpostM =  ExpPpost./sum(ExpPpost); 
            PostProb = [PostProb ; PpostM ]; % added to track the class probabilities
            ClassP = [ClassP find(Ppost==max(Ppost))]; % assigns class based on max likelihood
                                                       % is "correct" -
                                                       % matches Yf
            ClassPM = [ClassPM find(PpostM==max(PpostM))]; % assigns class based on class probability
                                                           % is not correct
                                                           % - opposite of
                                                           % Yf
    end           
            err = sum((ClassP-1-(Yf)')~=0);
            numErrorV = [numErrorV err ];
            numVarV = [ numVarV sum(gammaf)];


        numError(1, :) = numErrorV;
        numVar(1, :) = numVarV ;

        display(numErrorV)
        display(numVarV)


    disp(' ')
    disp('------- number of miss-classified subjects...')
    disp(' ')
    display(numError)

    disp(' ')
    disp('------- Miss-classification rate...')
    disp(' ')
    display(numError/length(Yf))

% [likelis PostProb (ClassP-1)'  Yf]
tpr = sum((ClassP-1) == 1 & Yf' == 1) / sum(Yf' == 1) ;
fpr = sum((ClassP-1) == 1 & Yf' == 0) / sum(Yf' == 0) ;

% %% trace plotting 
% 
% % % code for plotting trace of variables selected
% % countofvars  = cellfun(@(x) length(str2num(x)), GammaBI);
% % plot(1:73334 , countofvars(1:73334), 1:73334, countofvars(73335:146668), 1:73334, countofvars(146669:220002))
% % 
% % % code for plotting trace of first mean parameter
% % plot(1:93333, Beta01_gamBI(1, 1:93333), 1:93333, Beta01_gamBI(1, 93334:186666), 1:93333, Beta01_gamBI(1, 186667:279999))
% % plot(1:93333, Beta02_gamBI(1, 1:93333), 1:93333, Beta02_gamBI(1, 93334:186666), 1:93333, Beta02_gamBI(1, 186667:279999))
% 
% %% LASSO
% 
[B fitInfo] = lasso(X,Y,'CV',10);
  lambda1SE = fitInfo.Lambda1SE;
 idxLambda1SE = fitInfo.Index1SE;
%  
 [B,FitInfo] = lassoglm(X,Y,'binomial','CV',3);
 idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
 B0 = FitInfo.Intercept(idxLambdaMinDeviance);
 coef = [B0; B(:,idxLambdaMinDeviance)] ;
%  
% 
% coef = B(:,idxLambda1SE);
%lasso_coef = coef; 
 %coef0 = fitInfo.Intercept(idxLambda1SE);
%  
      yhat2 = exp(Xf*B(:,idxLambdaMinDeviance) + B0);
      yhat = glmval(coef,Xf,'logit');
      predYf = yhat > 0.5;
      misclas_lasso = length(Yf) - sum(predYf == Yf)
% %% ROI Curve
% % code for ROI curve
% % 
% % thresh_vec = 0:0.001:1;
% % sens = sum(MargGam(1:4) > thresh_vec)/4;
% % oneminusspec = 1 - sum(MargGam(5:104) < thresh_vec) / 100;
% % plot( oneminusspec, sens, '-o')
% % xlim([0 1])
% % ylim([0 1])
% 
% %% gelman rubin?
% % Done with 
% % clear 
% % datasetseed = 102389;
% % rng(datasetseed);
% 
% % addpath('/Users/katherineshoemaker/Documents/CompStatsToolboxV2/')
% % open '/Users/katherineshoemaker/Documents/CompStatsToolboxV2/Contents'
% % 
% % nu1 = [mean(mu_1_mat1'); mean(mu_1_mat2'); mean(mu_1_mat3')];
% % nu2 = [mean(mu_2_mat1'); mean(mu_2_mat2'); mean(mu_2_mat3')];
% % 
% % GR1 = csgelrub(nu1(:,1:4))
% % GR2 = csgelrub(nu2)
% 
% % save('predictingHPVwprior.mat')
% 
% confusionchart((Yf)', ClassP-1)
% 
% sum(N(ROI_sel) > median(N))
% mean(N(setdiff(1:length(N), ROI_sel)))
% 
% 
% %%%%%%%%%%%%%%%
% % do the SVM
% 
% train_data = [X Y];
% test_data = [Xf Yf];
% 
% svm_test = fitcsvm(X,Y);
% 
% sv = svm_test.SupportVectors;
% figure
% gscatter(X(:,1),X(:,2),Y)
% hold on
% plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
% legend('not HPV','HPV','Support Vector')
% hold off
% 
% 
% 
% svm_predict = predict(svm_test, Xf); 
% confusionchart(Yf, svm_predict)
