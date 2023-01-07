
cd '/Users/katherineshoemaker/Documents/GitHub/ReliableRadiomicsVS/SimulationOutputFiles/informative'

cd '/Users/katherineshoemaker/Documents/GitHub/ReliableRadiomicsVS/SimulationOutputFiles/noninformative'


d = dir('*.mat')

numfiles = length(d);

misclass_vec = zeros(1, numfiles);
misclass_lasso_vec = zeros(1, numfiles);
misclass_svm_vec = zeros(1, numfiles);
tpr_class_vec = zeros(1, numfiles);
tpr_class_lasso_vec = zeros(1,numfiles);
fpr_class_vec = zeros(1, numfiles);
fpr_class_lasso_vec = zeros(1, numfiles);

for i = 1:numfiles
    load(d(i).name);
    misclass_vec(i) = misclas;
    misclass_lasso_vec(i) = misclas_lasso;
    misclass_svm_vec(i) = misclas_svm;
    tpr_class_vec(i) = tpr_class;
    tpr_class_lasso_vec(i) =  tpr_class_lasso;
    fpr_class_vec(i) = fpr_class;
    fpr_class_lasso_vec(i) = fpr_class_lasso;
end 

disp(['Misclass Average: ' num2str(mean(misclass_vec))])
disp(['Misclass Average - LASSO: ' num2str(mean(misclass_lasso_vec))])
disp(['Class loss estimate - SVM: ' num2str(mean(misclass_svm_vec))])
disp(['TPR - Average: ' num2str(mean(tpr_class_vec))])
disp(['TPR LASSO - Average: ' num2str(mean(tpr_class_lasso_vec))])
disp(['FPR - Average: ' num2str(mean(fpr_class_vec))])
disp(['FPR LASSO - Average: ' num2str(mean(fpr_class_lasso_vec))])


