%% demo  HN_PA  dataset %%%%%%%%%%%%%
clc
clear all

% Input data
load 'HN_PA'   %Load data file
load 'HN_PASamplecategory' %Load the sample category file
X=HN_PA; %data name HN_PA
X=mapminmax(X,0,1);  %Data Preparation
%parameter setting
lambda1=10^ -4;
lambda2=10^1;
lambda3=10^-3;
epsilon1 = 1e-6;
epsilon2 = 1e-2; 
maxIter = 50;
% Obtaining low-rank symmetric matrix.
[Z_hat, ~, ~, ~, ~] =NSLRG(X, lambda1, lambda2, lambda3, epsilon1, epsilon2, maxIter);
% low-rank symmetric matrix matrix is used as input of Score function to select important features.
[gene_select_result, gene_number, NSLRG_desccending_order,select_number,der2_value]=NSLRG_S(X,Z_hat);
% Kmeans is taken to cluster samples by using selected feature genes.
kmeans_number=50;
[total_label,total_res,result ] =K_means_and_Measurement_Metrics  (Samplecategory,kmeans_number,gene_select_result);

save result.mat result NSLRG_desccending_order select_number 
beep;