clear;
%% compile different approaches
mex QRCDM_cversion.cpp
mex QRCDM_AP_cversion.cpp
mex PDHG_QDSFM_cversion.cpp
mex Subgradient_QDSFM_cversion.cpp

%% generate data
% number of vertices
N = 1000;
% for each vertex, din many hyperedges within the same cluster are generated (typically <=1) 
din = 1; 
% for each vertex, dout many hyperedges within the same cluster are generated (typically <=1) 
dout = 1;
% size of hyperedge
K = 20;
% rate of observed label 
observed_rate= 0.002;
a = zeros(1,N);
% assign observed labels
a(randperm(N/2, N/2*observed_rate)) = 1;
a(randperm(N/2, N/2*observed_rate)+N/2) = -1;
% count the number of hyperedges
count = 0;
% incidence_list{i} contains the vertices that hyperedge i covers
incidence_list = {};
% parameter_homo_list{i} contains the weight of hyperedge i
parameter_homo_list = {};
% submodular_type{i} contains the submodular type of the hyperedge i:
% standard hypergraphs use 'h'
submodular_type = {};
% generate hyperedges within the first cluster
for k = 1:N/2,
    if rand(1)> din,
        continue;
    end
    count = count + 1;
    list = randperm(N/2, K-1);
    while sum(list==k)>0,
        list = randperm(N/2, K-1);
    end
    incidence_list{count} = [k list];
    parameter_homo_list{count} = 1;
    submodular_type{count} = 'h';
end
% generate hyperedges within the second cluster
for k = 1:N/2,
    if rand(1)> din,
        continue;
    end
    count = count + 1;
    list = randperm(N/2, K-1);
    while sum(list==k)>0,
        list = randperm(N/2, K-1);
    end
    incidence_list{count} = [k list]+N/2;
    parameter_homo_list{count} = 1;
    submodular_type{count} = 'h';
end
% generate hyperedges across the two clusters
for k = 1:(N*dout),
    count = count + 1;
    list = randperm(N, K);
    while min(list) > N/2 || max(list) <= N/2,
        list = randperm(N, K);
    end
    incidence_list{count} = list;
    parameter_homo_list{count} = 1;
    submodular_type{count} = 'h';
end
incidence_list = incidence_list(:);
parameter_homo_list = parameter_homo_list(:);
submodular_type = submodular_type(:);

% total number of hyperedges
R =count;

% degrees of different vertices
degree_vec = degree_stat_homo(incidence_list, parameter_homo_list, N, R);
bias_vec = a./degree_vec;

%%semisupervised learning
% For all QDSFM problems, different approaches are to solve

% min_x ||x-bias_vec||_W^2 + sum_r [f_r(x)]^2

%lambda for QDSFM problems 
lambda_QDSFM = 0.02;
%weighted matrix for norm 
W = lambda_QDSFM*degree_vec;
%number of iterations
T = 300;
record_dis = T/30;
tic;
[x_PDHG, record, final_gap] = PDHG_QDSFM_cversion(incidence_list, parameter_homo_list, submodular_type, bias_vec, W, N, R, T,record_dis);
time_PDHG = toc;
[clustering_err, conductance, thre] = result_analysis_homo(incidence_list, parameter_homo_list, x_PDHG, degree_vec, N, R);
gap_PDHG = record;
fprintf('PDHG-based QDSFM result:\n #incorrect clustered vertices:%d\n conductance:%f\n cputime:%f\n', clustering_err, conductance, time_PDHG);

%weighted matrix for norm 
W = lambda_QDSFM*degree_vec;
%number of iterations
T = 300*R;
record_dis = T/30;
tic;
[x_QRCDM, record, final_gap] = QRCDM_cversion(incidence_list, parameter_homo_list, submodular_type, bias_vec, W, N, R, T,record_dis);
time_QRCDM = toc;
[clustering_err, conductance, thre] = result_analysis_homo(incidence_list, parameter_homo_list, x_QRCDM, degree_vec, N, R);
fprintf('RCD-based QDSFM result:\n #incorrect clustered vertices:%d\n conductance:%f\n cputime:%f\n', clustering_err, conductance, time_QRCDM);

%weighted matrix for norm
W = lambda_QDSFM*degree_vec;
%number of iterations
T = 300;
record_dis = T/30;
tic;
[x_AP, record, final_gap] = QRCDM_AP_cversion(incidence_list, parameter_homo_list, submodular_type, bias_vec, W, N, R, T,record_dis);
time_AP = toc;
[clustering_err, conductance, thre] = result_analysis_homo(incidence_list, parameter_homo_list, x_AP, degree_vec, N, R);
fprintf('AP-based QDSFM result:\n #incorrect clustered vertices:%d\n conductance:%f\n cputime:%f\n', clustering_err, conductance, time_AP);

%weighted matrix for norm
W = lambda_QDSFM*degree_vec;
%number of iterations
T = 15000;
record_dis = T/30;
tic;
[x_Subgradient, record] = Subgradient_QDSFM_cversion(incidence_list, parameter_homo_list, submodular_type, bias_vec, W, N, R, T,record_dis);
time_Subgradient = toc;
[clustering_err, conductance, thre] = result_analysis_homo(incidence_list, parameter_homo_list, x_Subgradient, degree_vec, N, R);
fprintf('Subgradint-based QDSFM result:\n #incorrect clustered vertices:%d\n conductance:%f\n cputime:%f\n', clustering_err, conductance, time_Subgradient);

%lambda for clique-expansion method 
lambda_CE = 0.001;
%weighted matrix for norm
W = lambda_CE*degree_vec;
tic;
x_CE = clique_expansion(incidence_list, parameter_homo_list, a, lambda_CE, N, R);
time_CE = toc;
[clustering_err, conduction, thre] = result_analysis_homo(incidence_list, parameter_homo_list, x_CE, degree_vec, N, R);
fprintf('clique-expansion result:\n #incorrect clustered vertices:%d\n conductance:%f\n cputime:%f\n', clustering_err, conductance, time_CE);






