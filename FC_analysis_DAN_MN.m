%% Functional connectivity analysis between networks
%  model inspired by Tavor et al. to predict inter-individual differences 
% in motor tasks exploiting MEG-based resting-state FC matrices. The idea
% is to consider links between two different networks to capture the
% inter-individual differences in connectivity patterns.
%
% Author: Anna Vettoruzzo
% Septemeber 2021


clear all;
close all;
clc;

%% Select band, ntw, task
% the idea is to first consider the entire matrix (164x164), 
% apply degree-normalization and finally restrict the analysis 
% only to the portion of the matrix of interest

datadir = ['D:\Tesi\Data\'];

band = 'betahigh'; %betahigh, betalow
ntw = 'AllNtw';
ntw1 = 'Dorsal Attention Network';
ntw2 = 'Motor Network';

task = 2; %1 = LH, 2 = LF, 4 = RH, 5 = RF, 6 = meanH, 7 = meanF 
task_name = {'LH', 'LF', '', 'RH', 'RF', 'meanH', 'meanF'};

subjects=textread([datadir, '\subjects.txt'],'%s');

fprintf('NETWORKS: %s - %s\n', ntw1, ntw2);
fprintf('BAND: %s\n', band);
fprintf('TASK: %s\n\r', task_name{task});


%% LOAD MATRICES
%FC_matrices is an array #subjects x 5. The columns represent:
%the name of the subjects (1), the rest matrices (2, 3, 4) and the task one (5).

fprintf('LOAD FC MATRICES\n\r');
FC_matrices = cell(length(subjects),5);
for i=1:length(subjects)
    FC_matrices{i,1} = subjects{i};
    
    %REST
    fname = [datadir, subjects{i}, '\Results\FunctionalConnectivity\Static\Rest\3-Restin\', band, '\', subjects{i}, '_MEG_3-Restin_icablpcorr_', band, 'conn.mat'];
    if isfile(fname)
        m = load(fname);  
        %find nodes belonging to the network selected
        if ntw == "AllNtw"
            nodes = [1:size(m.conn.complete,1)].';
        else
            nodes = m.conn.P2N{(find(strcmp(m.conn.P2N, ntw))),3};
        end
        FC_matrices{i,2} = m.conn.complete(nodes,nodes);
    else
        fprintf(['Subjects ', subjects{i}, ' without REST1 matrix\n']); 
        FC_matrices{i,2} = [];
    end
    
    fname = [datadir, subjects{i}, '\Results\FunctionalConnectivity\Static\Rest\4-Restin\', band, '\', subjects{i}, '_MEG_4-Restin_icablpcorr_', band, 'conn.mat'];
    if isfile(fname)
        m = load(fname);
        FC_matrices{i,3} = m.conn.complete(nodes,nodes);
    else
        fprintf(['Subjects ', subjects{i}, ' without REST2 matrix\n']); 
        FC_matrices{i,3} = [];
    end
    
    fname = [datadir, subjects{i}, '\Results\FunctionalConnectivity\Static\Rest\5-Restin\', band, '\', subjects{i}, '_MEG_5-Restin_icablpcorr_', band, 'conn.mat'];
    if isfile(fname)
        m = load(fname);
        FC_matrices{i,4} = (m.conn.complete(nodes,nodes));
    else
        fprintf(['Subjects ', subjects{i}, ' without REST3 matrix\n']); 
        FC_matrices{i,4}=[];
    end
    
    %TASK
    if task == 1 %LH
        fname_task = [datadir, subjects{i}, '\Results\FunctionalConnectivity\Static\Task\class1\', band, '\', subjects{i}, '_MEG_Motort-10-11(class1)_icablpcorr_', band, 'conn.mat'];
        m = load(fname_task);
        FC_matrices{i,5} = m.conn.complete(nodes,nodes); 
    elseif task == 2 %LF
        fname_task = [datadir, subjects{i}, '\Results\FunctionalConnectivity\Static\Task\class2\', band, '\', subjects{i}, '_MEG_Motort-10-11(class2)_icablpcorr_', band, 'conn.mat'];
        m = load(fname_task);
        FC_matrices{i,5} = m.conn.complete(nodes,nodes);    
    elseif task == 4 %RH
        fname_task = [datadir, subjects{i}, '\Results\FunctionalConnectivity\Static\Task\class4\', band, '\', subjects{i}, '_MEG_Motort-10-11(class4)_icablpcorr_', band, 'conn.mat'];    
        m = load(fname_task);
        FC_matrices{i,5} = m.conn.complete(nodes,nodes);     
    elseif task == 5 %RF
        fname_task = [datadir, subjects{i}, '\Results\FunctionalConnectivity\Static\Task\class5\', band, '\', subjects{i}, '_MEG_Motort-10-11(class5)_icablpcorr_', band, 'conn.mat'];    
        m = load(fname_task);
        FC_matrices{i,5} = m.conn.complete(nodes,nodes);     
    elseif task == 6 %meanH
        fname_task = [datadir, subjects{i}, '\Results\FunctionalConnectivity\Static\Task\class1\', band, '\', subjects{i}, '_MEG_Motort-10-11(class1)_icablpcorr_', band, 'conn.mat'];
        m1 = load(fname_task);

        fname_task = [datadir, subjects{i}, '\Results\FunctionalConnectivity\Static\Task\class4\', band, '\', subjects{i}, '_MEG_Motort-10-11(class4)_icablpcorr_', band, 'conn.mat'];
        m2 = load(fname_task);

    	FC_matrices{i,5} = (m1.conn.complete(nodes,nodes)+m2.conn.complete(nodes,nodes))/2; 
    elseif task == 7 %meanF
        fname_task = [datadir, subjects{i}, '\Results\FunctionalConnectivity\Static\Task\class2\', band, '\', subjects{i}, '_MEG_Motort-10-11(class2)_icablpcorr_', band, 'conn.mat'];
        m1 = load(fname_task);

        fname_task = [datadir, subjects{i}, '\Results\FunctionalConnectivity\Static\Task\class5\', band, '\', subjects{i}, '_MEG_Motort-10-11(class5)_icablpcorr_', band, 'conn.mat'];
        m2 = load(fname_task);

    	FC_matrices{i,5} = (m1.conn.complete(nodes,nodes)+m2.conn.complete(nodes,nodes))/2;  
    end   
    
end


%arrays only with links between the two networks
nodes1 = m.conn.P2N{(find(strcmp( m.conn.P2N, ntw1))),3};
nodes2 = m.conn.P2N{(find(strcmp( m.conn.P2N, ntw2))),3};
FC_arrays = cell(length(subjects),size(FC_matrices,2));
for i = 1:length(subjects)
    subj = subjects{i};
    FC_arrays{i,1} = subj;
    
    for j=2:size(FC_matrices,2)
        if size(FC_matrices{i,j}) == [0,0]
            FC_arrays{i,j} = [];
        else
            FC_arrays{i,j} = reshape(FC_matrices{i,j}(nodes1,nodes2).',[],1);
        end
    end
end

ntwName1 = strtrim(m.conn.P2N{(find(strcmp( m.conn.P2N, ntw1))),2});
ntwName2 = strtrim(m.conn.P2N{(find(strcmp( m.conn.P2N, ntw2))),2});
ntwName = [ntwName1 '-' ntwName2];
%% CONSISTENCY between rest matrices (using cosine similarity measure)
%before degree normalization
fprintf('CONSISTENCY BETWEEN REST MATRICES\n');
fprintf('Cosine Similarity\n');
CS=zeros(length(subjects),3);

%rest1 vs. rest2
for i = 1:length(subjects)
    a = FC_arrays{i,2};
    b = FC_arrays{i,3};
    CS(i,1) = cosine_similarity(a,b);       
end

%rest1 vs. rest3
for i = 1:length(subjects)
    a = FC_arrays{i,2};
    b = FC_arrays{i,4};
    CS(i,2) = cosine_similarity(a,b);
end

%rest2 vs. rest3
for i = 1:length(subjects)
    a = FC_arrays{i,3};
    b = FC_arrays{i,4};
    CS(i,3) = cosine_similarity(a,b);     
end

CS_tot = mean(CS, 2, 'omitnan');
for i=1:length(subjects)
    fprintf('Consistency subj %s', subjects{i});
    fprintf(': %3.3f\n', CS_tot(i));
end
avgCS = mean(CS_tot);
fprintf('Average Consistency: %3.3f\n\r', avgCS);



%Consistency between rest matrices with added WHITE NOISE
%add white noise with mean 5 and var 15 to a random rest matrix
%mantain the same position of NaN values
rnd_idx = randsample(3,1)+1;
FC_arrays_wnoise = FC_arrays;
for i = 1:size(FC_arrays,1)   
    FC_arrays_wnoise{i,rnd_idx} = FC_arrays{i,rnd_idx} + 15*randn(size(FC_arrays{i,rnd_idx}))+5;
end

fprintf('Cosine Similarity with Noisy Data\n');
CS_wnoise=zeros(length(subjects),3);

%rest1 vs. rest2
for i = 1:length(subjects)
    a = FC_arrays_wnoise{i,2};
    b = FC_arrays_wnoise{i,3};
    CS_wnoise(i,1) = cosine_similarity(a,b);       
end

%rest1 vs. rest3
for i = 1:length(subjects)
    a = FC_arrays_wnoise{i,2};
    b = FC_arrays_wnoise{i,4};
    CS_wnoise(i,2) = cosine_similarity(a,b);
end

%rest2 vs. rest3
for i = 1:length(subjects)
    a = FC_arrays_wnoise{i,3};
    b = FC_arrays_wnoise{i,4};
    CS_wnoise(i,3) = cosine_similarity(a,b);     
end

CS_tot_wnoise = mean(CS_wnoise, 2, 'omitnan');
avgCS_wnoise = mean(CS_tot_wnoise);
fprintf('Average Consistency with Noise: %3.3f\n\r', avgCS_wnoise);

%% APPLY DEGREE-NORMALIZATION and save only the upper triangular part

fprintf('DEGREE-NORMALIZATION\n\r');

for i=1:length(subjects)  
    for j=2:size(FC_matrices,2)
        sub_m = degreeNorm(FC_matrices{i,j}); %degree-norm applied to matrix with only some nodes
        FC_matrices{i,j} = triu(sub_m);
    end
end

%% Consider only links between the two networks

for i=1:length(subjects)
    for j = 2:size(FC_matrices,2)
        if size(FC_matrices{i,j}) == [0,0]
            FC_matrices{i,j} = [];
            FC_arrays{i,j} = [];
        else
            FC_matrices{i,j} = FC_matrices{i,j}(nodes1,nodes2);
            FC_arrays{i,j} = reshape(FC_matrices{i,j}.',[],1);
        end
    end
end
%% Verify position of NaN values

%verify if subjects have NaN in the same positions in all matrices
fprintf('POSITION OF NaN \n')
for i=1:length(subjects)
    disp(FC_matrices{i,1});
    k = 2;
    if size(FC_matrices{i,k}) == [0,0]
        k = k+1;
    end
    [r,c] = find(isnan(FC_matrices{i,k}));
    for j=(k+1):size(FC_matrices, 2)
        [r_tmp,c_tmp] = find(isnan(FC_matrices{i,j}));
        if ~(all(r_tmp == r) && all(c_tmp == c))
            disp('NaN in different positions'); 
        end
    end
    disp('OK');
end

%how many NaN in each subject
nNaN = zeros(length(subjects),1);
for i=1:length(subjects)
    nNaN(i) = sum(sum(isnan(FC_matrices{i,5})));
end

%verify structure similarities between matrices
%NB: it is the same for all bands and for all tasks.
similarity_matrix = zeros(length(subjects), length(subjects));
for i = 1:length(subjects)
    mi = FC_arrays{i,5};
    mi(~isnan(mi)) = 1;
    mi(isnan(mi)) = 0;
    for j = 1:length(subjects)
        mj = FC_arrays{j,5};
        mj(~isnan(mj)) = 1;
        mj(isnan(mj)) = 0;
        
        mij = mi==mj;
        similarity_matrix(i,j) = sum(mij);
    end
end

s_matrix = triu(similarity_matrix);
s_matrix(s_matrix==0) = NaN;
s_matrix = (s_matrix / length(FC_arrays{1,2})) * 100;

%show figure
figure('name', 'Similarity Matrix');
h = heatmap(s_matrix);
color = generatecolormapthreshold([70 80 90 95 99 100],[0.5430 0 0; 1 0 0; 1 0.6 0; 1 0.8 0; 1 1 1]);
h.Colormap = color;
h.ColorLimits = [70 100];
h.FontSize = 7;
h.Title = ['Similarity matrix - ',ntwName];
h.XLabel = 'Subjects';
h.YLabel = 'Subjects';

%savefigure
folderName = ['Results/', ntwName, '/'];
if ~exist(folderName, 'dir')
   mkdir(folderName)
end
figName = ['similarityMatrix_',ntwName];
saveas(h, [folderName , figName, '.png']);

fprintf('Minimum value of similarity: %3.3f\n\r', min(s_matrix,[],'all'));

%% Regression models 
%a model for each subject and each link

%compute beta coefficients for each subjects
fprintf('START LEARNING\n\r');
betas = cell(length(subjects), 1);
for i=1:length(subjects)
   betas{i,1} = subjects{i};
   x = mean(cat(2, FC_arrays{i,2}, FC_arrays{i,3}, FC_arrays{i,4}),2); %REST matrices
   
   y = FC_arrays{i,5}; %TASK matrix (can be left, right or mean)
   
   %consider only FC values not NaN
   ind = find(all(~isnan(x),2)); 
   ind_nan = find(all(isnan(x),2)); 
    
   %one regression model for each link
   beta = zeros(length(x),1);
   r = zeros(length(x),1);
   for j = 1:length(ind)
       idx = ind(j);
       [beta(idx),~,r(idx)] = regress(y(idx),x(idx));
   end
   beta(ind_nan) = NaN;
   betas{i,2} = beta;  
end


fprintf('LEAVE-ONE-OUT PREDICTION\n\r');
FC_pred = cell(length(subjects),2);
for i = 1:length(subjects)
    LOOsubject = subjects{i};
    FC_pred{i} = LOOsubject;
    
    %compute avg beta 
    loo_betas = zeros(size(betas{i,2},1),length(subjects)-1);
    cnt = 1;
    for j=setdiff(1:length(subjects),i)
        loo_betas(:,cnt) = betas{j,2};
        cnt = cnt+1;
    end
    avg_beta = mean(loo_betas,2);   
    
    %prediction
    pred= [];
    x = mean(cat(2, FC_arrays{i,2}, FC_arrays{i,3}, FC_arrays{i,4}),2);

    pred = x.*avg_beta;
    FC_pred{i, 2} = pred;
end    

%% Save true and predicted elements both as arrays and matrices
%save FC_pred as a matrix
n1 = length(nodes1);
n2 = length(nodes2);
for i=1:length(subjects)
    FC_pred{i,3} = reshape(FC_pred{i,2}, n1, n2);
end

%save FC_true as a matrix
FC_true = cell(length(subjects),3);
for i=1:length(subjects)
    FC_true{i,1} = subjects{i};
    FC_true{i,2} = FC_arrays{i,5};
    FC_true{i,3} = FC_matrices{i,5};
end

%% Goodness of fit REST-TASK model
%MSE between FC_pred and FC_true for each subject

mse = zeros(length(subjects),1);
diff_nan = zeros(length(subjects),1);
for i = 1:length(subjects)
    pred = FC_pred{i,2};
    ind_pred = find(all(~isnan(pred),2));
    true = FC_true{i,2};
    ind_true = find(all(~isnan(true),2));
    diff_nan(i) = abs(length(ind_pred)-length(ind_true));
    
    c = intersect(ind_pred, ind_true); %find common values in both arrays
    mse(i) = immse(pred(c),true(c));
end

avg_mse = mean(mse);
std_mse = std(mse);

fprintf('PERFORMANCE REST-TASK\n');
fprintf('Average MSE: (%3.3f ', avg_mse*10^5);
fprintf('± %3.3f)e-05 \n\n', std_mse*10^5);

%% Correlation (with Pearson) - Identifiability matrix

fprintf('IDENTIFIABILITY MATRIX\n');
IDmatrix = zeros(length(subjects));
for i = 1:length(subjects)
    for j = 1:length(subjects)
        IDmatrix(i,j) = corr(FC_pred{i,2},FC_true{j,2}, 'rows', 'complete');
    end         
end

%show figure (elements with negative correlation are shown in white)
figName = [task_name{task}, '_', band, '_', ntwName];
figure('name', figName);
tmp_IDmatrix = IDmatrix;
h = heatmap(tmp_IDmatrix);
h.FontSize = 7;
%h.MissingDataColor = [1,1,1];
h.caxis([-0.2,1]); %to set the same scale
cmap = jet(20);
cmap = cmap(8:19,:);
h.Colormap = cmap;
h.Title = ['IDmatrix', ' - ',task_name{task},' ',band,' ',ntwName];
h.xlabel('Subjects (true)');
h.ylabel('Subjects (predicted)');
%save figure
saveas(h, [folderName , figName, '.png']);

fprintf('Minimum value in the identifiability matrix: %3.3f\n\r', min(tmp_IDmatrix,[],'all'));

%% Performance of the model: 

fprintf('PERFORMANCE OF THE MODEL\n');
%PREDICTION SCORE (ref. Tavor et al.)
mean_diag = mean(diag(IDmatrix));
std_diag = std(diag(IDmatrix));
fprintf('Prediction score : %3.3f ', mean_diag);
fprintf('± %3.3f \n\n', std_diag);

%IDENTIFICATION RATE (ref. Finn et al.)
count = 0;
for i=1:length(subjects)
    corr = IDmatrix(i,:);
    [m, idx] = max(corr);
    if idx == i
        count = count+1;
    end
end
IDrate = count/length(subjects)*100;
fprintf('Identification Rate: %3.3f\n', IDrate);


%% Write in a file
fileID = fopen([folderName ,figName,'.txt'],'w');
%cosine similarity
fprintf(fileID, 'CS for each subject: \n');
fmt = [repmat('%3.3f   ', 1, length(subjects)), '%3.3f \n'];
fprintf(fileID, fmt, CS_tot);
fprintf(fileID, '\nAverage CS: %3.3f \n', avgCS);
%Structure similarity
fprintf(fileID, 'Min value of similarity: %3.3f \n\n', min(s_matrix,[],'all'));
%MSE
fprintf(fileID, 'MSE for each subject: \n');
fmt = [repmat('%3.3f e-05   ', 1, length(subjects)), '%3.3f \n'];
fprintf(fileID, fmt, mse*10^5);
fprintf(fileID, '\nAverage MSE: (%3.3f ', avg_mse*10^5);
fprintf(fileID, '± %3.3f)e-05 \n\n', std_mse*10^5);
%Pscore
fprintf(fileID, 'Prediction score : %3.3f ', mean_diag);
fprintf(fileID, '± %3.3f \n\n', std_diag);
%IDrate
fprintf(fileID, 'Identification rate: %3.3f%%\n\n', IDrate);
%Min value IDmatrix
fprintf(fileID, 'Min value in IDmatrix: %3.3f\n', min(tmp_IDmatrix,[],'all'));
fclose(fileID);

%% Normalization: degree-normalization (ref. Chiem)
%give in input sub_m=matrix with only selected nodes
function Mnorm = degreeNorm(M)
    M(isnan(M)) = 0;
    M = abs(M);
    d = sum(M, 2,'omitnan');    %degree of each single node
    D = diag(d);    %create a diagonal matrix
    %the inverse of D is obtained changing every value in the diagonal with
    %its reciprocal
    Mnorm = (D^(-0.5)) * M * (D^(-0.5));
    Mnorm(Mnorm == 0) = NaN;
end

%% Cosine similarity (ref. Menon2019)
%is given by the inner product of the two vectors divided by their norms.
function y = cosine_similarity(a,b)
    a(isnan(a)) = 0;    %change NaN entries with 0
    b(isnan(b)) = 0;
    c=0;					
    for k=1:length(a)		
        c = c + a(k)*b(k);	% update c by the k-th product in inner product
    end			   
    a_norm = norm(a);
    b_norm = norm(b);
    y = c/(a_norm*b_norm);
end