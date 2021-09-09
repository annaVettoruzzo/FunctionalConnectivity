%% Functional connectivity analysis
%  model inspired by Tavor et al. to predict inter-individual differences 
% in motor tasks exploiting MEG-based resting-state FC matrices.
% Identify how the performance of the model changes depending on the size
% of the dataset.
%
% Author: Anna Vettoruzzo
% Septemeber 2021

clear all;
close all;
clc;

%% Select band, ntw, task

datadir = ['D:\Tesi\Data\'];

band = 'betalow'; %betahigh, betalow
ntw = 'Motor Network'; %'AllNtw' o altre reti

task = 7; %1 = LH, 2 = LF, 4 = RH, 5 = RF, 6 = meanH, 7 = meanF 
task_name = {'LH', 'LF', '', 'RH', 'RF', 'meanH', 'meanF'};

totSubjects=textread([datadir, '\subjects.txt'],'%s');

fprintf('NETWORK: %s\n', ntw);
fprintf('BAND: %s\n', band);
fprintf('TASK: %s\n\r', task_name{task});

%% Consider a random subset of subjects,
% increase the size of the dataset at each iteration,
% save the performance and repeat for 20 trials to avoid the bias in the
% random selection.

trial = 20;
performanceTot = cell(trial,2);
for v = 1:trial    %prove con subsets diversi
    rng(v)
    s = rng;

    subjectsIdx = [1:length(totSubjects)];

    performancePscore = zeros(10,1);
    performanceIDrate = zeros(10,1);
    randomIdx = [];
    for z = 1:10     
        %take at random 5 subjects
        idx = randsample(subjectsIdx, 5);
        subjectsIdx = setdiff(subjectsIdx, idx);
        randomIdx = [randomIdx, idx];

        subjects = {};
        for i = 1:length(randomIdx)
            subjects{i,1} = totSubjects{randomIdx(i)};   %subset of subjects
        end
        disp(size(subjects,1));
        %% LOAD MATRICES
        %FC_matrices is an array #subjects x 5. The columns represent:
        %the name of the subjects (1), the rest matrices (2, 3, 4) and the task one (5).
        %The task matrix is computed as the mean matrix in the case of mH o mF.
        %NB: solve cases in which some matrices are not present

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

        %array with FC information (upper triangular part of each FC matrix)
        FC_arrays = cell(length(subjects),size(FC_matrices,2));
        for i = 1:length(subjects)
            FC_arrays{i,1} = subjects{i};

            for j=2:size(FC_matrices,2)
                R = triu(FC_matrices{i, j}); %upper triangular matrix 
                %convert into array mantaining the indexes (read the upper triangle of the matrix from left to
                %right and from top to bottom)
                Rt = R.';
                temp  = (1:size(Rt,1)).' >= (1:size(Rt,2));
                FC_arrays{i,j} = Rt(temp);        
            end
        end

        %% APPLY DEGREE-NORMALIZATION and save only the upper triangular part

        fprintf('DEGREE-NORMALIZATION\n\r');

        for i=1:length(subjects)  
            for j=2:size(FC_matrices,2)
                sub_m = degreeNorm(FC_matrices{i,j}); %degree-norm applied to matrix with only some nodes
                FC_matrices{i,j} = triu(sub_m);
            end
        end

        %array with FC information
        FC_arrays = cell(length(subjects),size(FC_matrices,2));
        for i = 1:length(subjects)
            subj = subjects{i};
            FC_arrays{i,1} = subj;

            for j=2:size(FC_matrices,2)
                R = FC_matrices{i, j}; %upper triangular matrix 
                %convert into array mantaining the indexes (read the upper triangle of the matrix from left to
                %right and from top to bottom)
                Rt = R.';
                temp  = (1:size(Rt,1)).' >= (1:size(Rt,2));
                FC_arrays{i,j} = Rt(temp);        
            end
        end


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
        n = length(nodes);
        ind = find(tril(ones(n,n)));
        for i=1:length(subjects)
            p = FC_pred{i,2};
            P = zeros(n);
            P(ind) = p;
            FC_pred{i,3} = P';
        end

        %save FC_true as a matrix
        FC_true = cell(length(subjects),3);
        for i=1:length(subjects)
            FC_true{i,1} = subjects{i};
            FC_true{i,2} = FC_arrays{i,5};
            FC_true{i,3} = FC_matrices{i,5};
        end

        %% Correlation (with Pearson) - Identifiability matrix

        fprintf('IDENTIFIABILITY MATRIX\n');
        IDmatrix = zeros(length(subjects));
        for i = 1:length(subjects)
            for j = 1:length(subjects)
                IDmatrix(i,j) = corr(FC_pred{i,2},FC_true{j,2}, 'rows', 'complete');
            end         
        end


        if ntw == "AllNtw"
            ntwName = 'AllNtw';
        else
            ntwName = m.conn.P2N{(find(strcmp( m.conn.P2N, ntw))),2};
            ntwName = strtrim(ntwName);
        end

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
            correlation = IDmatrix(i,:);
            [m, idx] = max(correlation);
            if idx == i
                %disp(i);
                count = count+1;
            end
        end
        IDrate = count/length(subjects)*100;
        fprintf('Identification Rate: %3.3f\n', IDrate);

        performancePscore(z) = mean_diag;
        performanceIDrate(z) = IDrate;
    end
    performanceTot{v,1} = performanceIDrate;
    performanceTot{v,2} = performancePscore;
end

%% Save
if ntwName == "CON"
    folderName = 'Results/CONT/'; 
else
    folderName = ['Results/', ntwName, '/'];
end

if ~exist(folderName, 'dir')
   mkdir(folderName)
end

fileName = ['PerfvsComplexity_', task_name{task}, '_', band, '_', ntwName];


%% Compute mean and std between the trials
IDrate_tot = [];
Pscore_tot = [];
for i=1:trial
    IDrate_tot = [IDrate_tot, performanceTot{i,1}];
    Pscore_tot = [Pscore_tot, performanceTot{i,2}];
end
mean_IDrate = mean(IDrate_tot,2);
std_IDrate = std(IDrate_tot,0,2);
mean_Pscore = mean(Pscore_tot,2);
std_Pscore = std(Pscore_tot,0,2);

%% Write in a file
fileID = fopen([folderName ,fileName,'.txt'],'w');
fprintf(fileID, 'Mean P_score: \n');
fmt = [repmat('%2.2f   ', 1, length(subjects)), '%2.2f \n'];
fprintf(fileID, fmt, mean_Pscore);
fprintf(fileID, '\nStd P_score: \n');
fmt = [repmat('%2.2f   ', 1, length(subjects)), '%2.2f \n'];
fprintf(fileID, fmt, std_Pscore);

fprintf(fileID, '\n\nMean ID_rate: \n');
fmt = [repmat('%2.2f   ', 1, length(subjects)), '%2.2f \n'];
fprintf(fileID, fmt, mean_IDrate);
fprintf(fileID, '\nStd IDrate: \n');
fmt = [repmat('%2.2f   ', 1, length(subjects)), '%2.2f \n'];
fprintf(fileID, fmt, std_IDrate);
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

function y = cosineSimilarity(a,b)
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