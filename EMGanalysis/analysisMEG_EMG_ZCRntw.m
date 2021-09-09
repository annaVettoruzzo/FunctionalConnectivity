%% EMG analysis
%predict inter-individual differences in muscle activation duing motor
%tasks relating graph metrics extracted from FC matrices at rest with ZCR feature extracted from EMG data.
%
% Author: Anna Vettoruzzo
% September 2021

clear all;
close all;
clc;

%% Select band, ntw and task

datadirMEG = ['D:\Tesi\Data\'];
datadirEMG = ['D:\Tesi\Codice\EMG_analysis\Data'];

band = 'betalow'; %betahigh, betalow
ntw = 'Motor Network'; %'AllNtw' o altre reti

task = 4; %1 = LH, 2 = LF, 4 = RH, 5 = RF
task_name = {'LH', 'LF', '', 'RH', 'RF'};
task_namelong = {'Left Hand', 'Left Foot', '', 'Right Hand', 'Right Foot'};

subjects=textread([datadirMEG, '\subjectsTMP.txt'],'%s');  %read only the subjects with both MEG and EMG data

fprintf('NETWORK: %s\n', ntw);
fprintf('BAND: %s\n', band);
fprintf('TASK: %s\n\r', task_name{task});

%% Load FC data during resting state

fprintf('LOAD FC MATRICES\n\r');
FC_matrices = cell(length(subjects),4);
for i=1:length(subjects)
    FC_matrices{i,1} = subjects{i};
    
    %REST
    fname = [datadirMEG, subjects{i}, '\Results\FunctionalConnectivity\Static\Rest\3-Restin\', band, '\', subjects{i}, '_MEG_3-Restin_icablpcorr_', band, 'conn.mat'];
    m = load(fname);
    if ntw == "AllNtw"
        nodes = [1:size(m.conn.complete,1)].';
    else
        nodes = m.conn.P2N{(find(strcmp( m.conn.P2N, ntw))),3}; %find nodes belonging to the network selected
    end
    FC_matrices{i,2} = degreeNorm(m.conn.complete(nodes,nodes));
   
    fname = [datadirMEG, subjects{i}, '\Results\FunctionalConnectivity\Static\Rest\4-Restin\', band, '\', subjects{i}, '_MEG_4-Restin_icablpcorr_', band, 'conn.mat'];
    m = load(fname);
    FC_matrices{i,3} = degreeNorm(m.conn.complete(nodes,nodes));
    
    fname = [datadirMEG, subjects{i}, '\Results\FunctionalConnectivity\Static\Rest\5-Restin\', band, '\', subjects{i}, '_MEG_5-Restin_icablpcorr_', band, 'conn.mat'];
    m = load(fname);
    FC_matrices{i,4} = degreeNorm(m.conn.complete(nodes,nodes));
end

%array with FC information
FC_arrays = cell(length(subjects),size(FC_matrices,2));
for i = 1:length(subjects)
    subj = subjects{i};
    FC_arrays{i,1} = subj;
    
    for j=2:size(FC_matrices,2)
        R = triu(FC_matrices{i, j}); %upper triangular matrix 
        %convert into array mantaining the indexes (read the upper triangle of the matrix from left to
        %right and from top to bottom)
        Rt = R.';
        temp  = (1:size(Rt,1)).' >= (1:size(Rt,2));
        FC_arrays{i,j} = Rt(temp);        
    end
end

%% Load EMG data
%EMG_data is a cell array with 3 columns to contain: the name of the
%subject, the average RMS in 10-Motort, the average RMS in 11-Motort
EMG_data = cell(length(subjects),3);
for i=1:length(subjects)
    EMG_data{i,1} = subjects{i};
    
    fname10 = [datadirEMG, '\', subjects{i}, '\preProcessing\10-Motort\EMGtaskdata.mat'];
    emg = struct2cell(load(fname10).emgData)';
    idx_task = find(strcmp(emg(:,1), task_namelong{task})); %find trial associated with the task
    idx_label = find(emg{1,3} == task_name{task});
    
    zc_v = zeros(1,length(idx_task));
    for j=1:length(idx_task)
        s = emg{idx_task(j),2};
        s = s(idx_label,:);
        zc_v(j) = mean(abs(diff(sign(s))));  %compute ZC for each segment
    end
    EMG_data{i,2} = mean(zc_v); %average all RMS
    
    fname11 = [datadirEMG, '\', subjects{i}, '\preProcessing\11-Motort\EMGtaskdata.mat'];
    emg = struct2cell(load(fname11).emgData)';
    idx_task = find(strcmp(emg(:,1), task_namelong{task})); %find trial associated with the task
    idx_label = find(emg{1,3} == task_name{task});
    for j=1:length(idx_task)
        s = emg{idx_task(j),2};
        s = s(idx_label,:);
        zc_v(j) = mean(abs(diff(sign(s))));  %compute ZC for each segment
    end
    EMG_data{i,3} = mean(zc_v); %average all RMS
end

%average the RMS from the two runs
EMG_ZC = cell(length(subjects),2);
for i=1:length(subjects)
    EMG_ZC{i,1} = subjects{i};
    EMG_ZC{i,2} = (EMG_data{i,2}+EMG_data{i,3})/2;
end

%% Extract features from FC matrices
FC_features = cell(length(subjects), 1);

for i = 1:length(subjects)
    FC_features{i,1} = subjects{i};
    FC1 = FC_matrices{i,2};
    FC2 = FC_matrices{i,3};
    FC3 = FC_matrices{i,4};
    
    %weighted degree
    wdegree1 = sum(FC1,2, 'omitnan')./sum(~isnan(FC1),2);
    wdegree2 = sum(FC2,2, 'omitnan')./sum(~isnan(FC2),2);
    wdegree3 = sum(FC3,2, 'omitnan')./sum(~isnan(FC3),2);
    FC_features{i,2} = mean([mean(wdegree1);mean(wdegree2);mean(wdegree3)]);
    disp(["WD mean = "+ FC_features{i,2}]);
    disp(["WD std = " + std([mean(wdegree1);mean(wdegree2);mean(wdegree3)])]);
    
    FC1(isnan(FC1))=0;
    G1 = graph(FC1, 'upper');
    FC2(isnan(FC2))=0;
    G2 = graph(FC2, 'upper');
    FC3(isnan(FC3))=0;
    G3 = graph(FC3, 'upper');
    %average path length
    shortestPath1 = avgPathLength(G1);
    shortestPath2 = avgPathLength(G2);
    shortestPath3 = avgPathLength(G3);
    
    FC_features{i,3} = mean([shortestPath1; shortestPath2; shortestPath3]);
    disp(["PL mean = "+ FC_features{i,3}]);
    disp(["PL std = " + std([shortestPath1; shortestPath2; shortestPath3])]);
    
    %efficiency
    E1 = efficiency(G1);
    E2 = efficiency(G2);
    E3 = efficiency(G3);
    
    FC_features{i,4} = mean([E1;E2;E3]);
    disp(["E mean = "+ FC_features{i,4}]);
    disp(["E std = " + std([E1; E2; E3])]);
    
end
%% Regression models 
%a model for each subject and each link

%compute beta coefficients for each subjects
fprintf('START LEARNING\n\r');
betas = cell(length(subjects), 1);
for i=1:length(subjects)
   betas{i,1} = subjects{i};
   
   %consider only FC values not NaN
   x = [FC_features{i,2}, FC_features{i,3}, FC_features{i,4}];
   %add the term for the intercept
   X = [ones(size(x,1),1), x]; 

   y = EMG_ZC{i,2}; %feature extracted from EMG data 

   beta = pinv(X)*y;
   betas{i,2} = beta;  
end


fprintf('LEAVE-ONE-OUT PREDICTION\n\r');
EMG_ZC_pred = cell(length(subjects),2);
for i = 1:length(subjects)
    LOOsubject = subjects{i};
    EMG_ZC_pred{i} = LOOsubject;
    
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
    x = [FC_features{i,2}, FC_features{i,3}, FC_features{i,4}];
    X = [ones(size(x,1),1), x]; 
    
    pred = X*avg_beta;
    EMG_ZC_pred{i, 2} = mean(pred);
end    

%% Correlation (with Pearson) - Identifiability matrix

fprintf('IDENTIFIABILITY MATRIX\n');
IDmatrix = zeros(length(subjects));
for i = 1:length(subjects)
    for j = 1:length(subjects)
        IDmatrix(i,j) = mse(EMG_ZC_pred{i,2},EMG_ZC{j,2});
    end         
end

if ntw == "AllNtw"
    ntwName = 'AllNtw';
else
    ntwName = m.conn.P2N{(find(strcmp( m.conn.P2N, ntw))),2};
    ntwName = strtrim(ntwName);
end

%normalization IDmatrix in [0 1]
maxElem = max(IDmatrix,[],'all');
minElem = min(IDmatrix,[],'all');
norm_IDmatrix = (IDmatrix - minElem)/(maxElem-minElem);
%tmp_IDmatrix = IDmatrix/max(IDmatrix,[],'all');

%show figure (elements with negative correlation are shown in white)
figName = [task_name{task}, '_', band, '_', ntwName];
figure('name', figName);
h = heatmap(norm_IDmatrix);
%h.MissingDataColor = [1,1,1];
h.caxis([0,1]); %to set the same scale
cmap = flipud(jet(20));
cmap = cmap(2:11,:);
h.Colormap = cmap;
h.FontSize = 12;
h.Title = ['Normalized IDmatrix with ZCR', ' - ',task_name{task},' ',band,' ',ntwName];
h.xlabel('Subjects (true)');
h.ylabel('Subjects (predicted)');
%save figure
if ntwName == "CON"
    folderName = 'D:\Tesi\Codice\EMG_analysis\ResultsEMG\ZCR\NtwApproach\CONT/'; 
else
    folderName = ['D:\Tesi\Codice\EMG_analysis\ResultsEMG\ZCR\NtwApproach\', ntwName, '/'];
end

if ~exist(folderName, 'dir')
   mkdir(folderName)
end
saveas(h, [folderName , figName, '.png']);

%% Performance of the model: 

fprintf('PERFORMANCE OF THE MODEL\n');
%PREDICTION SCORE (ref. Tavor et al.)
mean_diag = mean(diag(IDmatrix));
std_diag = std(diag(IDmatrix));
disp('Prediction score : ');
disp( mean_diag);
disp('Std prediction score : ');
disp( std_diag);


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

%% Write in a file
totWD = [];
totPL = [];
totE = [];
for i = 1:length(subjects)
    totWD = [totWD; FC_features{i,2}];
    totPL = [totPL; FC_features{i,3}];
    totE = [totE; FC_features{i,4}];
end
meanWD = mean(totWD);
stdWD = std(totWD);
meanPL = mean(totPL);
stdPL = std(totPL);
meanE = mean(totE);
stdE = std(totE);

%% Write in a file
fileID = fopen([folderName ,band, '_', ntwName,'.txt'],'w');
fprintf(fileID, 'Weighted Degree : %3.3e ', meanWD);
fprintf(fileID, '± %3.3e \n\n', stdWD);
fprintf(fileID, 'Path Length : %3.3e ', meanPL);
fprintf(fileID, '± %3.3e \n\n', stdPL);
fprintf(fileID, 'Efficiency : %3.3e ', meanE);
fprintf(fileID, '± %3.3e \n\n', stdE);
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

%% Average path length
function pathLength = avgPathLength(G)
    n = size(G.Nodes,1);
    x = [];
    for c = 1:n
        for dest=1:n
            if(c~=dest)
                [~,d] = shortestpath(G,c,dest);
                x = [x, d];
            end
        end
    end
    pathLength = mean(x);
end

%% Efficiency
function E = efficiency(G)
    n = size(G.Nodes,1);
    d = distances(G);
    d = 1./d;
    d(d == Inf) =0;
    E = 1/(n*(n-1))*sum(d, 'all');
end