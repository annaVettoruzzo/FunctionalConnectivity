%% EMG analysis with a broadband or a band-limited approach
%predict inter-individual differences in muscle activation duing motor
%tasks relating FC values with RMS feature extracted from EMG data.
%
% Author: Anna Vettoruzzo
% September 2021


clear all;
close all;
clc;

%% Select ntw, task and band

datadirMEG = ['D:\Tesi\Data\'];
datadirEMG = ['D:\Tesi\Codice\EMG_analysis\Data'];

band = 'betalow'; %betahigh, betalow
ntw = 'Motor Network'; %'AllNtw' o altre reti

task = 5; %1 = LH, 2 = LF, 4 = RH, 5 = RF
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

%% Load EMG data (RMS) and sekect a broadband or a band-limited approach
%EMG_data is a cell array with 3 columns to contain: the name of the
%subject, the average RMS in 10-Motort, the average RMS in 11-Motort

%perform the analysis only selecting the frequency band of interest also
%for EMG signals (band-limited)
%NB: the processing doesn't allow to select alpha band
selectBand = "False";    
fs = 2034.5;

EMG_data = cell(length(subjects),3);
for i=1:length(subjects)
    EMG_data{i,1} = subjects{i};
    
    fname10 = [datadirEMG, '\', subjects{i}, '\preProcessing\10-Motort\EMGtaskdata.mat']; %load data after cleaning
    emg = struct2cell(load(fname10).emgData)';
    idx_task = find(strcmp(emg(:,1), task_namelong{task})); %find trial associated with the task
    idx_label = find(emg{1,3} == task_name{task});
    
    rms_v = zeros(1,length(idx_task));
    for j=1:length(idx_task)
        s = emg{idx_task(j),2};
        s = s(idx_label,:);
        if selectBand == "True"
            s = filtering(s, band, fs);
        end
        rms_v(j) = rms(s);  %compute RMS for each segment
    end
    EMG_data{i,2} = mean(rms_v); %average all RMS
    
    fname11 = [datadirEMG, '\', subjects{i}, '\preProcessing\11-Motort\EMGtaskdata.mat'];
    emg = struct2cell(load(fname11).emgData)';
    idx_task = find(strcmp(emg(:,1), task_namelong{task})); %find trial associated with the task
    idx_label = find(emg{1,3} == task_name{task});
    for j=1:length(idx_task)
        s = emg{idx_task(j),2};
        s = s(idx_label,:);
        if selectBand == "True"
            s = filtering(s, band, fs);
        end
        rms_v(j) = rms(s);  %compute RMS for each segment
    end
    EMG_data{i,3} = mean(rms_v); %average all RMS
end

%average the RMS from the two runs
EMG_rms = cell(length(subjects),2);
for i=1:length(subjects)
    EMG_rms{i,1} = subjects{i};
    EMG_rms{i,2} = (EMG_data{i,2}+EMG_data{i,3})/2;
end

%% Regression models 
%a model for each subject 

%compute beta coefficients for each subjects
fprintf('START LEARNING\n\r');
betas = cell(length(subjects), 1);
for i=1:length(subjects)
   betas{i,1} = subjects{i};
   
   %consider only FC values not NaN
   x = mean(cat(2, FC_arrays{i,2}, FC_arrays{i,3}, FC_arrays{i,4}),2)'; %REST 
   %add the term for the intercept
   X = [ones(size(x,1),1), x];  
   ind = find(all(~isnan(X),1)); 
   ind_nan = find(all(isnan(X),1)); 
   X = X(:,ind);
   
   y = EMG_rms{i,2}; %feature extracted from EMG data 
   
   invX = pinv(X);
   beta(ind) = invX*y;
   beta(ind_nan) = NaN;
   betas{i,2} = beta';  
end


fprintf('LEAVE-ONE-OUT PREDICTION\n\r');
EMG_rms_pred = cell(length(subjects),2);
for i = 1:length(subjects)
    LOOsubject = subjects{i};
    EMG_rms_pred{i} = LOOsubject;
    
    %compute avg beta 
    loo_betas = zeros(size(betas{i,2},1),length(subjects)-1);
    cnt = 1;
    for j=setdiff(1:length(subjects),i)
        loo_betas(:,cnt) = betas{j,2};
        cnt = cnt+1;
    end
    avg_beta = mean(loo_betas,2);   %mean only in the position where all elements are values
    
    %prediction
    pred= [];
    x = mean(cat(2, FC_arrays{i,2}, FC_arrays{i,3}, FC_arrays{i,4}),2)';
    X = [ones(size(x,1),1), x];  
     
    ind = find(all(~isnan(X),1)); 
    ind2 = find(all(~isnan(avg_beta),2)); 
    [val,pos]=intersect(ind, ind2);
    pred = X(:,val)*avg_beta(val,:);
    EMG_rms_pred{i, 2} = pred;
end    

%% Correlation (with MSE) - Identifiability matrix

fprintf('IDENTIFIABILITY MATRIX\n');
IDmatrix = zeros(length(subjects));
for i = 1:length(subjects)
    for j = 1:length(subjects)
        IDmatrix(i,j) = mse(EMG_rms_pred{i,2},EMG_rms{j,2});
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

%show figure 
figName = [task_name{task}, '_', band, '_', ntwName];
figure('name', figName);
h = heatmap(norm_IDmatrix);
h.caxis([0,1]); %to set the same scale
h.FontSize = 12;
cmap = (jet(20));
cmap = flipud(cmap(10:19,:));
h.Colormap = cmap;
h.Title = ['Normalized IDmatrix with RMS', ' - ',task_name{task},' ',band,' ',ntwName];
h.xlabel('Subjects (true)');
h.ylabel('Subjects (predicted)');

%save figure
if selectBand == "True"
    if ntwName == "CON"
        folderName = 'D:\Tesi\Codice\EMG_analysis\ResultsEMG\RMS\FirstApproach\bandLimited/CONT/'; 
    else
        folderName = ['D:\Tesi\Codice\EMG_analysis\ResultsEMG\RMS\FirstApproach\bandLimited/', ntwName, '/'];
    end
else
    if ntwName == "CON"
        folderName = 'D:\Tesi\Codice\EMG_analysis\ResultsEMG\RMS\FirstApproach/CONT/'; 
    else
        folderName = ['D:\Tesi\Codice\EMG_analysis\ResultsEMG\RMS\FirstApproach/', ntwName, '/'];
    end
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
    error = IDmatrix(i,:);
    [~, idx] = min(error);
    if idx == i
        disp(i);
        count = count+1;
    end
end
IDrate = count/length(subjects)*100;
fprintf('Identification Rate: %3.3f\n', IDrate);


%% Write in a file
fileID = fopen([folderName ,figName,'.txt'],'w');
fprintf(fileID, 'Min value of similarity: %3.3e \n\n', min(IDmatrix,[],'all'));
%Pscore
fprintf(fileID, 'Prediction score : %3.3e ', mean_diag);
fprintf(fileID, '± %3.3e \n\n', std_diag);
%IDrate
fprintf(fileID, 'Identification rate: %3.3f%%\n\n', IDrate);
%Min value IDmatrix
fprintf(fileID, 'Min value in IDmatrix: %3.3e\n', min(IDmatrix,[],'all'));
%Max value IDmatrix
fprintf(fileID, 'Max value in IDmatrix: %3.3e\n', max(IDmatrix,[],'all'));
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

%% Apply filtering operation to select the correct frequency band
function filtered_signal = filtering(s, f_band, fs)
    if string(f_band) == "betalow"
        Lowpass = 26;
        Highpass = 15;
    elseif string(f_band) == "betahigh"
        Lowpass = 35;
        Highpass = 26;
    else 
        disp('ERROR: only betalow and betahigh bands are ok!');
        return;
    end
    
    % Lowpass filter 
    NyqFreq = fs/2;
    Wn = Lowpass/NyqFreq;
    [B,A] = butter (4,Wn,'low');
    % run filter
    low_filtered = filtfilt(B,A,s);
    % Highpass filter 
    Wo = Highpass/NyqFreq;
    [D,C] = butter (4,Wo,'high');
    % run filter
    filtered_signal = filtfilt(D,C,low_filtered);
end
