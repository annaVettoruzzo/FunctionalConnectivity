%% Remove bad epochs and bad channels from EMG data
%
% Author: Anna Vettoruzzo
% September 2021


clc;
close all;
clear all;

%% Set paths
datadirEMG = ['D:\Tesi\Codice\EMG_analysis\Data\'];
datadirMEG = ['D:\Tesi\Data\'];

subjects=textread([datadirMEG, '\subjectsTMP.txt'],'%s');
run = '11-Motort';

%% Load and clean data
disp('LOAD DATA');
for k = 1:length(subjects)
    subj = subjects{k};

    %Read EMG data
    filename = [datadirEMG, subj,'\preProcessing\',run,'\EMGfiltered.mat'];
    EMGdata = load(filename).EMGfiltered;   %LH, LF, RH, RF 
    tabella = load(filename).tabella;

    fs = 2.0345e+03;
    MS_TRIALLENGTH = 0.8;
    trl=round(MS_TRIALLENGTH*fs); %number of samples in each trial

    %Clean EMG signal from bad segments
    %read the file and create an array with the bad segments
    disp(['CLEAN DATA subject ', subj]);
    
    badsegments_fileName = [datadirMEG, subj,'\',subj,'_MEG_',run,'_baddata_badsegments.txt'];
    fileID = fopen(badsegments_fileName, 'r');
    tline = fgetl(fileID);
    while contains(tline, 'badsegment.all') ~= 1
        tline = fgetl(fileID);
    end

    bad_segments = [];
    while tline ~= -1 %~feof(fileID)
        num = regexp(tline,'\d*','Match');
        interval = [str2double(num{1}), str2double(num{2})];
        bad_segments = [bad_segments; interval];
        tline = fgetl(fileID);
    end
    fclose(fileID);

    %remove bad segments from EMG signal
    EMGdatacleaned = EMGdata;
    for i = length(bad_segments):-1:1
        start = bad_segments(i,1);
        stop = bad_segments(i,2);
        EMGdatacleaned(:, start:stop) = [nan];
    end

    %Note: some trials can be empty or without some samples because they
    %corresponds to bad samples.
    emgData=[];
    classConn={'Left Hand','Left Foot','Right Hand','Right Foot'};
    cl=[1,2,4,5];
    for itemg1=1:size(tabella,1)
       junk2=struct;
       junk2.type=classConn{find(cl==tabella(itemg1,2))};
       junk2.EMG=EMGdata(:,tabella(itemg1,4):tabella(itemg1,4)+trl-1);
       junk2.EMG(:,find(isnan(junk2.EMG(1,:)))) = [];
       junk2.label=["LH"; "LF"; "RH"; "RF"];
       junk2.first=tabella(itemg1,6);
       emgData=[emgData;junk2];
    end

    %% Save EMG data

    outpath=[datadirEMG, subj, '/preProcessing/', run];
    if ~exist(outpath, 'dir')
        mkdir(outpath);
    end
    save([outpath,'/','EMGtaskdata.mat'],'emgData');
end