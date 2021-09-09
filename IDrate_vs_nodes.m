%% Functional connectivity analysis
%  model inspired by Tavor et al. to predict inter-individual differences 
% in motor tasks exploiting MEG-based resting-state FC matrices.
% Identify if there is a linear relation between IDrate and #nodes of the
% networks.
%
% Author: Anna Vettoruzzo
% Septemeber 2021

clear all;
close all;
clc;

%% Select band and task
datadir = ['D:\Tesi\Data\'];

band = 'betalow'; %betahigh, betalow, alpha
task = 6; %1 = LH, 2 = LF, 4 = RH, 5 = RF, 6 = meanH, 7 = meanF 
task_name = {'LH', 'LF', '', 'RH', 'RF', 'meanH', 'meanF'};

subjects=textread([datadir, '\subjects.txt'],'%s');

fprintf('BAND: %s\n', band);
fprintf('TASK: %s\n\r', task_name{task});

%% Nodes and IDrate for each network

%NODES
fname = [datadir, subjects{1}, '\Results\FunctionalConnectivity\Static\Rest\3-Restin\', band, '\', subjects{1}, '_MEG_3-Restin_icablpcorr_', band, 'conn.mat'];
m = load(fname);
ntwData = cell(11,3); %10 ntws plus AllNtw
ntwData{1,1} = 'AllNtw';
ntwData{1,2} = size(m.conn.complete,1);
for i = 1:10
    ntwData{i+1,1} = strtrim(m.conn.P2N{i,2});
    ntwData{i+1,2} = length(m.conn.P2N{i,3});
end

%IDRATE
for i = 1:size(ntwData,1)
   ntwName = ntwData{i,1}; 
   if ntwName == "CON"
       dirID = ['Results\CONT\', task_name{task}, '_', band, '_', ntwName, '.txt']; 
   else
       dirID = ['Results\', ntwName, '\', task_name{task}, '_', band, '_', ntwName, '.txt']; 
   end
   fid = fopen(dirID,'r');
   while ~feof(fid)
       st = fgetl(fid);
       if  contains(st,'Identification rate: ')
           stread = textscan(st,'%s %s %f');
           ntwData{i,3} = stread{3};
       end
   end
   fclose(fid);
end


%% Linear relation between # nodes and ID rate in the networks?
%To do so compute the Pearson correlation coefficient and the pvalue:
%if pvalue < 0.05 --> data are correlated
%otherwise not.

nodes = cell2mat(ntwData(:,2));
IDrate = cell2mat(ntwData(:,3));
[r,p] = corr(IDrate,nodes);

figName = ['IDrate vs #Nodes - ', task_name{task}, '_', band];
figure('name', figName);
c = 1:11;
s=scatter(nodes, IDrate, 25, c, 'filled');
colormap(gca,'jet')
xlim([0 40]);
xlabel('#Nodes');
ylabel('IDrate');
title(['IDrate vs #Nodes - ', task_name{task}, ' ', band]);
grid on;
labels = strings(size(ntwData,1),1);
for i=1:size(ntwData,1)
    labels(i) = ntwData{i,1};
end 
labelpoints(nodes,IDrate,labels, 'FontSize', 10);
% hold on
% x = 0:max(nodes);
% plot(x, r*x);
% hold off;

%% Save image and write in a file
folderName = ['Results\IDRATEvsNODES\', task_name{task},'\'];
if ~exist(folderName, 'dir')
   mkdir(folderName)
end
figName = ['IDrateNodes_',task_name{task}, '_', band];
saveas(s, [folderName , figName, '.png']);

fileID = fopen([folderName,'_', band, '.txt'],'w');
fprintf(fileID, 'IDrate vs #Nodes\n');
fprintf(fileID, 'Correlation coefficient: %3.3f\n', r);
fprintf(fileID, 'P-value: %3.3f ', p);
fclose(fileID);