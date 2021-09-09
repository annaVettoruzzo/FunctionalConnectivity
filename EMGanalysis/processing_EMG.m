%% Processing EMG signals to remove noise and artifacts.
%
% Author: Anna Vettoruzzo
% September 2021


clc;
close all;
clear all;

%% Set paths
addpath D:\Tesi\Codice\fieldtrip-20210311
ft_defaults %sets the defaults and configures up the minimal required path settings. 

datadirMEG = ['D:\Tesi\Data\'];
datadirEMG = ['D:\Tesi\Codice\EMG_analysis\Data'];
run = '10-Motort';

subjects=textread([datadirMEG, '\subjectsTMP.txt'],'%s');

%% EMG data (time domain)

for k=1:length(subjects)
    subj = subjects{k};

    %Read EMG data
    classConn={'Left Hand','Left Foot','Right Hand','Right Foot'};
    cl=[1,2,4,5];

    %c,rfDC : c = continuously recorded data, rfDC = no filtering at all
    %MEG data is composed of many channels and many time points. It contains a
    %sample =single number for every ChannelxTime point.
    filename = [datadirMEG, subj,'/Experiments/',subj,'_MEG/Scans/',run,'/Resources/4D/c,rfDC'];

    event = ft_read_event(filename); %returns a structure with trigger info

    %This function also exploit theinfo contained in the file 'config' in the
    %same directory as c,rfDC
    hdr   = ft_read_header(filename); %returns a structure with the header info

    %NchansxNsample matrix
    %data = ft_read_data(filename); %returns a 2D array with the data

    cfg          = [];
    cfg.datafile = filename;
    cfg.channel = 'EMG';
    tmpdata     = ft_preprocessing(cfg);

    %TRIGGER of interest
    %LH (1) = 256 + 22; RH (4) = 256 + 70; LF (2) = 256 + 38; RF (5) = 256 + 134
    buff=[];
    last=0;
    for itemg1=1:max(size(event))
       if(event(itemg1).value==278) %LEFT HAND
           cc=double((last~=event(itemg1).value)); %1 se l'evento itemg1 è diverso da last
           last=event(itemg1).value;
           buff=[buff;0,1,0,event(itemg1).sample,0,cc];
       elseif(event(itemg1).value==390) %RIGHT FOOT
           cc=double((last~=event(itemg1).value));
           last=event(itemg1).value;
           buff=[buff;0,5,0,event(itemg1).sample,0,cc];
       elseif(event(itemg1).value==326) %RIGHT HAND
           cc=double((last~=event(itemg1).value));
           last=event(itemg1).value;
           buff=[buff;0,4,0,event(itemg1).sample,0,cc];
       elseif(event(itemg1).value==294) %LEFT FOOT
           cc=double((last~=event(itemg1).value));
           last=event(itemg1).value;
           buff=[buff;0,2,0,event(itemg1).sample,0,cc];        
       end
    end
    tabella=buff;

    %% Plot time domain
    % Sampling frequency
    fs = hdr.Fs; % [Hz]
    % Number of data samples
    N = hdr.nSamples; 
    % Signal duration
    T = N/fs; % [s]

    % Time vector defined once fs is known
    t = 0:1/fs:T-1/fs;  %t=tmpdata.time{1}
    % Samples vector
    n = 1:N;

    %understand tmpdata: it contains 4 signals (LF, LH, RF, RH) and each one
    %has 1420391 samples. In every signal 8 peaks reflect the 8 tasks in each run.
     figure(1);
     subplot(4,1,1)
      lf = tmpdata.trial{1}(1,:);
     plot(tmpdata.time{1}, lf); %LF
     title('Left Foot');
     grid on;
     ylim([-2e-3, 2e-3]);
     xlim([170,200]);
     ylabel('[V]');
     subplot(4,1,2)
      rf = tmpdata.trial{1}(3,:);
     plot(tmpdata.time{1}, rf); %RF
     title('Right Foot');
     grid on;
     ylim([-2e-3, 2e-3]);
     xlim([170,200]);
     ylabel('[V]');
     subplot(4,1,3)
      lh = tmpdata.trial{1}(2,:);
     plot(tmpdata.time{1}, lh); %LH
     title('Left Hand');
     grid on;
     ylim([-2e-3, 2e-3]);
     xlim([170,200]);
     ylabel('[V]');
     subplot(4,1,4)
      rh = tmpdata.trial{1}(4,:);
     plot(tmpdata.time{1}, rh); %RH
     title('Right Hand');
     grid on;
     ylim([-2e-3, 2e-3]);
     xlim([170,200]);
     xlabel('Time [s]');
     ylabel('[V]');
     sgtitle('Raw signals (time domain)') ;
     


    %% EMG data (frequency domain)
    LF = 1/fs*fft(tmpdata.trial{1}(1,:));   % apply the DFT operator
    N_S = length(LF);    % number of frequency samples
    F = fs/N_S;         % frequency resolution
    f = F*(0:N_S-1);    % frequency vector
    
    RF = 1/fs*fft(tmpdata.trial{1}(3,:));   % apply the DFT operator
    LH = 1/fs*fft(tmpdata.trial{1}(2,:));   % apply the DFT operator
    RH = 1/fs*fft(tmpdata.trial{1}(4,:));   % apply the DFT operator

    % Plot the EMG signal in the frequency domain
    figure(3)
    subplot(2,2,1)
    plot(f,abs(LF).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 100 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Left Foot')
    subplot(2,2,2)
    plot(f,abs(RF).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 500 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Right Foot')
    subplot(2,2,3)
    plot(f,abs(LH).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 500 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Left Hand')
    subplot(2,2,4)
    plot(f,abs(RH).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 500 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Right Hand')
    sgtitle('Raw signals (frequency domain)');
    
    %% Apply high-pass and low-pass filters
    % Lowpass filter 450Hz
    NyqFreq = fs/2;
    Lowpass = 450;
    Wn = Lowpass/NyqFreq;
    [B,A] = butter (4,Wn,'low');
    % run filter
    lf_low_filtered = filtfilt(B,A,lf);
    rf_low_filtered = filtfilt(B,A,rf);
    lh_low_filtered = filtfilt(B,A,lh);
    rh_low_filtered = filtfilt(B,A,rh);
    % Highpass filter 15Hz
    Highpass = 15;
    Wo = Highpass/NyqFreq;
    [D,C] = butter (4,Wo,'high');
    % run filter
    high_pass_lf_filtdata = filtfilt(D,C,lf_low_filtered);
    high_pass_rf_filtdata = filtfilt(D,C,rf_low_filtered);
    high_pass_lh_filtdata = filtfilt(D,C,lh_low_filtered);
    high_pass_rh_filtdata = filtfilt(D,C,rh_low_filtered);

    % Plot filtered signals
    figure(4);
    subplot(4,1,1)
    plot(tmpdata.time{1}, high_pass_lf_filtdata); %LF
    title('Left Foot');
    grid on;
    ylim([-2e-3, 2e-3]);
    ylabel('[V]');
    subplot(4,1,2)
    plot(tmpdata.time{1}, high_pass_rf_filtdata); %RF
    title('Right Foot');
    grid on;
    ylim([-2e-3, 2e-3]);
    ylabel('[V]');
    subplot(4,1,3)
    plot(tmpdata.time{1}, high_pass_lh_filtdata); %LH
    title('Left Hand');
    grid on;
    ylim([-2e-3, 2e-3]);
    ylabel('[V]');
    subplot(4,1,4)
    plot(tmpdata.time{1}, high_pass_rh_filtdata); %RH
    title('Right Hand');
    grid on;
    ylim([-2e-3, 2e-3]);
    xlabel('Time [s]');
    ylabel('[V]');
    sgtitle('Band-pass filtered signals (time domain)')
%     
    %FREQUENCY DOMAIN
    %LF
    LF_filtered = 1/fs*fft(high_pass_lf_filtdata);   % apply the DFT operator
    N_S = length(LF_filtered);    % number of frequency samples
    F = fs/N_S;         % frequency resolution
    f = F*(0:N_S-1);    % frequency vector
    %RF
    RF_filtered = 1/fs*fft(high_pass_rf_filtdata);   % apply the DFT operator
    %LH
    LH_filtered = 1/fs*fft(high_pass_lh_filtdata);   % apply the DFT operator
    %RH
    RH_filtered = 1/fs*fft(high_pass_rh_filtdata);   % apply the DFT operator

    % Plot the EMG signal in the frequency domain
    figure(5)
    subplot(2,2,1)
    plot(f,abs(LF_filtered).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 500 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Left Foot')
    subplot(2,2,2)
    plot(f,abs(RF_filtered).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 500 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Right Foot')
    subplot(2,2,3)
    plot(f,abs(LH_filtered).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 500 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Left Hand')
    subplot(2,2,4)
    plot(f,abs(RH_filtered).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 500 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Right Hand')
    sgtitle('Band-pass filtered signals (frequency domain)')

    %% Notch filtering for power line noise removal (60 Hz) and its third harmonic
    Fnotch = 60; % Notch Frequency
    q = 35; %quality factor
    w0 = Fnotch/(fs/2);
    bw = w0/35;
    Apass = 1; % Bandwidth Attenuation

    [b, a] = iirnotch (w0, bw, Apass); 
    Hd = dfilt.df2 (b, a);
    %fvtool(b,a);

    lf_filtered = filter(Hd,high_pass_lf_filtdata);
    rf_filtered = filter(Hd,high_pass_rf_filtdata);
    lh_filtered = filter(Hd,high_pass_lh_filtdata);
    rh_filtered = filter(Hd,high_pass_rh_filtdata);
    
    Fnotch = 180; % Notch Frequency
    q = 35; %quality factor
    w0 = Fnotch/(fs/2);
    bw = w0/35;
    Apass = 1; % Bandwidth Attenuation

    [b, a] = iirnotch (w0, bw, Apass); 
    Hd = dfilt.df2 (b, a);
    %fvtool(b,a);

    lf_filtered = filter(Hd,lf_filtered);
    rf_filtered = filter(Hd,rf_filtered);
    lh_filtered = filter(Hd,lh_filtered);
    rh_filtered = filter(Hd,rh_filtered);

    % Plot filtered signals
    figure(6);
    subplot(4,1,1)
    plot(tmpdata.time{1}, lf_filtered); %LF
    title('Left Foot');
    ylabel('[V]');
    grid on;
    ylim([-2e-3, 2e-3]);
    subplot(4,1,2)
    plot(tmpdata.time{1}, rf_filtered); %RF
    title('Right Foot');
    ylabel('[V]');
    grid on;
    ylim([-2e-3, 2e-3]);
    subplot(4,1,3)
    plot(tmpdata.time{1}, lh_filtered); %LH
    title('Left Hand');
    ylabel('[V]');
    grid on;
    ylim([-2e-3, 2e-3]);
    subplot(4,1,4)
    plot(tmpdata.time{1}, rh_filtered); %RH
    title('Right Hand');
    ylabel('[V]');
    grid on;
    ylim([-2e-3, 2e-3]);
    xlabel('Time [s]');
    ylabel('[V]');
    sgtitle('Filtered signals (time domain)')

    %FREQUENCY DOMAIN
    %LF
    LF_filtered = 1/fs*fft(lf_filtered);   % apply the DFT operator
    N_S = length(LF_filtered);    % number of frequency samples
    F = fs/N_S;         % frequency resolution
    f = F*(0:N_S-1);    % frequency vector
    %RF
    RF_filtered = 1/fs*fft(rf_filtered);   % apply the DFT operator
    %LH
    LH_filtered = 1/fs*fft(lh_filtered);   % apply the DFT operator
    %RH
    RH_filtered = 1/fs*fft(rh_filtered);   % apply the DFT operator

    % Plot the EMG signal in the frequency domain
    figure(8)
    subplot(2,2,1)
    plot(f,abs(LF_filtered).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 500 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Left Foot')
    subplot(2,2,2)
    plot(f,abs(RF_filtered).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 500 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Right Foot')
    subplot(2,2,3)
    plot(f,abs(LH_filtered).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 500 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Left Hand')
    subplot(2,2,4)
    plot(f,abs(RH_filtered).^2, 'k', 'LineWidth', 1) % plot the power spectrum
    grid on;                               
    axis([0 500 0 1e-7])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('Right Hand')
    sgtitle('Filtered signals (frequency domain)');

    %% Plot original signals and filtered ones in the same plot
    
    fig = figure(9);
    ax1 = subplot(4,1,1);
    plot(ax1, tmpdata.time{1}, lf, 'r'); %LF
    ylim([-2e-3, 2e-3]);
    title('Left Foot');
    ylabel('[V]');
    grid on;
    ax2 = subplot(4,1,2);
    plot(ax2, tmpdata.time{1}, rf, 'r'); %RF
    title('Right Foot');
    ylim([-2e-3, 2e-3]);
    ylabel('[V]');
    grid on;
    ax3 = subplot(4,1,3);
    plot(ax3, tmpdata.time{1}, lh, 'r'); %LH
    title('Left Hand');
    ylim([-2e-3, 2e-3]);
    ylabel('[V]');
    grid on;
    ax4 = subplot(4,1,4);
    plot(ax4, tmpdata.time{1}, rh, 'r'); %RH
    title('Right Hand');
    ylim([-2e-3, 2e-3]);
    ylabel('[V]');
    grid on;
    xlabel('Time [s]');
    sgtitle('Signals (time domain)')

    hold(ax1,'on');
    plot(ax1, tmpdata.time{1}, lf_filtered, 'b', 'linewidth', 1.5); %LF
    hold(ax1,'off');
    hold(ax2,'on');
    plot(ax2, tmpdata.time{1}, rf_filtered, 'b', 'linewidth', 1.5); %RF
    hold(ax2,'off');
    hold(ax3,'on');
    plot(ax3, tmpdata.time{1}, lh_filtered, 'b', 'linewidth', 1.5); %LH
    hold(ax3,'off');
    hold(ax4,'on');
    plot(ax4, tmpdata.time{1}, rh_filtered, 'b', 'linewidth', 1.5); %RH
    hold(ax4,'off');
    linkaxes([ax1 ax2 ax3 ax4], 'xy');
    hL = legend('raw', 'filtered');
    newPosition = [0.81 0.87 0.05 0.05];
    newUnits = 'normalized';
    set(hL,'Position', newPosition, 'Units', newUnits);
    
    %save figure
    outpath="D:/Tesi/Codice/EMG_analysis/Data/" +subj+ "/preProcessing/"+run;
    if ~exist(outpath, 'dir')
        mkdir(outpath);
    end
    name = [outpath + 'filtered.png'];
    saveas(fig, name)
    
%     %%FIGURA DI PROVA PER VEDERE PSD DEI SEGNALI ZOOMATI
%     start = interp1(t,1:length(t),180,'nearest');
%     e = interp1(t,1:length(t),200,'nearest');
%     ZLF = 1/fs*fft(lf_filtered(start:e));
%     ZRF = 1/fs*fft(rf_filtered(start:e));
%     ZLH = 1/fs*fft(lh_filtered(start:e));
%     ZRH = 1/fs*fft(rh_filtered(start:e));
%     N_S = length(ZLF);    % number of frequency samples
%     F = fs/N_S;         % frequency resolution
%     f = F*(0:N_S-1);    % frequency vector
%     % Plot the EEG signal in the frequency domain
%     figure(11)
%     subplot(2,2,1)
%     plot(f,abs(ZLF).^2, 'k', 'LineWidth', 1) % plot the power spectrum
%     grid on;                               
%     axis([0 500 0 max(abs(ZLF).^2)])
%     xlabel('Frequency [Hz]')
%     ylabel('P(f)')
%     title('LF')
%     subplot(2,2,2)
%     plot(f,abs(ZRF).^2, 'k', 'LineWidth', 1) % plot the power spectrum
%     grid on;                               
%     axis([0 500 0 max(abs(ZRF).^2)])
%     xlabel('Frequency [Hz]')
%     ylabel('P(f)')
%     title('RF')
%     subplot(2,2,3)
%     plot(f,abs(ZLH).^2, 'k', 'LineWidth', 1) % plot the power spectrum
%     grid on;                               
%     axis([0 500 0 max(abs(ZLH).^2)])
%     xlabel('Frequency [Hz]')
%     ylabel('P(f)')
%     title('LH')
%     subplot(2,2,4)
%     plot(f,abs(ZRH).^2, 'k', 'LineWidth', 1) % plot the power spectrum
%     grid on;                               
%     axis([0 500 0 max(abs(ZRH).^2)])
%     xlabel('Frequency [Hz]')
%     ylabel('P(f)')
%     title('RH')
%     sgtitle('Zoom Filtered signals (frequency domain)')
    
    %% Spettro di potenza medio 
    %Nsample = fs*t
    % every movement block is composed by 10 moviments; every moviments is
    % composed by 150 ms with a visual cue and 1050 ms to do the movement.
    %150 ms = 306 samples
    %1050ms = 2.1362e+03 samples

    %find all contraction segments (LF)
    lf_idx = find(tabella(:,2)==2);
    lf_start = tabella(lf_idx,4);
    lf_end = lf_start + round(0.8*fs);
    lf_fft = cell(length(lf_start),1);
    lf_fft_filt = cell(length(lf_start),1);
    %compute fft for each contraction and control if the frequency resolution
    %is always the same
    for i=1:length(lf_start)
        lf_signal = lf(:,lf_start(i):lf_end(i));
        lf_fft{i,1} = 1/fs*fft(lf_signal);   % apply the DFT operator  

        lf_signal_filt = lf_filtered(:,lf_start(i):lf_end(i));
        lf_fft_filt{i,1} = 1/fs*fft(lf_signal_filt);   % apply the DFT operator
        %lf_fft_filt{i,2} = fs/length(lf_fft_filt{i,1}); %frequency resolution    
    end
    %compute the mean of fft
    lf_fft_matrix = cell2mat(lf_fft);
    lf_mean_power = mean(lf_fft_matrix, 1);
    lf_fft_filt_matrix = cell2mat(lf_fft_filt);
    lf_mean_power_filt = mean(lf_fft_filt_matrix, 1);

   %find all contraction segments (RF)
    rf_idx = find(tabella(:,2)==5);
    rf_start = tabella(rf_idx,4);
    rf_end = rf_start + round(0.8*fs);
    rf_fft = cell(length(rf_start),1);
    rf_fft_filt = cell(length(rf_start),1);
    %compute fft for each contraction and control if the frequency resolution
    %is always the same
    for i=1:length(rf_start)
        rf_signal = rf(:,rf_start(i):rf_end(i));
        rf_fft{i,1} = 1/fs*fft(rf_signal);   % apply the DFT operator  

        rf_signal_filt = rf_filtered(:,rf_start(i):rf_end(i));
        rf_fft_filt{i,1} = 1/fs*fft(rf_signal_filt);   % apply the DFT operator
        %lf_fft_filt{i,2} = fs/length(lf_fft_filt{i,1}); %frequency resolution    
    end
    %compute the mean of fft
    rf_fft_matrix = cell2mat(rf_fft);
    rf_mean_power = mean(rf_fft_matrix, 1);
    rf_fft_filt_matrix = cell2mat(rf_fft_filt);
    rf_mean_power_filt = mean(rf_fft_filt_matrix, 1);
    
    %find all contraction segments (LH)
    lh_idx = find(tabella(:,2)==1);
    lh_start = tabella(lh_idx,4);
    lh_end = lh_start + round(0.8*fs);
    lh_fft = cell(length(lh_start),1);
    lh_fft_filt = cell(length(lh_start),1);
    %compute fft for each contraction and control if the frequency resolution
    %is always the same
    for i=1:length(lh_start)
        lh_signal = lh(:,lh_start(i):lh_end(i));
        lh_fft{i,1} = 1/fs*fft(lh_signal);   % apply the DFT operator  

        lh_signal_filt = lh_filtered(:,lh_start(i):lh_end(i));
        lh_fft_filt{i,1} = 1/fs*fft(lh_signal_filt);   % apply the DFT operator
        %lf_fft_filt{i,2} = fs/length(lf_fft_filt{i,1}); %frequency resolution    
    end
    %compute the mean of fft
    lh_fft_matrix = cell2mat(lh_fft);
    lh_mean_power = mean(lh_fft_matrix, 1);
    lh_fft_filt_matrix = cell2mat(lh_fft_filt);
    lh_mean_power_filt = mean(lh_fft_filt_matrix, 1);
    
    %find all contraction segments (RH)
    rh_idx = find(tabella(:,2)==4);
    rh_start = tabella(rh_idx,4);
    rh_end = rh_start + round(0.8*fs);
    rh_fft = cell(length(rh_start),1);
    rh_fft_filt = cell(length(rh_start),1);
    %compute fft for each contraction and control if the frequency resolution
    %is always the same
    for i=1:length(rh_start)
        rh_signal = rh(:,rh_start(i):rh_end(i));
        rh_fft{i,1} = 1/fs*fft(rh_signal);   % apply the DFT operator  

        rh_signal_filt = rh_filtered(:,rh_start(i):rh_end(i));
        rh_fft_filt{i,1} = 1/fs*fft(rh_signal_filt);   % apply the DFT operator
        %lf_fft_filt{i,2} = fs/length(lf_fft_filt{i,1}); %frequency resolution    
    end
    %compute the mean of fft
    rh_fft_matrix = cell2mat(rh_fft);
    rh_mean_power = mean(rh_fft_matrix, 1);
    rh_fft_filt_matrix = cell2mat(rh_fft_filt);
    rh_mean_power_filt = mean(rh_fft_filt_matrix, 1);    
    
    
    N_S = size(rh_mean_power,2);    % number of frequency samples
    F = fs/N_S;         % frequency resolution
    f = F*(0:N_S-1);    % frequency vector
    fig = figure(11);
    ax1 = subplot(4,1,1);
    plot(ax1, f,abs(lf_mean_power).^2, 'r') % plot the power spectrum
    hold(ax1,'on');
    plot(ax1, f,abs(lf_mean_power_filt).^2, 'b', 'LineWidth', 1.5) % plot the power spectrum
    hold(ax1,'off');
    grid on;                               
    axis([0 500 0 max(abs(lf_mean_power_filt).^2)])
    ylabel('P(f) [V^2/Hz]')
    title('LeftFoot');  
    
    ax2 = subplot(4,1,2);
    plot(ax2, f,abs(rf_mean_power).^2, 'r') % plot the power spectrum
    hold(ax2,'on');
    plot(ax2, f,abs(rf_mean_power_filt).^2, 'b', 'LineWidth', 1.5) % plot the power spectrum
    hold(ax2,'off');
    grid on;                               
    axis([0 500 0 max(abs(lf_mean_power_filt).^2)])
    ylabel('P(f) [V^2/Hz]')
    title('RightFoot');
    
    ax3 = subplot(4,1,3);
    plot(ax3, f,abs(lh_mean_power).^2, 'r') % plot the power spectrum
    hold(ax3,'on');
    plot(ax3, f,abs(lh_mean_power_filt).^2, 'b', 'LineWidth', 1.5) % plot the power spectrum
    hold(ax3,'off');
    grid on;                               
    axis([0 500 0 max(abs(lf_mean_power_filt).^2)])
    ylabel('P(f) [V^2/Hz]')
    title('LeftHand');
    
    ax4 = subplot(4,1,4);
    plot(ax4, f,abs(rh_mean_power).^2, 'r') % plot the power spectrum
    hold(ax4,'on');
    plot(ax4, f,abs(rh_mean_power_filt).^2, 'b', 'LineWidth', 1.5) % plot the power spectrum
    hold(ax4,'off');
    grid on;                               
    axis([0 500 0 max(abs(lf_mean_power_filt).^2)])
    xlabel('Frequency [Hz]')
    ylabel('P(f) [V^2/Hz]')
    title('RightHand');
    hL = legend('raw', 'filtered');
    newPosition = [0.81 0.87 0.05 0.05];
    newUnits = 'normalized';
    set(hL,'Position', newPosition, 'Units', newUnits);

    sgtitle('Mean power spectrum');
    name = [outpath + 'meanSpectrum.png'];
    saveas(fig, name)
    
    %% Save filtered signals
    %put signals in a matrix
    %LH, LF, RH, RF 
    EMGfiltered = [lh_filtered; lf_filtered; rh_filtered; rf_filtered];
    label = ["LH", "LF", "RH", "RF"];
    save([outpath+'/'+'EMGfiltered.mat'],'EMGfiltered', 'label', 'tabella');   
    
    %% Rectification
    % Full wave rectification of EMG
    lf_rectified = abs(lf_filtered);
    rf_rectified = abs(rf_filtered);
    lh_rectified = abs(lh_filtered);
    rh_rectified = abs(rh_filtered);

%     % Plot rectified signals
%     figure(12);
%     subplot(4,1,1)
%     plot(tmpdata.time{1}, lf_rectified); %LF
%     title('Left Foot');
%     grid on;
%     ylim([0, 2e-3]);
%     subplot(4,1,2)
%     plot(tmpdata.time{1}, rf_rectified); %RF
%     title('Right Foot');
%     grid on;
%     ylim([0, 2e-3]);
%     subplot(4,1,3)
%     plot(tmpdata.time{1}, lh_rectified); %LH
%     title('Left Hand');
%     grid on;
%     ylim([0, 2e-3]);
%     subplot(4,1,4)
%     plot(tmpdata.time{1}, rh_rectified); %RH
%     title('Right Hand');
%     grid on;
%     ylim([0, 2e-3]);
%     xlabel('Time [s]');
%     sgtitle('Rectified signals (time domain)')

    %% Smoothing - Linear envelope (6Hz lowpass)
    LP = 6;
    Wp = LP/NyqFreq;
    [F,E] = butter (4,Wp,'low');
    lf_linear = filtfilt(F,E,lf_rectified);
    rf_linear = filtfilt(F,E,rf_rectified);
    lh_linear = filtfilt(F,E,lh_rectified);
    rh_linear = filtfilt(F,E,rh_rectified);

%     %Plot smoothed signals
%     figure(14);
%     subplot(4,1,1)
%     plot(tmpdata.time{1}, lf_linear); %LF
%     title('Left Foot');
%     grid on;
%     ylim([0, 5e-4]);
%     subplot(4,1,2)
%     plot(tmpdata.time{1}, rf_linear); %RF
%     title('Right Foot');
%     grid on;
%     ylim([0, 5e-4]);
%     subplot(4,1,3)
%     plot(tmpdata.time{1}, lh_linear); %LH
%     title('Left Hand');
%     grid on;
%     ylim([0, 5e-4]);
%     subplot(4,1,4)
%     plot(tmpdata.time{1}, rh_linear); %RH
%     title('Right Hand');
%     grid on;
%     ylim([0, 5e-4]);
%     xlabel('Time [s]');
%     sgtitle('Smoothed signals (time domain)')

    %LF
    LF_linear = 1/fs*fft(lf_linear);   % apply the DFT operator
    N_S = length(LF_linear);    % number of frequency samples
    F = fs/N_S;         % frequency resolution
    f = F*(0:N_S-1);    % frequency vector
    %RF
    RF_linear = 1/fs*fft(rf_linear);   % apply the DFT operator
    %LH
    LH_linear = 1/fs*fft(lh_linear);   % apply the DFT operator
    %RH
    RH_linear = 1/fs*fft(rh_linear);   % apply the DFT operator

%     % Plot smoothed signal in the frequency domain
%     figure(16)
%     subplot(2,2,1)
%     plot(f,abs(LF_linear).^2, 'k', 'LineWidth', 1) % plot the power spectrum
%     grid on;                               
%     axis([0 10 0 max(abs(LF_linear).^2)])
%     xlabel('Frequency [Hz]')
%     ylabel('P(f)')
%     title('LF')
%     subplot(2,2,2)
%     plot(f,abs(RF_linear).^2, 'k', 'LineWidth', 1) % plot the power spectrum
%     grid on;                               
%     axis([0 10 0 max(abs(RF_linear).^2)])
%     xlabel('Frequency [Hz]')
%     ylabel('P(f)')
%     title('RF')
%     subplot(2,2,3)
%     plot(f,abs(LH_linear).^2, 'k', 'LineWidth', 1) % plot the power spectrum
%     grid on;                               
%     axis([0 10 0 max(abs(LH_linear).^2)])
%     xlabel('Frequency [Hz]')
%     ylabel('P(f)')
%     title('LH')
%     subplot(2,2,4)
%     plot(f,abs(RH_linear).^2, 'k', 'LineWidth', 1) % plot the power spectrum
%     grid on;                               
%     axis([0 10 0 max(abs(RH_linear).^2)])
%     xlabel('Frequency [Hz]')
%     ylabel('P(f)')
%     title('RH')
%     sgtitle('Smoothed signals (frequency domain)')
    
    %% Save smoothed signals
    %put signals in a matrix
    %LH, LF, RH, RF 
    EMGsmoothed = [lh_linear; lf_linear; rh_linear; rf_linear];
    label = ["LH", "LF", "RH", "RF"];
    save([outpath+'/'+'EMGsmoothed.mat'],'EMGsmoothed', 'label', 'tabella');   
end