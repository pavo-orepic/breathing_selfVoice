function [data] = ecg(data, nSub, minPeakHeight, ... 
    minPeakDistance, lowPassFr, highPassFr, samplingRate)
    % Receives a row of a subject's data and computes continuous phase of
    % heartbeat [0, 2Pi] and 2 versions of the binary phase, based on 2
    % different papers - Hyeong (binaryPhase1) and Kunzendorf
    % (binaryPhase2). Computes heartRate and heartRatevariability as well.
    % binary phase: 1 = systole, 2 = diastole, 0 = neither
    % 
    % inputs are: mean peak values, distance between peaks and bandpass
    % filter frequencies. They should be:
    %     minPeakHeight = 0.23;
    %     minPeakDistance = 800; %samples
    %     samplingRate = 2000;
    %     lowPassFr = 1;
    %     highPassFr = 40;
    
    % bandpass filter
    d =  designfilt('bandpassiir','FilterOrder',8, ...
    'HalfPowerFrequency1',lowPassFr,'HalfPowerFrequency2',highPassFr, ...
    'SampleRate',samplingRate);

    for s = 1:nSub
        %ecg_data = dataSubject.heartbeat;
        data(s).bp_heartbeat = ft_preproc_bandpassfilter(data(s).heartbeat, samplingRate, [lowPassFr highPassFr],2,'but');


        % Find the peaks of the ecg signal
        [qrspeaks,peak_locations] = findpeaks(data(s).bp_heartbeat,'MinPeakHeight',minPeakHeight,...
            'MinPeakDistance',minPeakDistance);

        trigger_index = data(s).T2_index;
        %click_index = data(s).T3_index;

        continuousPhase = [];
        binaryPhase1 = [];
        binaryPhase2 = [];
        %click_phase = [];

        for i = 1:length(trigger_index)
            % for all he intervals between 2 peaks.
            for j = 1:length(peak_locations)-1
                % If the trigger is between them, compute the phase
                if trigger_index(i) >= peak_locations(j) && trigger_index(i) < peak_locations(j+1)
                    peakT1 = peak_locations(j);
                    peakT2 = peak_locations(j+1);

                    % Compute the phase at the trigger
                    continuousPhase(i) = 2*pi*((trigger_index(i)-peakT1)/(peakT2-peakT1));

                    % actual time after R peak (added later)
                    timeAfterRpeak(i) = (trigger_index(i) - peakT1) / samplingRate; 


                    % ----------------%
                    %  Hyeong's paper %
                    % ----------------%
                    binaryPhase1(i) = 0; % neither systole nor diastole

                    % -------
                    % systole 
                    % -------
                    % Systole: —200 before R peak to +50ms after R peak 

                    sys_beg_time1 = 0; 
                    sys_end_time1 = 0.050;

                    sys_beg1 = sys_beg_time1 * samplingRate + peakT1;
                    sys_end1 = sys_end_time1 * samplingRate + peakT1;

                    sys_beg_time2 = 0.2; 
                    %sys_end_time2 = 1;

                    sys_beg2 = peakT2 - sys_beg_time2 * samplingRate;
                    sys_end2 = peakT2;


                    % --------
                    % diastole
                    % ---------
                    % Diastole:  +150 ms to +450 ms after R peak


                    dia_beg_time = 0.150; 
                    dia_end_time = 0.450;

                    dia_beg = dia_beg_time * samplingRate + peakT1;
                    dia_end = dia_end_time * samplingRate + peakT1;

                    if trigger_index(i) >= sys_beg1 && trigger_index(i) <= sys_end1
                        binaryPhase1(i) = 1; % systole 
                    elseif trigger_index(i) >= sys_beg2 && trigger_index(i) <= sys_end2
                        binaryPhase1(i) = 1; % systole 
                    elseif trigger_index(i) >= dia_beg && trigger_index(i) <= dia_end
                        binaryPhase1(i) = 2; % diastole 
                    end 


                    % --------------------%
                    %  Kunzendorf's paper %
                    % --------------------%
                    binaryPhase2(i) = 0; % neither systole nor diastole

                    % -------
                    % systole 
                    % -------
                    % begining:  53.0–82.4 ms (M = 65.1, SD = 6.25) after the R peak
                    % end: 296–375  ms (M = 329, SD = 22.2) after the R peak


                    sys_beg_time_2 = 0.065; 
                    sys_end_time_2 = 0.329;

                    sys_beg_2 = sys_beg_time_2 * samplingRate + peakT1;
                    sys_end_2 = sys_end_time_2 * samplingRate + peakT1;


                    % --------
                    % diastole
                    % ---------
                    % beginning: 346–425 ms (M = 379, SD = 22.2) after the R peak 
                    % end: 652–1101 ms (M = 805, SD = 129) after the R peak

                    dia_beg_time_2 = 0.379; 
                    dia_end_time_2 = 0.805;

                    dia_beg_2 = dia_beg_time_2 * samplingRate + peakT1;
                    dia_end_2 = dia_end_time_2 * samplingRate + peakT1;

                    if trigger_index(i) >= sys_beg_2 && trigger_index(i) <= sys_end_2
                        binaryPhase2(i) = 1; % systole 
                    elseif trigger_index(i) >= dia_beg_2 && trigger_index(i) <= dia_end_2
                        binaryPhase2(i) = 2; % diastole 
                    end
                    
                end
            end
        end

        %data(s).bp_heartbeat = ecg_data_bpf_1_40;
        data(s).ecgContinuousPhase = continuousPhase;
        data(s).ecgBinaryPhase1 = binaryPhase1;
        data(s).ecgBinaryPhase2 = binaryPhase2;
        data(s).timeAfterRpeak = timeAfterRpeak;

    %     % heart rate
    %     duration = length(data(s).heartbeat)/samplingRate;
    %     nr_peaks = length(qrspeaks);
    %     data(s).ecgRate = nr_peaks/duration*60; %beats per minute
    %     
    %     % heartrate variability
    %     % (Implementation by Arthur Barakat, a student who did a project with me)
    %     addpath D:\student_lab_immersion\heartbeat\BiopacAnalysis
    %     addpath D:\student_lab_immersion\heartbeat\BiopacAnalysis\ECG_preprocess\biosig4octmat-3.0.7
    %     addpath D:\student_lab_immersion\heartbeat\BiopacAnalysis\ECG_preprocess\biosig4octmat-3.0.7\biosig\t300_FeatureExtraction
    %     addpath D:\student_lab_immersion\heartbeat\BiopacAnalysis\ECG_preprocess\biosig4octmat-3.0.7\biosig\t250_ArtifactPreProcessingQualityControl
    %     addpath D:\student_lab_immersion\heartbeat\BiopacAnalysis\HRV_Stuff\HRVAS
    %     addpath D:\student_lab_immersion\heartbeat\BiopacAnalysis\ECG_preprocess\biosig4octmat-3.0.7\biosig\t200_FileAccess
    %     addpath(genpath('ECG_preprocess/biosig4octmat-3.0.7'));
    %     addpath D:\student_lab_immersion\heartbeat\BiopacAnalysis\ECG_preprocess\biosig4octmat-3.0.7\NaN\inst
    %     addpath D:\student_lab_immersion\heartbeat\BiopacAnalysis\ECG_preprocess\biosig4octmat-3.0.7\tsa\inst
    %     
    %     % Compute the HRV. RMSSD is used generally as a score for HRV
    %     % NOTE : the biosig_HRV function is modified to add the biopack toolbox
    %     % path and remove it after use to prevent shadowing of other functions.
    %     [HRV_struct, ~,~,~]=biosig_HRV(data(s).heartbeat,samplingRate);
    %     data(s).ecgRateVariability = HRV_struct.RMSSD;
    end 
end