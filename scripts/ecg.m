function [data] = ecg(data, nSub, minPeakHeight, ... 
    minPeakDistance, lowPassFr, highPassFr, samplingRate)
	% receives data structure that contains raw heartbeat recordings (from Biopac MP36R system) and trigger locations
	%
    % computes heartbeat continuous and binary phase,
    % as well as heartbeat rate and rate variability
    % continuous phase: [0, 2*pi]
    % binary phase:  based on Kunzendorf, et al., "Active information sampling varies across the cardiac cycle"
	%
    % inputs are: minimum R-peak values, min distance between R-peaks and bandpass filter frequencies.
    % They should be:
    %     minPeakHeight = 0.23;
    %     minPeakDistance = 800; %samples
    %     samplingRate = 2000;
    %     lowPassFr = 1;
    %     highPassFr = 40;
    
    for s = 1:nSub
		% bandpass filter
        data(s).bp_heartbeat = ft_preproc_bandpassfilter(data(s).heartbeat, samplingRate, [lowPassFr highPassFr],2,'but');


        % Find the peaks of the ecg signal
        [qrspeaks,peak_locations] = findpeaks(data(s).bp_heartbeat,'MinPeakHeight',minPeakHeight,...
            'MinPeakDistance',minPeakDistance);

        trigger_index = data(s).trigger;

        continuousPhase = [];
        binaryPhase = [];
		timeAfterRpeak = [];


        for i = 1:length(trigger_index)
            % for all the intervals between 2 peaks.
            for j = 1:length(peak_locations)-1
                % If the trigger is between them, compute the phase
                if trigger_index(i) >= peak_locations(j) && trigger_index(i) < peak_locations(j+1)
                    peakT1 = peak_locations(j);
                    peakT2 = peak_locations(j+1);

                    % Compute the phase at the trigger
                    continuousPhase(i) = 2*pi*((trigger_index(i)-peakT1)/(peakT2-peakT1));

                    % actual time after R peak (added later)
                    timeAfterRpeak(i) = (trigger_index(i) - peakT1) / samplingRate; 


                    binaryPhase(i) = 0; % neither systole nor diastole

                    % -------
                    % systole 
                    % -------
                    % begining:  53.0–82.4 ms (M = 65.1, SD = 6.25) after the R peak
                    % end: 296–375  ms (M = 329, SD = 22.2) after the R peak


                    sys_beg_time = 0.065; 
                    sys_end_time = 0.329;

                    sys_beg = sys_beg_time * samplingRate + peakT1;
                    sys_end = sys_end_time * samplingRate + peakT1;


                    % --------
                    % diastole
                    % ---------
                    % beginning: 346–425 ms (M = 379, SD = 22.2) after the R peak 
                    % end: 652–1101 ms (M = 805, SD = 129) after the R peak

                    dia_beg_time = 0.379; 
                    dia_end_time = 0.805;

                    dia_beg = dia_beg_time * samplingRate + peakT1;
                    dia_end = dia_end_time * samplingRate + peakT1;

                    if trigger_index(i) >= sys_beg && trigger_index(i) <= sys_end
                        binaryPhase2(i) = 1; % systole 
                    elseif trigger_index(i) >= dia_beg && trigger_index(i) <= dia_end
                        binaryPhase2(i) = 2; % diastole 
                    end
                    
                end
            end
        end

        data(s).ecgContinuousPhase = continuousPhase;
        data(s).ecgBinaryPhase = binaryPhase;
        data(s).timeAfterRpeak = timeAfterRpeak;

        % heart rate
        duration = length(data(s).heartbeat)/samplingRate;
        nr_peaks = length(qrspeaks);
        data(s).ecgRate = nr_peaks/duration*60; %beats per minute
        
        % heartrate variability
		addpath(genpath('ECG_preprocess/biosig4octmat-3.0.7'));
		detect_struct = qrsdetect(data(s).heartbeat,samplingRate,2);
		EVENT = detect_struct.EVENT;
		r_onsets = EVENT.POS(EVENT.TYP==hex2dec('0501'))/EVENT.SampleRate;
		r_on = EVENT.POS(EVENT.TYP==hex2dec('0501'));
		r_diffs = diff(r_onsets);
		temp_on=r_onsets(1);
		r_onsets(1)=[];
		ibi = [r_onsets, r_diffs];

		% remove ectopic ibi with HRVAS methods
		outliers = locateOutliers(ibi(:,1),ibi(:,2),'sd',3);
		ibi(:,2)=replaceOutliers(ibi(:,1),ibi(:,2),outliers,'spline',nan);
			
		new_r_onsets = cumsum([temp_on;ibi(:,2)]);
		EVENT.POS(EVENT.TYP==hex2dec('0501')) = new_r_onsets*EVENT.SampleRate;
		detect_struct.EVENT=EVENT;
		HRV_struct = heartratevariability(detect_struct);
		data(s).ecgRateVariability = HRV_struct.RMSSD;

		rmpath(genpath('ECG_preprocess/biosig4octmat-3.0.7'));
    end 
end