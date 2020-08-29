function [data] = respiration(data, nSub, samplingRate, dataType, ... 
    lowPassBreathing, highPassBreathing)
	% receives data structure that contains raw respiration recordings (from Biopac MP36R system) and trigger locations
	% 
    % computes breathing continuous and binary phase,
    % as well as breathing rate and rate variability
    % continuous phase: [-pi, pi]
    % binary phase:  1 = inspiration = [-pi, 0)
    %                0 = expiration  = [0, pi]
    %
    % inputs should be:
    %       dataType = 'humanAirflow';
    %       samplingRate = 2000;
    %       lowPassBreathing = 0.2;
    %       highPassBreathing = 0.8;
	

    for i = 1:nSub
            % bandpass filter
            data(i).bp_respiration = ft_preproc_bandpassfilter(data(i).respiration, samplingRate, [lowPassBreathing highPassBreathing],2,'but');

            % compute phase using hilbert transform
            data(i).signalBreathingPhase = ft_preproc_hilbert(data(i).bp_respiration, 'angle');

            % get phase of respiration regarding the stimulus onset
            data(i).breathingContinuousPhase = data(i).signalBreathingPhase(data(i).trigger);
            data(i).breathingBinaryPhase = (data(i).breathingContinuousPhase < 0);

            % compute breathing rate and rate variability
            addpath C:\toolbox\breathmetrics-master
            addpath C:\toolbox\breathmetrics-master\breathmetrics_functions
            
            bm = breathmetrics(data(i).bp_respiration, samplingRate, dataType);
            bm.estimateAllFeatures('sliding', 0, 0);
            sfeat = getSecondaryRespiratoryFeatures(bm);
            data(i).breathingRate = sfeat('Breathing Rate');
            data(i).breathingRateVariability = sfeat('Coefficient of Variation of Breathing Rate');
    end
end