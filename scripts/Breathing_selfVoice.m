% psychotoolbox matlab file used to run the experiment described in
% "Breathing affects self-other voice discrimination in a bodily state with increased otherness"
% by Orepic, P., Park, H.D., Rognini, G., Faivre, N., and Blanke, O. 

function [Info] = Breathing_selfVoice(subjectID, pairedWithID, subject, task, condition)
% subjectID: identifier given to the subject
% pairedWithID: identifier of another subject with whom this subject's voice was morphed
% subject: which subject is the 'self': A or B, denoting the end of the self-other voice continuum between the 2 subjects
% 		 	as the same files are used, 15% self-voice of subject A is 85% self-voice of subject B
% task: "L" - loudness, "V" - voice, "N" - no task - questionnaire blocks in which participants passively listen to the voices from the self-other task
% condition: "S" - synchronous, or "A" - asynchronous
%
% example usage:
% Breathing_selfVoice('1', '2', 'A', 'L', 'S');

robotSeconds = 60;
arduinoBoard = 'Micro';
arduinoPort = 'COM9';

% configure arduino connection for triggers
ard = arduino(arduinoPort, arduinoBoard);
configurePin(ard,'D2', 'DigitalOutput');      % trial start
configurePin(ard,'D3', 'DigitalOutput');      % sounds
configurePin(ard,'D4', 'DigitalOutput');      % response

fileName = ['Breathing_selfVoice_S' num2str(subjectID) '_' subject '_' task '_' condition ]; % output .csv file

try 

    % ----------------------
    % Check input arguments:
    % ----------------------
    
    if ~(robotSeconds>3)
        disp('robotSeconds has to be bigger than 3! ');
        Info = robotSeconds; 
        return
    end

    if ~(subject == 'A' || subject == 'B')                          % the voice continuum rises or falls
		disp('Subject has to be either ''A'' or ''B''! ');
		Info = subject;
		return
    end    
    if ~(task == 'L' || task == 'V' || task == 'N')                  % loudness, voice or no task
		disp('Task has to be ''L'', ''V'' or ''N''! ');
		Info = task;
		return
    end    
	if ~(condition == 'S' || condition == 'A')                      % sync or async
		disp('Condition has to be either ''S'' or ''A''! ')
		Info = condition;
		return
	end
	
	
	% -----------------------
    % Parameters Modifiables:
    % -----------------------
	
	if (task == 'L')
		P.task              = [10 11 11.5 12.5 13 14];           % loudness (level) in dBFS from the quietest to the loudest 
    else                                                         % (10 = -14 dBFS, 12 = -12 dBFS, 14 = -10 dBFS)
        P.task              = [15 30 45 55 70 85];               % percentage of subject's A voice (It is reversed for the subject B!)
	end
    P.word                  = [1 2 3 4 5 6 7 8 9 10];            % sound files



	% ------
    % Clicks:
    % ------
	% check GetMouse() for details
	leftClick 			= 1;								
	rightClick 			= 3;
	
   
    % ----------------------------------------
    % Randomization of the trial presentation: 
    % ----------------------------------------
    % Define trials
    Info   = struct;
    Info.SubjectID = subjectID;
    Info.PairedWithID = pairedWithID;
    Info.subject = subject;
	Info.condition = condition;
    Info.task = task;    
    Info.T = [];

    %-------------------------------------------------------------
    % Info.This loop determines the locations of relevant stimuli:
    %-------------------------------------------------------------  
    itrial = 0;
    for i_task = 1:length(P.task)             
        for i_word = 1:length(P.word)
            itrial = itrial + 1;
            % Which task?
            Info.T(itrial).task = P.task(i_task);
            % which word
            Info.T(itrial).word = P.word(i_word);
        end            
    end


    % Now randomize the trials by shuffling the order.
    nt = length(Info.T);
    Info.T(1:nt)= Info.T(Shuffle(1:nt));
    

    %-------------------------------------------------------------
    %                  Start The Experiment:                     
    %-------------------------------------------------------------
    disp('Experiment starts in 3 seconds! ');
    WaitSecs(3);                                                                    	% just so that it is not too sudden
    disp('Experiment started! ');
    disp(['Manipulate the robot for ' num2str(robotSeconds) ' seconds! ']);
    [beepSound, beepFreq] = MakeBeep(500,0.5); 
    sound(beepSound, beepFreq);
    WaitSecs(robotSeconds-3);                                                                    	% start with the robot
    disp('Sounds start in 3 seconds! '); 
    sound(beepSound, beepFreq);
    WaitSecs(3);   
   
		for  t=1:length(Info.T)        
			disp(['Start of trial ' num2str(t) '.']);
            % set trigger for the trial start
            writeDigitalPin(ard, 'D2', 1);
            WaitSecs(0.01);
            writeDigitalPin(ard, 'D2', 0);
            
            %---------------------------------------------------------
            % Creating jitter between trials and words within a trial:
            %---------------------------------------------------------
            P.timeBetweenWords					= 0.5;                                  % fixed delay between the reference and target voice
            P.ITI         						= 1 + 0.5*rand(1);  					% intertrial inteval jitter between 1.0 and 1.5 seconds
            Info.T(t).ITI 						= P.ITI;              					% save the ITI used
			Info.T(t).timeBetweenWords 			= P.timeBetweenWords;              		% save the timeBetweenWords used

			
            %-----------------------------
            % Open the corresponding file:
            %-----------------------------			
			folder = 'voice';
			if task == 'L'
				folder = 'loudness';
			end			
			
            path = [ num2str(Info.T(t).word) filesep folder filesep num2str(Info.T(t).task) '.wav'];
            
            [soundData, soundFreq] = audioread(path);
			soundData = soundData';
			
			
			%---------------------------------------
			% Open the corresponding reference file:
			%---------------------------------------
			ref = 50;
			if task == 'L'
				ref = 12;
			end		
	
			pathReference = [num2str(Info.T(t).word) filesep folder filesep num2str(ref) '.wav'];
			
			[referenceData, referenceFreq] = audioread(pathReference);
			referenceData = referenceData';
			
			
			%-------------------
            % Play the Sounds:
            %-------------------
			
			sound(referenceData, referenceFreq); 				
            writeDigitalPin(ard, 'D3', 1);                      % the sound trigger goes up after the first sound
            
            WaitSecs(P.timeBetweenWords); 
            
			sound(soundData, soundFreq);
            writeDigitalPin(ard, 'D3', 0);                      % the sound trigger goes down after the second sound
			
			%------------------
			% Get the response:
			%------------------
             
			responseTime = 0;
			buttons = 0;
			startTime = GetSecs;								
            if ~(task == 'N')                  	
                while ~buttons										% waiting for the click
                    [~, ~, buttons] = GetMouse();
                    responseTime = GetSecs - startTime;                    
                end
                
                writeDigitalPin(ard, 'D4', 1);                  	% the trigger for response
                WaitSecs(0.01);
                writeDigitalPin(ard, 'D4', 0);
                
                Info.T(t).responseTime = responseTime;
                buttonsIndex = find(buttons);                       % in case the subject presses both
                Info.T(t).response = buttonsIndex(1);
            else
                Info.T(t).responseTime = 0;                     	% in questionnaire blocks, where participants passively listen to voices 
                Info.T(t).response = 0;								% there is no response
            end
            

			%----------------------------------------------------- 
			% Calculate if response was correct or not & feedback:
			%-----------------------------------------------------
			
			% the reference sound is always played first, therefore corresponds always to the left click
			% the played sound is the second one to be played, therefore corresponds always to the right click
			% the question is always "Which one do you perceive bigger (louder or 'more you' respectively)?"
			%
			% lower than reference means the reference is bigger and therefore left button should be pressed
			% bigger than reference means the right click should be pressed
			%
            % for the subject A the voice percentages are reversed
            % (e.g. Info.T(t).task = 30 means 30% for the subject B and 70% for subject A)
            %
            %    task | 15  |  30  |  45  |  50  | 55  |  70  |  85  | 
            %       A | 85  |  70  |  55  |  50  | 45  |  30  |  15  |
            % correct |  R  |   R  |   R  |   x  |  L  |   L  |   L  |
            %       B | 15  |  30  |  45  |  50  | 55  |  70  |  85  |
            % correct |  L  |   L  |   L  |   x  |  R  |   R  |   R  |
			
            % disp([num2str(ref) ' , ' num2str(Info.T(t).task) ' , '  num2str(Info.T(t).response) ]);
            Info.T(t).ResponseCorrect = 0;                                              % for the (task == N) case
            if (task == 'L' || (task == 'V' && subject == 'B'))
                if (Info.T(t).task < ref) && (Info.T(t).response == leftClick)
                    Info.T(t).ResponseCorrect = 1;
                    
                elseif (Info.T(t).task > ref) && (Info.T(t).response == rightClick)
                    Info.T(t).ResponseCorrect = 1;
				
                else 
                	Info.T(t).ResponseCorrect = 0;				
                end
                
            elseif (task == 'V' && subject == 'A')    
                 if (Info.T(t).task < ref) && (Info.T(t).response == rightClick)
                    Info.T(t).ResponseCorrect = 1;
			
                elseif (Info.T(t).task > ref) && (Info.T(t).response == leftClick)
                    Info.T(t).ResponseCorrect = 1;
				
                else 
                	Info.T(t).ResponseCorrect = 0;				
                end               
            end
			
            if (task == 'N')
                WaitSecs(Info.T(t).ITI + 1.2);            % wait a bit longer to compensate the response time
            else
                WaitSecs(Info.T(t).ITI);
            end
            
            			         

		end

    disp('Experiment ended! ');
    save(fileName, 'Info', 'P');
	disp(['File: ' fileName '.mat saved!']);


catch ME
    
    rethrow(ME);

end

%------------------
% Create csv files:
%------------------  

% the heard stimuli is adjusted to subject B, so for subject A it
% should be reversed (e.g. when heard 30, it is 30% of subject B, 70% of A)
% the columns of the data should therefore be swapped for each
% word-percentage pair (only in task V and subject A)
%
% subject | 15  |  30  |  45  |  50  | 55  |  70  |  85  | 
%       A | 85  |  70  |  55  |  50  | 45  |  30  |  15  |
%       B | 15  |  30  |  45  |  50  | 55  |  70  |  85  | 

% extract the data
disp('Creating the csv file!');
load([fileName '.mat']);
d = cell2mat(struct2cell(Info.T));
d=reshape(d,[7,60]);
d = d';

% dataOld = d;
% reverse 15-85, 30-70, 45-55 for subject A
if (Info.task == 'V' && Info.subject == 'A')
   [rows, ~] = size(d);
   
   for w=1:10 
       % find the corresponding rows (word - percentage)
       row15 = d(d(:,2)==w & d(:, 1)==15, :);
       row85 = d(d(:,2)==w & d(:, 1)==85, :);
       
       row30 = d(d(:,2)==w & d(:, 1)==30, :);
       row70 = d(d(:,2)==w & d(:, 1)==70, :);
       
       row45 = d(d(:,2)==w & d(:, 1)==45, :);
       row55 = d(d(:,2)==w & d(:, 1)==55, :);
       
       i15 = 0;
       i85 = 0;
       i30 = 0;
       i70 = 0;
       i45 = 0;
       i55 = 0;    
       % find the corresponding  row indices 
       for i=1:rows
           if size(find(d(i, :) == row15), 2) == size(row15, 2)
               i15 = i;
           end
           if size(find(d(i, :) == row85), 2) == size(row85, 2)
               i85 = i;
           end
           
           if size(find(d(i, :) == row30), 2) == size(row30, 2)
               i30 = i;
           end
           if size(find(d(i, :) == row70), 2) == size(row70, 2)
               i70 = i;
           end
           
           if size(find(d(i, :) == row45), 2) == size(row45, 2)
               i45 = i;
           end
           if size(find(d(i, :) == row55), 2) == size(row55, 2)
               i55 = i;
           end
       end
       % swap the columns
       d([i15 i85], 3:7) = d([i85 i15], 3:7);
       d([i30 i70], 3:7) = d([i70 i30], 3:7);
       d([i45 i55], 3:7) = d([i55 i45], 3:7);
   end
   
end

subjectID = Info.SubjectID;
pairedWithID = Info.PairedWithID;

% code subject order
if (Info.subject == 'A')
    subjectFirst = 1;
else
    subjectFirst = 2;
end

% code task type
if (Info.task == 'L')
    taskID = 1;
elseif (Info.task == 'V')
    taskID = 2;
else
    taskID = 3;
end

% code condition
if (Info.condition == 'A')
    conditionID = 1;
else
    conditionID = 2;
end

% create columns
subject = ones(60, 1)*subjectID;
pairedWith = ones(60, 1)*pairedWithID;
subjectOrder = ones(60, 1)*subjectFirst;
task = ones(60, 1)*taskID;
condition = ones(60, 1)*conditionID;
word = d(:, 2);
taskValue = d(:, 1);
response = d(:, 6);
responseCorrect = d(:, 7);
iti = d(:, 3);
responseTime = d(:, 5);

% create output matrix
data = horzcat(subject, pairedWith, subjectOrder, task, condition, word, taskValue, response, responseCorrect, iti, responseTime);

csvwrite([ fileName '.csv'], data);
disp(['File:' fileName '.csv saved!']);

end

