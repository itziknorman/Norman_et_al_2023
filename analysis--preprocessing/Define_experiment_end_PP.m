% Define the borders of the entire experimental duration in each subject:
% (from file onset until the end of the relevant task)

switch subjid
    case 'PP01'
          [experiment_end, ~] = time2indices({'14:14:39'},{'14:48:25'},0,EEG);
          warning('Patient excluded. hippocampal lesion.');
    case 'PP02'
          [experiment_end, ~] = time2indices({'16:32:35'},{'17:02:40'},0,EEG);
    case 'PP03'
          [experiment_end, ~] = time2indices({'13:42:13'},{'14:14:10'},0,EEG); 
    case 'PP04'
          [experiment_end, ~] = time2indices({'08:47:12'},{'09:34:15'},0,EEG);             
    case 'PP05'
          [experiment_end, ~] = time2indices({'14:28:37'},{'14:54:40'},0,EEG);
    case 'PP06'
          [experiment_end, ~] = time2indices({'10:37:40'},{'11:14:30'},0,EEG);        
    case 'PP07'
          [experiment_end, ~] = time2indices({'15:14:23'},{'15:43:50'},0,EEG);  
    case 'PP08'
          [experiment_end, ~] = time2indices({'16:08:42'},{'16:48:40'},0,EEG);      
    case 'PP09'
          [experiment_end, ~] = time2indices({'08:23:07'},{'09:03:20'},0,EEG);
    case 'PP10'
          [experiment_end, ~] = time2indices({'15:46:00'},{'16:23:00'},0,EEG);
    case 'PP11'
          [experiment_end, ~] = time2indices({'08:46:16'},{'09:16:40'},0,EEG);
    case 'PP12'
         warning('Patient excluded. no hippocampal electrodes.');
    case 'PP13'
         [experiment_end, ~] = time2indices({'16:16:40'},{'16:54:00'},0,EEG);
         warning('Patient excluded. no hippocampal electrodes.');
    case 'PP14'
        [experiment_end, ~] = time2indices({'15:04:05'},{'15:34:00'},0,EEG);
    case 'PP15'
        [experiment_end, ~] = time2indices({'14:43:18'},{'15:14:45'},0,EEG);
    case 'PP16'
        [experiment_end, ~] = time2indices({'16:10:56'},{'16:46:20'},0,EEG);
    case 'PP17'
        [experiment_end, ~] = time2indices({'15:30:18'},{'16:06:45'},0,EEG);
    case 'PP18'
        [experiment_end, ~] = time2indices({'14:06:37'},{'14:34:50'},0,EEG);
    case 'PP19'
        [experiment_end, ~] = time2indices({'11:49:37'},{'12:25:00'},0,EEG);
        warning('Patient excluded. Seizure at the end of the memory test');
end

