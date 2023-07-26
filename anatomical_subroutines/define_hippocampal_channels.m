function [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag,select_by_anatomy_flag)
% Defining the Hippocampus Channel: (based on anatomy)

if nargin<2, ref_flag = 2; end
if nargin<3, select_by_anatomy_flag = 1; end

method = 1;
% ref_flag: 1 = common-reference (CREF); 2 = bipolar montage (BP)
switch subjid
    case 'PP01' % excluded
        if ref_flag==1
            hippocampus ='TB5'; % selected CA1 electrode
            hippocampus_all_channels = {'TB4','TB5','TB6'};
            WM_ref={'TB8'};
        elseif ref_flag==2
            hippocampus = 'TB5-TB8';
            hippocampus_all_channels = {'TB4-TB8','TB5-TB8'};
        end
    case 'PP02'
        if ref_flag==1
            hippocampus ='HA6'; % selected CA1 electrode
            hippocampus_all_channels = {'HA1','HA2','HA3','HA4','HA5','HA6','TB1'};
            WM_ref={'HA8','TB2'}; 
        elseif ref_flag==2
            hippocampus = 'HA6-HA8';
            hippocampus_all_channels = {'HA1-HA8','HA2-HA8','HA3-HA8','HA4-HA8','HA5-HA8','HA6-HA8','TB1-TB2'};
            
        end
    case 'PP03'
        if ref_flag==1
            hippocampus ='NA2'; % selected CA1 electrode
            hippocampus_all_channels = {'NA2','NA1','HM3','HM2','HM1'};
            WM_ref={'NA9','HM8'};
        elseif ref_flag==2
            hippocampus = 'NA2-NA9';
            hippocampus_all_channels = {'NA2-NA9','NA1-NA9','HM3-HM8','HM2-HM8','HM1-HM8'};
        end
        
    case 'PP04'
        if ref_flag==1
            hippocampus ='HA3'; % selected CA1 electrode
            hippocampus_all_channels = {'HP2','HA2','HA3','HA4'};
            WM_ref={'HA8','HP5'}; 
        elseif ref_flag==2
            hippocampus = 'HA3-HA8';
            hippocampus_all_channels = {'HP2-HP5','HA2-HA8','HA3-HA8','HA4-HA8'};
        end
        
    case 'PP05'
        if ref_flag==1
            hippocampus ='HM4'; % selected CA1 electrode
            hippocampus_all_channels = {'HM1','HM2','HM4'}; % 'HM3' was excluded during preproc because of noise
            WM_ref={'HM8'};
        elseif ref_flag==2
            hippocampus = 'HM4-HM8'; 
            hippocampus_all_channels = {'HM1-HM8','HM2-HM8','HM4-HM8'};
        end
        
    case 'PP06'
        if ref_flag==1
            hippocampus ='HA3'; % selected CA1 electrode
            hippocampus_all_channels = {'HA3','HA2','HP4','HP3','HP2'};
            WM_ref={'HP11','HA10'};
        elseif ref_flag==2
            hippocampus = 'HA3-HA10'; 
            hippocampus_all_channels = {'HA3-HA10','HA2-HA10','HP4-HP11','HP3-HP11','HP2-HP11'};
        end
        
    case 'PP07'
        if ref_flag==1
            hippocampus ='HA5'; % selected CA1 electrode (HA4 as well)
            WM_ref={'HA8','NA8','HP5'};            
            hippocampus_all_channels = {'HP3','HP2','HP1','NA5','NA4','NA3','NA2','NA1','HA6','HA5','HA4','HA3','HA2','HA1'}; 
        elseif ref_flag==2
            hippocampus = 'HA5-HA8'; 
            hippocampus_all_channels = {'HP3-HP5','HP2-HP5','HP1-HP5','NA5-NA8','NA4-NA8','NA3-NA8','NA2-NA8','NA1-NA8',...
                                        'HA6-HA8','HA5-HA8','HA4-HA8','HA3-HA8','HA2-HA8','HA1-HA8'};
            
        end
    case 'PP08'
        if ref_flag==1
            hippocampus ='HM3';
            hippocampus_all_channels = {'HM3'}; % not connected: 'HM5','HM4'
            WM_ref={'TS4'}; 
        elseif ref_flag==2
            hippocampus = 'HM3-TS4';
            hippocampus_all_channels = {'HM3-TS4'};
        end
    case 'PP09'
        if ref_flag==1
            hippocampus ='HA3';            
            hippocampus_all_channels = {'HP3','HP2','TB1','HA4','HA3','HA2','HA1'}; % lesion in the amygdala and the uncus
            WM_ref={'HA7','HP5','TB5'};
        elseif ref_flag==2
            hippocampus = 'HA3-HA7';
            hippocampus_all_channels = {'HP3-HP5','HP2-HP5','TB1-TB5','HA4-HA7','HA3-HA7','HA2-HA7','HA1-HA7'};
        end
    case 'PP10' 
        if ref_flag==1
            hippocampus ='HM4';
            hippocampus_all_channels = {'HM2','HM3','HM4','TB3','TB2'};
            WM_ref={'HM6','TB4'};
        elseif ref_flag==2
            hippocampus = 'HM4-HM6';
            hippocampus_all_channels = {'HM2-HM6', 'HM3-HM6', 'HM4-HM6','TB3-TB4','TB2-TB4'};
        end
    case 'PP11'
        if ref_flag==1
            hippocampus ='HA4';
            hippocampus_all_channels = { 'HA5','HA4','HA3','HA2','HA1','NA1'};         
            WM_ref={'HA9','NA6'}; 
        elseif ref_flag==2
            hippocampus = 'HA4-HA9';
            hippocampus_all_channels = {'HA5-HA9','HA4-HA9','HA3-HA9','HA2-HA9','HA1-HA9','NA1-NA6'};
        end
    case 'PP14'
        if ref_flag==1
            hippocampus ='HA3';
            hippocampus_all_channels = {'HA4','HA3','HA2','HP6','HP5','HP4','HP3','HP2','TO2'};
            WM_ref={'HA8','HP8','TO8'};
        elseif ref_flag==2
            hippocampus = 'HA3-HA8';
            hippocampus_all_channels = {'HA2-HA8','HA3-HA8','HA4-HA8','HP2-HP8','HP3-HP8','HP4-HP8','HP6-HP8','HP5-HP8','TO2-TO8'};
        end
        
    case 'PP15'
        if ref_flag==1
            hippocampus ='HM4';
            hippocampus_all_channels = {'HM1','HM2', 'HM3', 'HM4','HM5'};
            WM_ref={'HM7'};
        elseif ref_flag==2
            hippocampus = 'HM4-HM7';
            hippocampus_all_channels = {'HM1-HM7','HM2-HM7','HM3-HM7','HM4-HM7','HM5-HM7'};
        end
        
    case 'PP16'
        if ref_flag==1
            hippocampus ='HM4';
            hippocampus_all_channels = {'HM2', 'HM3', 'HM4'}; % not connected: 'HM5'
            WM_ref={'HM7'};
        elseif ref_flag==2
            hippocampus = 'HM4-HM7';
            hippocampus_all_channels = {'HM2-HM7','HM3-HM7','HM4-HM7'};
        end
        
    case 'PP17'
        if ref_flag==1
            hippocampus ='HA1'; % HA3 is closer to CA1, but shows abnormal activity
            hippocampus_all_channels = {'HA1', 'HA2', 'HA3'};
            WM_ref={'HA7'};
        elseif ref_flag==2
            hippocampus = 'HA1-HA7';  % HA3-HA7 is closer to CA1, but shows abnormal activity
            hippocampus_all_channels = {'HA1-HA7', 'HA2-HA7', 'HA3-HA7'};
        end
        
    case 'PP18'
        if ref_flag==1
            hippocampus ='HA4';
            hippocampus_all_channels = {'TB4','TB3','TB2','HA5','HA4','HA3','HA2','HA1','TI4','TI3','TI2'};
            WM_ref={'HA7','TB7','TI9'};
        elseif ref_flag==2
            hippocampus = 'HA4-HA7';
            hippocampus_all_channels = {'TB4-TB7', 'TB3-TB7', 'TB2-TB7', 'HA5-HA7','HA4-HA7','HA3-HA7','HA2-HA7','HA1-HA7','TI4-TI9','TI3-TI9','TI2-TI9'};
        end
        
    case 'SUB01ASSUTA'
        if ref_flag==1
            hippocampus ='HAR3';
            hippocampus_all_channels = {'PHL2','PHL1','HAL2','HAL1','HAR3','HAR2','HAR1'};
            WM_ref={'HAL6','HAR7','PHL7'};
        elseif ref_flag==2
            hippocampus = 'HAR3-HAR7';
            hippocampus_all_channels = {'PHL2-PHL7','PHL1-PHL7','HAL2-HAL6','HAL1-HAL6','HAR3-HAR7','HAR2-HAR7','HAR1-HAR7'};
        end
end

if ~exist('WM_ref','var'); WM_ref={}; end
if select_by_anatomy_flag, hippocampus = select_electrode_by_anatomy(subjid,method,hippocampus_all_channels); end
end

function selected_electrode = select_electrode_by_anatomy(subjid,method,hippocampus_all_channels)
% The function returns the electrode closest to CA1 based on freesurfer segmentation

    fprintf('\n SELECTING HIPPOCAMPAL ELECTRODE BY ANATOMY');
    
    % set the relevant path:
    [filepath,~] = fileparts(matlab.desktop.editor.getActiveFilename); addpath(fullfile(filepath,'..'));
    [parentfolder,path_to_toolboxes] = defineParentFolderAndPath;

    FsDir = fullfile(parentfolder,'Freesurfer');
    
    % process hippocampal electrodes:
    L = join([hippocampus_all_channels; repmat({subjid},size(hippocampus_all_channels))],'_',1);
    [elecinfo,Hmap] = plotHippocampalElectrodesSingleSubject(L,{subjid},0,0);

    %load(fullfile(FsDir,'hippocampalElectrodesInfo.mat'));
    % Find the best CA1 electrode in this subjects:
    subjdir=fullfile(parentfolder, subjid);
    fid = fopen(fullfile(subjdir,'subjectname')); C = textscan(fid,'%s'); fclose(fid);
    FSsubjid = cell2mat(C{1}); % subject FS code
    fprintf('\n *** processing subject %s (%s) *** \n',subjid, FSsubjid);
    currentSubElec = find(contains(elecinfo.groupLabels ,FSsubjid));
    switch method
        case 1 % CA1 closest electrode
            [~,ind] = min(elecinfo.groupDist2sf.CA1(currentSubElec,:));
        case 2 % best wrapped electrode
            [~,ind] = max(elecinfo.groupNumOfAdjVox.CA1(currentSubElec,:));
    end
    elecAnatLabel = elecinfo.groupLabels{currentSubElec(ind)};  tmp = split(elecAnatLabel,'-');
    selected_electrode = hippocampus_all_channels{contains(hippocampus_all_channels,tmp(2))};
    fprintf('-CA1 electrode: %s (%s) \n',selected_electrode,elecAnatLabel);
end


