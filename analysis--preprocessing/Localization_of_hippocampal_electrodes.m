%% Locating the hippocampal electrodes based on distance form CA1/CA2/CA3/Subiculum

clear all
close all
clc;

% set the relevant path:
[filepath,~] = fileparts(matlab.desktop.editor.getActiveFilename); addpath(fullfile(filepath,'..'));
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath;

clc; warning('off');
% set the relevant path:
global path_to_toolboxes
global globalFsDir
globalFsDir = fullfile(parentfolder,'Freesurfer');

% Set Path:
addpath(fullfile(path_to_toolboxes,'eeglab2021.1'));
rmpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));
addpath(genpath(fullfile(path_to_toolboxes,'iELVis-master')));
addpath(genpath(fullfile(path_to_toolboxes,'afni_matlab')));

[ALLEEG, EEG, CURRENTSET] = eeglab;
warning('off')

subjects={'PP02','PP03','PP04','PP05','PP06','PP07','PP08','PP09','PP10','PP11','PP15','PP16','PP17','PP18'};

plot_single_subjects = 1;

 for iSub = 1:length(subjects)
        subjid = subjects{iSub};

        fprintf('Working on patient %s\n',subjid);
        maindir=fullfile(parentfolder, subjid);
        fid = fopen(fullfile(maindir,'subjectname'));
        C = textscan(fid,'%s'); fclose(fid);        
        FSsubjid = cell2mat(C{1}); % subject FS code
        if isempty(dir(fullfile(globalFsDir,FSsubjid,'elec_recon','*.LEPTO')))
            makeIniLocTxtFile(FSsubjid)
            yangWangElecPjct(FSsubjid);           
        end
        cfg=[]; radius = 2;
        
        % compute distances from hippocampal subfield and plot all electrodes within "radius" mm from CA1/CA2/CA3/SUB:        
        [elecCoord, elecNames, isLeft, isSubdural, isSFparc, dist2sf, numOfAdjVox, LUT, H1, H2] = measureDistanceToHippocampalSubfields(FSsubjid,cfg,radius,subjid,plot_single_subjects);
       
        % Save:
        % set figure folder:
        outdir=fullfile(globalFsDir,FSsubjid,'hippocampal_channels');
        if ~exist(outdir,'dir'), mkdir(outdir); end
        %=========================================================================
        if plot_single_subjects 
            for F = [H1 H2]
                figure(F);
                if F==H2
                    set_font_size_and_type;
                    export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-pdf','-painters','-nofontswap');
                else
                    export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-pdf','-opengl','-nofontswap');
                end
                export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-jpg','-painters','-r300','-nofontswap')
            end
        end
 end

