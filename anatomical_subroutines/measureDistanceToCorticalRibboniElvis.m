function [elecCoord, elecNames, isLeft, isSubdural, dist2CR] = measureDistanceToCorticalRibboniElvis(subj,cfg)

%function [avgCoord, elecNames, isLeft] = measureDistanceToCorticalRibboniElvis(subj,cfg)
%
% This function filters out electrodes that are farther than D mm from the 
% cortical ribbon (in the patient's brain).
% For subdural electrodes, it takes RAS "pial" coordinates (snapped 
% to the pial surface). For depth electrodes it uses the "CT" RAS
% coordinates.
%
% Input:
%   subj = FreeSurfer subject name
%
% Optional Inputs: passed as fields in a configuration structure
%   plotEm = 1 or 0.  If nonzero, a figure is created illustrating
%            the distribution of electrode distances. {default: 1}
%   elecCoord = N-by-3 numeric array with RAS electrode coordinates. 
%               {default: not used; the function looks into the subject's 
%               Freesurfer folder for electrode coordinate file instead}
%   elecNames = cell array of strings with electrode names, corresponding
%               to the rows of elecCoord. This argument is required 
%               if elecCoord is used. {default: not used; the function
%               looks into the subject's Freesurfer folder for electrode
%               name file instead}
%   isLeft    = N-dimensional binary vector where N is the # of electrodes.
%               0 indicates that an electrode is on/in the right hemisphere.
%               1 indicates a left hemisphere location. This argument is
%               required if elecCoord is used. Otherwise this information
%               is read from the participant's electrodeNames file.
%   isSubdural= N-dimensional binary vector where N is the # of electrodes.
%               0 indicates that an electrode is a depth electrode. 1
%               indicates a subdural electrode. This argument is only used
%               if elecNames is used. Otherwise this information is read
%               from the participant's electrodeNames file {default:
%               all electrodes are assumed to be subdural}
%   rmDepths = 1 or 0. If nonzero, depth electrodes are ignored. {default: 0}
%
% Outputs:
%   avgCoords = Electrode coordinates on FreeSurfer avg brain pial surface
%                (RAS coordinates)
%   elecNames = Channel names with the participant's name appended to the
%               beginning (e.g., PT001-Gd1)
%   isLeft    = N-dimensional binary vector where N is the # of electrodes.
%               0 indicates that an electrode is on/in the right hemisphere.
%               1 indicates a left hemisphere location.
%   dist   = distance from freesurfer's cortical ribbon
%
%
% Author:
% Yitzhak Norman, Malach Lab 2020
% (based on the code of: David Groppe, Mehtalab)
%

% parse input parameters in cfg structure and set defaults
if  ~isfield(cfg,'plotEm'),         plotEm = 1;     else    plotEm = cfg.plotEm;            end
if  ~isfield(cfg,'elecCoord'),      elecCoord = []; else    elecCoord = cfg.elecCoord;      end
if  ~isfield(cfg,'elecNames'),      elecNames = []; else    elecNames = cfg.elecNames;      end
if  ~isfield(cfg,'isLeft'),        isLeft = [];   else    isLeft = cfg.isLeft;      end
if  ~isfield(cfg,'isSubdural'),     isSubdural = [];   else    isSubdural = cfg.isSubdural;      end
if  ~isfield(cfg,'rmDepths'),       rmDepths = 0;   else    rmDepths = cfg.rmDepths;      end
checkCfg(cfg,'sub2AvgBrain.m');

% FreeSurfer Subject Directory
fsDir=getFsurfSubDir();
avgDir=fullfile(fsDir,'fsaverage');
subDir=fullfile(fsDir,subj);


if ~exist(avgDir,'dir')
    error('Folder for fsaverage is not present in FreeSurfer subjects directory (%s). Download it from here https://osf.io/qe7pz/ and add it.', ...
        fsDir);
end
if ~exist(subDir,'dir')
    error('Folder for %s is not present in FreeSurfer subjects directory (%s).',subj,fsDir);
end


% Take care of electrode names, hemisphere, and type
if isempty(elecNames)
    % Import electrode names
    elecFname=fullfile(subDir,'elec_recon',[subj '.electrodeNames']);
    elecInfo=csv2Cell(elecFname,' ',2);
    elecNames=elecInfo(:,1);
    nElec=size(elecInfo,1);
    isLeft=zeros(nElec,1);
    isSubdural=zeros(nElec,1);
    for a=1:nElec,
        if ~strcmpi(elecInfo{a,2},'D')
            isSubdural(a)=1;
        end
        if strcmpi(elecInfo{a,3},'L')
            isLeft(a)=1;
        end
    end
else
    nElec=length(elecNames);
    if isempty(isLeft)
        error('You need to specify cfg.isLeft when using cfg.elecNames');
    else
        if length(isLeft)~=nElec,
            error('elecNames and isLeft do not have the same # of electrodes.');
        end
    end
    if isempty(isSubdural),
        % assume all electrodes are subdural
        isSubdural=ones(nElec,1);
    else
        if length(isSubdural)~=nElec,
            error('isSubdural and isLeft do not have the same # of electrodes.');
        end
    end
end


% Take care of electrode coordinates in participant space
if isempty(elecCoord) % no electrode coordinates have been passed in the function call:
    % Import electrode PIAL coordinates
    coordFname=fullfile(subDir,'elec_recon',[subj '.PIAL']);
    coordCsv=csv2Cell(coordFname,' ',2);
    elecCoord=zeros(nElec,3);
    for a=1:nElec,
        for b=1:3,
            elecCoord(a,b)=str2double(coordCsv{a,b});
        end
    end
else
    if size(elecCoord,1)~=nElec,
        error('Electrode coordinates need to have the same number of rows as electrode names.');
    end
end


% Remove depths (optional)
if universalYes(rmDepths),
    keepIds=find(isSubdural);
    isLeft=isLeft(keepIds);
    elecNames=elecNames(keepIds);
    elecCoord=elecCoord(keepIds,:);
    nElec=length(elecNames);
    isSubdural=isSubdural(keepIds);
else
    % for depth electrodes: import "CT" coordinates in patient space:
    ptntCoordFile=fullfile(subDir,'elec_recon',[subj '.CT']);
    tempCsv=csv2Cell(ptntCoordFile,' ',2);
    nElec=size(tempCsv,1);
    ctCoords=zeros(nElec,3);
    for a=1:nElec
        for b=1:3
            ctCoords(a,b)=str2double(tempCsv{a,b});
        end
    end   
    
    elecNamesFile=fullfile(subDir,'elec_recon',[subj '.electrodeNames']);
    tempCsv=csv2Cell(elecNamesFile,' ',2);
    elecNamesCT=tempCsv(:,1);
    isLeftCT=zeros(nElec,1);
    isDepthCT=zeros(nElec,1);
    for a=1:nElec,
        if tempCsv{a,3}=='L'
            isLeftCT(a)=1;
        end
        if tempCsv{a,2}=='D'
            isDepthCT(a)=1;
        end
    end
  
    % sanity check:
    assert(all(isLeftCT==isLeft));
    assert(all(isDepthCT==~isSubdural));
    assert(all(strcmpi(elecNamesCT,elecNames)));
    
    % Select only depths
    elecCoord(isDepthCT==1,:)=ctCoords(isDepthCT==1,:);
    nElec=sum(isDepthCT);
    fprintf('\n %d depth electrodes were successfully processed... \n',nElec);
end


%read smoothwm and pial surfaces into matlab:
surftype = 'smoothwm';
[smoothwm.lh.vert, smoothwm.lh.tri] = read_surf(fullfile(subDir,'surf',['lh.' surftype]));
[smoothwm.rh.vert, smoothwm.rh.tri] = read_surf(fullfile(subDir,'surf',['rh.' surftype]));

surftype = 'pial-outer-smoothed';
[smoothpial.lh.vert, pial.lh.tri] = read_surf(fullfile(subDir,'surf',['lh.' surftype]));
[smoothpial.rh.vert, pial.rh.tri] = read_surf(fullfile(subDir,'surf',['rh.' surftype]));


dist2CR=nan(size(elecCoord,1),1);
for iElec=1:size(elecCoord,1)
    curCoord=elecCoord(iElec,:);
    switch isLeft(iElec)
        case 0
            d1=min(sqrt(sum(bsxfun(@minus,smoothwm.rh.vert,curCoord).^2,2)));
            d2=min(sqrt(sum(bsxfun(@minus,smoothpial.rh.vert,curCoord).^2,2)));
            dist2CR(iElec) = min(d1,d2);
        case 1
            d1=min(sqrt(sum(bsxfun(@minus,smoothwm.lh.vert,curCoord).^2,2)));
            d2=min(sqrt(sum(bsxfun(@minus,smoothpial.lh.vert,curCoord).^2,2)));
            dist2CR(iElec) = min(d1,d2);
    end    
end

if plotEm   
       
    figure('name','distance from CR','color','w','position',[0 0 1000 300]); 
    h1 = subplot(1,3,1); hold on;   
    % plot mesh LH   
    surf2plot = struct;
    surf2plot.faces = smoothwm.lh.tri; 
    surf2plot.vertices = smoothwm.lh.vert; 
    surf2plot = reducepatch(surf2plot,0.5);
    volSurfHp = patch('Vertices',surf2plot.vertices,'Faces',surf2plot.faces,...
        'facecolor',[0.5 0.5 0.5],'edgecolor','none','facealpha',0.25);
    box off
    axis tight square off     
  
    % plot mesh LH   
    surf2plot = struct;
    surf2plot.faces = smoothwm.rh.tri; 
    surf2plot.vertices = smoothwm.rh.vert; 
    surf2plot = reducepatch(surf2plot,0.5);
    volSurfHp = patch('Vertices',surf2plot.vertices,'Faces',surf2plot.faces,...
        'facecolor',[0.5 0.5 0.5],'edgecolor','none','facealpha',0.25);
    box off
    axis tight square off 
    
    view(180,10);
           
    % plot electrodes          
    for iElec=1:size(elecCoord,1)
            hold all;
            plot3(elecCoord(iElec,1),elecCoord(iElec,2),elecCoord(iElec,3),'o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor','k','MarkerSize',5);
            hold off;
    end
    
    subplot(1,3,2:3); 
    bar(1:length(dist2CR),dist2CR,0.75,'k'); hold on;
    axis tight;
    ylabel('Distance (mm)');
    xlabel('Electrode #'); 
    title('Distance from cortical ribbon');
end

























% 

% %%
% FSLUT=readFSLUT;
% % show all areas:
% [~,ind]=ismember(unique(vol),FSLUT.code);
% FSLUT.name(ind)
% 
% subcorticalArea={'Left-Cerebellum-Cortex','Left-Thalamus-Proper','Left-Caudate','Left-Putamen','Left-Pallidum','Brain-Stem','Left-Hippocampus','Left-Amygdala','Left-Accumbens-area'
% grayMatterAreas=FSLUT.code(ismember(FSLUT.name,{'Left-Cerebral-Cortex','Right-Cerebral-Cortex'}));
% grayMatterMask=ismember(vol,grayMatterAreas);
% 
% % convert to gray matter coordinate matrix
% [i j k]=ind2sub(size(vol),find(grayMatterMask));
% 
% % for each electrode, find the nearest gray matter
% 
% % computer distances between current electrode and each gray matter
% % coordinatee
% e
% % pick minimal distance
% 
% % store distance and gray matter kind
% 
% 
% 
