 function [elecCoord, elecNames, isLeft, isSubdural, isSFparc, dist2sf, numOfAdjVox, LUT, H1, H2, included_elec] = measureDistanceToHippocampalSubfields(subj,cfg,radius,subjid,plot_single_subjects)

% function [elecCoord, elecNames, isLeft, isSubdural, isSFparc, dist2sf, numOfAdjVox, LUT, H1, H2] = measureDistanceToHippocampalSubfields(subj,cfg,radius,subjid,plot_single_subjects)
%
% This function measures the distance of electrodes to each hippocampal
% subfield.
%
% Author:
% Yitzhak Norman, Malach Lab 2020
% (based on the code of: David Groppe, Mehtalab)

% parse input parameters in cfg structure and set defaults
if  ~isfield(cfg,'plotEm'),         plotEm = 1;     else    plotEm = cfg.plotEm;            end
if  ~isfield(cfg,'elecCoord'),      elecCoord = []; else    elecCoord = cfg.elecCoord;      end
if  ~isfield(cfg,'elecNames'),      elecNames = []; else    elecNames = cfg.elecNames;      end
if  ~isfield(cfg,'isLeft'),        isLeft = [];   else    isLeft = cfg.isLeft;      end
if  ~isfield(cfg,'isSubdural'),     isSubdural = [];   else    isSubdural = cfg.isSubdural;      end
if  ~isfield(cfg,'rmDepths'),       rmDepths = 0;   else    rmDepths = cfg.rmDepths;      end

if ~exist('plot_single_subjects','var'), plot_single_subjects=0; end

dist2sf = [];
numOfAdjVox = [];
H1 = []; H2 = [];
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


% read hippocampal subfields into matlab:
HP = struct; 
lh_subfield_flag = 0; rh_subfield_flag = 0;
filename1 = dir(fullfile(subDir,'mri','lh.hippoSfLabels-T1*FSvoxelSpace.mgz'));
if ~isempty(filename1)
    HP.lh = MRIread(fullfile(filename1.folder,filename1.name)); lh_subfield_flag = 1;
else 
    fprintf('\n *** subfields file is missing!!! \n');  
    HP.lh = MRIread(fullfile(subDir,'mri','aseg.mgz')); HP.lh.vol = HP.lh.vol==17;
end 
filename2 = dir(fullfile(subDir,'mri','rh.hippoSfLabels-T1*FSvoxelSpace.mgz'));
if ~isempty(filename2)    
    HP.rh = MRIread(fullfile(filename2.folder,filename2.name)); rh_subfield_flag = 1;
else
    fprintf('\n *** subfields file is missing!!! \n');
    HP.rh = MRIread(fullfile(subDir,'mri','aseg.mgz')); HP.rh.vol = HP.rh.vol==53;
end
transform=struct;
% Define dicomIJK (CRS) -> dicomXYZ transform
transform.IJK2XYZ=affine3d(HP.lh.vox2ras');
% Define dicomIJK (CRS) -> fs mesh XYZ (tkReg) transform
transform.IJK2XYZ_tkreg=affine3d(HP.lh.tkrvox2ras');

% sf = unique(HP.lh.vol(HP.lh.vol~=0));
sf = [203:206, 208:210, 226];
[c, ll, ct]=read_fscolorlut(fullfile(fsDir,'FreeSurferColorLUT.txt'));

% load and edit look-up table (LUT):
fn = {}; colortable = []; code = [];
for k = 1:length(sf)
    tmp = ll(c==sf(k),:);
    fn{k,1} = tmp(isstrprop(tmp,'alphanum'));
    if strcmpi(fn{k,1},'CA3'),fn{k,1}='CA2_3'; end
    code(k,1) = sf(k);
    colortable(k,:) = [ct(c==sf(k),1:3)]./255;
end
LUT = table(fn,code,colortable,'VariableNames',{'subfield','code','color'});

for k=1:size(LUT,1)
    if ~rh_subfield_flag && lh_subfield_flag        
        binaryVolume = permute(bsxfun(@plus, HP.lh.vol==LUT.code(k), HP.rh.vol==1),[2,1,3]);
    elseif rh_subfield_flag && ~lh_subfield_flag        
        binaryVolume = permute(bsxfun(@plus, HP.lh.vol==1, HP.rh.vol==LUT.code(k)),[2,1,3]);
    elseif  ~rh_subfield_flag && ~lh_subfield_flag        
        binaryVolume = permute(bsxfun(@plus, HP.lh.vol==1, HP.rh.vol==1),[2,1,3]);
    else
        binaryVolume = permute(bsxfun(@plus, HP.lh.vol==LUT.code(k), HP.rh.vol==LUT.code(k)),[2,1,3]);
    end
    [numOfAdjVox(k,:),dist2sf(k,:)] = distToBinaryVolume(elecCoord,binaryVolume,transform,radius);
    fprintf('\n%s\n',LUT.subfield{k});
end
fprintf('\n');
isSFparc = true(size(dist2sf,2),1);
isSFparc(isLeftCT==1)=logical(lh_subfield_flag);
isSFparc(isLeftCT==0)=logical(rh_subfield_flag);
dist2sf = array2table(dist2sf','RowNames',strcat(subjid,'_',subj,'-',elecNames),'VariableNames',LUT.subfield);
numOfAdjVox = array2table(numOfAdjVox','RowNames',strcat(subjid,'_',subj,'-',elecNames),'VariableNames',LUT.subfield);
included_elec = dist2sf.CA1<=radius | dist2sf.CA2_3<=radius | dist2sf.subiculum<=radius;

if plot_single_subjects

    H1 = figure('name','hippocampal channels no labels','color','w','position',[0 0 500 200]);
    ind = strcmpi(LUT.subfield,'CA1');  
    vollh = HP.lh.vol;
    volrh = HP.rh.vol;
    assert(any(vollh(:))==1);
    assert(any(volrh(:))==1);
    
    % Use 3d smoothing:
    LHsurf = struct; RHsurf = struct;
    [LHsurf.faces,LHsurf.vertices] = isosurface(vollh,0.5);
    [RHsurf.faces,RHsurf.vertices] = isosurface(volrh,0.5);
    LHsurf = reducepatch(LHsurf,0.75);   % downsample the 3D surface and then smooth
    RHsurf = reducepatch(RHsurf,0.75);   % downsample the 3D surface
    try
        fprintf('\n\n *** trying smoothing *** \n\n');
        LHsurf = smoothpatch(LHsurf,0,1);  % smooth the 3D surface (tested only on MATLAB R2018b)
        RHsurf = smoothpatch(RHsurf,0,1);  % smooth the 3D surface (tested only on MATLAB R2018b)
    catch % if smoothing didn't work
        fprintf('\n\n *** skipping smoothing *** \n\n');
        LHsurf = reducepatch(LHsurf,0.75);   % downsample the 3D surface and then smooth
        RHsurf = reducepatch(RHsurf,0.75);   % downsample the 3D surface
        LHsurf = smoothpatch(LHsurf,0,1);  % smooth the 3D surface (tested only on MATLAB R2018b)
        RHsurf = smoothpatch(RHsurf,0,1);  % smooth the 3D surface (tested only on MATLAB R2018b)
    end
    % plot mesh 1:
    h_lh = subplot(1,2,2); hold on;
    
    %title('left hippocampus')
    surf2plot = struct;
    surf2plot.faces = LHsurf.faces;
    surf2plot.vertices = transform.IJK2XYZ_tkreg.transformPointsForward(LHsurf.vertices);
    volSurfHp = patch('Vertices',surf2plot.vertices,'Faces',surf2plot.faces,...
        'facecolor',[0.9 0.9 0.9],'edgecolor','none','facealpha',0.5);    
    textpos = [mean(surf2plot.vertices(:,1)) min(surf2plot.vertices(:,2))-2 max(surf2plot.vertices(:,3))];   

    view(-110,25);
    axis tight equal off; box off; lights_on; 
    
    pos = get(gca,'position'); s1 = pos(3)*1.15; s2 = pos(4)*1.15;
    set(gca,'position',[pos(1) pos(2) s1 s2])
       
    % plot mesh 2:
    h_rh = subplot(1,2,1); hold on;

    %title('right hippocampus')
    surf2plot = struct;
    surf2plot.faces = RHsurf.faces;
    surf2plot.vertices = transform.IJK2XYZ_tkreg.transformPointsForward(RHsurf.vertices);
    volSurfHp = patch('Vertices',surf2plot.vertices,'Faces',surf2plot.faces,...
        'facecolor',[0.9 0.9 0.9],'edgecolor','none','facealpha',0.5);
    view(110,25);
    textpos = [mean(surf2plot.vertices(:,1)) min(surf2plot.vertices(:,2))-2 max(surf2plot.vertices(:,3))];   
    axis tight equal off; box off; lights_on; 
    pos = get(gca,'position');  set(gca,'position',[pos(1) pos(2) s1 s2])

    % plot the electrodes:
    hipch = []; htext = []; 
    
    % Electrode in 3d Mesh:
    [xx,yy,zz]=sphere(20);
    R=1; % sphere radius in mm
    xx=xx'*R; yy=yy'*R; zz=zz'*R;
    roi = [205,206,207,208]; % CA1/CA2/CA3/SUB
    pullingExtent = 50;
   
    for iElec=1:size(elecCoord,1)
        xyz = [elecCoord(iElec,1),elecCoord(iElec,2),elecCoord(iElec,3)];        
        if dist2sf.CA1(iElec)<=radius || dist2sf.CA2_3(iElec)<=radius || dist2sf.subiculum(iElec)<=radius
         
            if isLeft(iElec) axes(h_lh); else axes(h_rh);  end                        
            h = surf(xx+xyz(1),yy+xyz(2),zz+xyz(3),'facecolor',[1,0,0],'edgecolor','none','FaceLighting','gouraud','SpecularStrength',0.2); 
                       %h = plot3(xyz(1),xyz(2),xyz(3),'o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor','none','MarkerSize',5); uistack(h,'top');
            h.Tag = elecNames{iElec}; 
            hipch = cat(1,hipch,iElec);   
            % add electrode names (optional):
%             xyz_pooled = pullElectrodesTowardsTheCamera(xyz, gca, pullingExtent); % pool coordinates toward the camera to add fully visible text            
%             htext(end+1) = text(xyz_pooled(1),xyz_pooled(2),xyz_pooled(3), ['   ' elecNames{iElec}],'fontsize',5,'VerticalAlignment','top','fontweight','bold');            
%             uistack(htext(end),'top');       
        end
    end
    
    % Bar chart of electrode distance to different subfields:
    H2 = figure('name','distance to subfields','color','w','position',[0 0 350 120]);
    D = [dist2sf.CA1(hipch)'; dist2sf.CA2_3(hipch)'; dist2sf.subiculum(hipch)';];
    bcolor = [LUT.color(strcmpi(LUT.subfield,'CA1'),:); LUT.color(strcmpi(LUT.subfield,'CA2_3'),:); LUT.color(strcmpi(LUT.subfield,'subiculum'),:)];
    n = length(hipch);
    xlabels = {elecNames{hipch}};
    if size(D,2)==1
        D = [D nan(3,1)];
        bcolor = repmat(bcolor,[1,1,2]);    
        n = n+1;
        xlabels = {elecNames{hipch},''};
    end    
    bcolor = permute(repmat(bcolor,[1,1,length(hipch)]),[3,1,2]);
    h=superbar([1:n]',D','BarFaceColor',bcolor); hold on;
    set(gca,'XTickLabelRotation',45,'Xtick',[1:n],'XtickLabel',{elecNames{hipch}},...
        'xlim',[0.25 max(n,10)+0.75],'ylim',[0 10],'ytick',[0:2:10]);
  
    ylabel('Distance (mm)');       
    
end