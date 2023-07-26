function [elecinfo, H1avg] = getHippocampalElectrodeCoord(channels, subjects, plot_single_subjects, useSaved)

% compute MNI coordinate for depth electrodes
% L1 = channel names
% L2 = subject IDs
% plotFlag = whether or not to plot 3D hippocampal map

narginchk(2,4);
if ~exist('useSaved','var'), useSaved=0; end
if ~exist('plot_single_subjects','var'), plot_single_subjects=0; end
%%

clc; warning('off');
% set the relevant path:
[parentfolder,path_to_toolboxes] = defineParentFolderAndPath;
% set the relevant path:
global path_to_toolboxes
global globalFsDir
globalFsDir = fullfile(parentfolder,'Freesurfer');
% Set Path:
addpath(fullfile(path_to_toolboxes,'eeglab14_1_2b'));
addpath(genpath(fullfile(parentfolder,'matlab_scripts')));
addpath(genpath(fullfile(path_to_toolboxes,'iELVis-master')));
addpath(genpath(fullfile(path_to_toolboxes,'afni_matlab\')));

if useSaved,  load(fullfile(globalFsDir,'hippocampalElectrodesInfo.mat'));
else     
    close all
    groupAvgCoords=[];
    groupNativeCoords=[];
    groupLabels=[];
    groupIsLeft=[];
    groupDist2sf = [];
    groupNumOfAdjVox = [];
    groupIsSFparc=[];

    
    % Compile smoothpatch's c-code functions (for plotting hippocampal anatomy):
    keepwd = pwd;
    cd(fullfile(parentfolder,'matlab_scripts','general_functions','smoothpatch_version1b'));
    eval(sprintf('mex %s %s','smoothpatch_curvature_double.c', '-compatibleArrayDims'))
    eval(sprintf('mex %s %s','smoothpatch_inversedistance_double.c', '-compatibleArrayDims'))
    eval(sprintf('mex %s %s','vertex_neighbours_double.c', '-compatibleArrayDims'))
    cd(keepwd);
    
    % channel list:
    tmp = split(channels,'_');
    L1 = tmp(:,1);
    L2 = tmp(:,2);
    tmp = split(L1,'-');
    L1 = tmp(:,1); clear tmp
    
    for iSub=1:length(subjects)
        subjid = subjects{iSub};
        currentCh = L1(strcmpi(L2,subjid));
        fprintf('Working on patient %s\nchannels: %s\n',subjid,cell2mat(join(currentCh,','))');
        maindir=fullfile(parentfolder, subjid);
        fid = fopen(fullfile(maindir,'subjectname'));
        C = textscan(fid,'%s'); fclose(fid);        
        FSsubjid = cell2mat(C{1}); % subject FS code
        cfg=[]; radius = 2; % radius (in mm) = maximal distance allowed from hippocampal subfiegetHippocampalElectrodeCoordlds CA1/2/3/SUB
        % compute distances from hippocampal subfield and plot all electrodes within "radius" mm from CA1/CA2/SUB:
        [elecCoord, elecNames, isLeft, isSubdural, isSFparc, dist2sf, numOfAdjVox, LUT, H1, H2, included_elec] = measureDistanceToHippocampalSubfields(FSsubjid,cfg,radius,subjid,plot_single_subjects);
        % Save figures:
        % set path:
        outdir=fullfile(globalFsDir,FSsubjid,'hippocampal_channels');
        if ~exist(outdir,'dir'), mkdir(outdir); end
        %=========================================================================
        if plot_single_subjects 
            for F = [H1 H2]
                figure(F);
                if F==H2
                    set_font_size_and_type;
                    export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-pdf','-painters','-nofontswap');
                end
                export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-jpg','-painters','-r300','-nofontswap')
            end
        end
        % =====================================================================
        selectedind = find(startsWith(elecNames,currentCh,'IgnoreCase',true) & endsWith(elecNames,currentCh,'IgnoreCase',true) ); 
        % use this if you want to exclude electrodes based on FS parcellation: 
        % selectedind = find(startsWith(elecNames,currentCh,'IgnoreCase',true) & endsWith(elecNames,currentCh,'IgnoreCase',true) & included_elec);
        % =====================================================================
        cfg.plotEm = 0;
        cfg.elecCoord = elecCoord(selectedind,:);
        cfg.elecNames = elecNames(selectedind);
        cfg.isLeft = isLeft(selectedind);
        cfg.isSubdural = ones(size(cfg.elecNames)); % make all "surface" elctrodes tagged "subdural"
        [avgCoords, elecNames, isLeft] = sub2AvgBrain(FSsubjid,cfg);
        groupAvgCoords = [groupAvgCoords; avgCoords];
        groupNativeCoords = [groupNativeCoords; elecCoord(selectedind,:)];
        groupLabels = [groupLabels; elecNames];
        groupIsLeft = [groupIsLeft; isLeft];
        groupDist2sf = [groupDist2sf; dist2sf(selectedind,:)];
        groupNumOfAdjVox = [groupNumOfAdjVox; numOfAdjVox(selectedind,:)];
        groupIsSFparc = [groupIsSFparc; isSFparc(selectedind,:)];
      
    end
    
    elecinfo=struct;
    elecinfo.groupAvgCoords = groupAvgCoords;
    elecinfo.groupNativeCoords = groupNativeCoords;
    elecinfo.groupLabels = groupLabels;
    elecinfo.groupIsSFparc = groupIsSFparc;
    elecinfo.groupIsLeft = logical(groupIsLeft);
    elecinfo.groupDist2sf = groupDist2sf;
    elecinfo.groupNumOfAdjVox = groupNumOfAdjVox;
    elecinfo.LUT = LUT;
    save(fullfile(globalFsDir,'hippocampalElectrodesInfo.mat'),'elecinfo');
    disp('data was saved')
    fprintf('\n filename: %s \n\n',fullfile(globalFsDir,'hippocampalElectrodesInfo.mat'));
end

%% Read fsaverage hippocampus into matlab:

HP = struct;
HP.lh = MRIread(fullfile(globalFsDir,'fsaverage','mri','aseg.mgz')); HP.lh.vol = HP.lh.vol==17;
HP.rh = MRIread(fullfile(globalFsDir,'fsaverage','mri','aseg.mgz')); HP.rh.vol = HP.rh.vol==53;
transform=struct;
transform.IJK2XYZ_tkreg=affine3d(HP.lh.tkrvox2ras'); % Define dicomIJK (CRS) -> fs mesh XYZ (tkReg) transform

set_figure_colors;
H1avg = figure('name','hippocampal channels','color','w','position',[0 0 500 200]);
H1avg.Tag='nopdf';
ind = strcmpi(elecinfo.LUT.subfield,'CA1');
vollh = HP.lh.vol;
volrh = HP.rh.vol;
assert(any(vollh(:))==1);
assert(any(volrh(:))==1);
% volCA1 = bsxfun(@plus, HP.lh.vol==LUT.code(ind), HP.rh.vol==LUT.code(ind));

% Use 3d smoothing:
smt = 3;
vollh = smooth3(vollh,'box',[smt smt smt]);
volrh = smooth3(volrh,'box',[smt smt smt]);
LHsurf = struct; RHsurf = struct;
[LHsurf.faces,LHsurf.vertices] = isosurface(vollh,0.5);
[RHsurf.faces,RHsurf.vertices] = isosurface(volrh,0.5);
LHsurf.vertices = transform.IJK2XYZ_tkreg.transformPointsForward(LHsurf.vertices);
RHsurf.vertices = transform.IJK2XYZ_tkreg.transformPointsForward(RHsurf.vertices);
%     LHsurf = smoothpatch(LHsurf);%,1,1,0.25);  % smooth the 3D surface
%     RHsurf = smoothpatch(RHsurf);%,1,1,0.25);  % smooth the 3D surface
str='tail';

% plot mesh 1:
h_lh = subplot(1,2,2); hold on;
title('left hippocampus')
surf2plot = struct;
surf2plot.faces = LHsurf.faces;
surf2plot.vertices = LHsurf.vertices;
volSurfHp = patch('Vertices',surf2plot.vertices,'Faces',surf2plot.faces,...
    'facecolor',COLOR.lightgray,'edgecolor','none','facealpha',0.1);
textpos = [mean(surf2plot.vertices(:,1)) min(surf2plot.vertices(:,2))-2 max(surf2plot.vertices(:,3))];
%text(textpos(1), textpos(2), textpos(3),str,'fontsize',10)
view(-110,25); axis tight equal off; box off; lights_on;
pos = get(gca,'position'); s1 = pos(3)*1.15; s2 = pos(4)*1.15;
set(gca,'position',[pos(1) pos(2) s1 s2])


% plot mesh 2:
h_rh = subplot(1,2,1); hold on;
title('right hippocampus')
surf2plot = struct;
surf2plot.faces = RHsurf.faces;
surf2plot.vertices = RHsurf.vertices;
volSurfHp = patch('Vertices',surf2plot.vertices,'Faces',surf2plot.faces,...
    'facecolor',COLOR.lightgray,'edgecolor','none','facealpha',0.1);

textpos = [mean(surf2plot.vertices(:,1)) min(surf2plot.vertices(:,2))-2 max(surf2plot.vertices(:,3))];
%text(textpos(1), textpos(2), textpos(3),str,'fontsize',10,'HorizontalAlignment','center')
view(110,25); axis tight equal off; box off; lights_on;
pos = get(gca,'position');  set(gca,'position',[pos(1) pos(2) s1 s2])

% plot electrodes
htext = [];
elecinfo.groupIsAnt = false(size(elecinfo.groupLabels));
D = [elecinfo.groupDist2sf.CA1'; elecinfo.groupDist2sf.CA2_3'; elecinfo.groupDist2sf.subiculum';];
% Electrode in 3d Mesh:
[xx,yy,zz]=sphere(20);
R=0.75; % sphere radius in mm
xx=xx'*R; yy=yy'*R; zz=zz'*R;
Helec = [];
pullingExtent = 10;
for iElec=1:size(elecinfo.groupAvgCoords,1)
    xyz = [elecinfo.groupAvgCoords(iElec,1),elecinfo.groupAvgCoords(iElec,2),elecinfo.groupAvgCoords(iElec,3)];
    if elecinfo.groupIsLeft(iElec), axes(h_lh); else axes(h_rh);  end
    d=[]; ind=[];
    [d(1),ind(1)]=min(sqrt(sum(bsxfun(@minus,LHsurf.vertices,xyz).^2,2)));
    [d(2),ind(2)]=min(sqrt(sum(bsxfun(@minus,RHsurf.vertices,xyz).^2,2)));
    [d,ii] = min(d);
    ind = ind(ii);
    % if the electrode is farther than 2 mm from the average MNI hipp, put it
    % closer to the mesh (for clearer visualization):
    if d>2
        disp(elecinfo.groupLabels{iElec});
        v = (xyz-LHsurf.vertices(ind,:))./norm(xyz-LHsurf.vertices(ind,:)); % norm direction (1mm length)
        if ii==1, xyz = LHsurf.vertices(ind,:)+v;
        else, xyz = RHsurf.vertices(ind,:)+v;
        end
    end     

    if xyz(2)<-22, col = COLOR.black;
    else, col = COLOR.black; elecinfo.groupIsAnt(iElec)=true; end
    Helec(iElec) = surf(xx+xyz(1),yy+xyz(2),zz+xyz(3),'facecolor',col,'edgecolor','none','FaceLighting','gouraud','SpecularStrength',0.2);
    set(Helec(iElec),'Tag',sprintf('%.4f,%.4f,%.4f,%d',xyz(1),xyz(2),xyz(3),ii));
%     xyz = pullElectrodesTowardsTheCamera(xyz, gca, pullingExtent);
%     Helec(iElec) =  plot3(xyz(1),xyz(2),xyz(3),'o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor','none','MarkerSize',5);
%     set(Helec(iElec),'Tag',elecinfo.groupLabels{iElec});
%     htext(end+1) = text(xyz(1),xyz(2),xyz(3),elecinfo.groupLabels{iElec},'fontsize',5);
end

elecinfo.Helec = Helec;
if false % change to true to add y=-22 plane (optional):    
    axes(h_rh); hold on;
    minx=min(RHsurf.vertices(:,1));
    maxx=max(RHsurf.vertices(:,1));
    minz=min(RHsurf.vertices(:,3));
    maxz=max(RHsurf.vertices(:,3));
    p = patch( [maxx minx minx maxx] , [-22 -22 -22 -22], [maxz maxz minz minz], COLOR.blue,'facealpha',0.1,'edgecolor','none');
    direction = [1 0 1];
    rotate(p,direction,-30)
end

end
