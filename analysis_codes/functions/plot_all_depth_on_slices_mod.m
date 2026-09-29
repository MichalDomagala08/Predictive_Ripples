function plot_all_depth_on_slices_mod(fsSub,elecInfoType,cfg)

% function plotAllDepthsOnSlices(fsSub,elecInfoType,cfg)
% Minimalna modyfikacja: overlay pokazuje TYLKO wskazany region (np. hippocampus).
% Inne elementy pozostają bez zmian.

% Creates a figure illustrating the location of each depth electrode contact 
% in a sagittal, coronal, and axial slice and indicates which part of
% the brain it is in. You need to first need to record electrode
% coordinates in iELVis conventions (e.g., using yangWangElecPjct.m or
% dykstraElecPjct.m), have an mgrid file of electrode information from
% BioImageSuite, or have analogous info in MNI space.
%
% Required Inputs:
%  fsSub - Patient's freesurfer directory name
%  elecInfoType - {'mgrid', 'BIDS-iEEG', or 'MNI'} a string indicating the
%                 format electrode information is stored in. 
%
% Optional cfg parameters:
%  mgridFname - mgrid filename and path. If empty, and
%               elecInfoType=='mgrid', we assume the mgrid file is in the
%               subject's elec_recon subfolder and named *.mgrid, where *
%               is the subject's FreeSurfer ID.
%  fullTitle  - If 1, the mgrid and mri voxel coordinates are displayed in
%               the figure title along with the electrode name and anatomical
%               location. {default: 0}
%  markerSize - The size of the dot in each slice used to represent an
%               electrode's location. {default: 30}
%  cntrst     - 0< number <=1 The lower this number, the lower the brightness
%               of the image (i.e., the lower the voxel value corresponding to
%               white). {default: 0.5}
%  anatOverlay -If 1, color is overlayed on the brain to show FreeSurfer's
%              automatic segmentation of brain areas (neocortex uses
%              Desikan-Killiany parcellation). Alternatively define the fullpath
%              to another parcellation file. {default: 0}
%  colorLUT    - fullpath to color lookup table if you would like to use
%                non default colors for your parcellation overlay.
%                {default: FreeSurferColorLUTnoFormat.txt}
%  pauseOn   - If 1, Matlab pauses after each figure is made and waits for
%              a keypress. {default: 0}
%  figOverwrite - If non-zero, only one new Matlab figure will be produced 
%              when this function is called. Each time a new electrode
%              needs to be visualized the figure is cleared. This is useful
%              when lots of depths have been used and you're printing the
%              figures. {default: 1}
%  printFigs - 1 or directory. If 1, each figure is output to a jpg file in the patient's
%              elec_recon/PICS folder and the figure is closed after the
%              jpg is created. This is particularly useful for implants with
%              a large number of depth contacts. If a directory, the jpg files
%              get saved there instead. {default: 0}
%  bidsDir   - Full path to an BIDS-iEEG root directory. If specified,
%              electrode location and pial surface files will be
%              imported from this directory. If not specified,
%              these data will be imported from the subject's FreeSurfer directory.
%              Note that all electrodes are colored red when you use this
%              option as BIDS-iEEG does not associate a default color with
%              each electrode.
%  bidsSes   - integer. The BIDS-iEEG session number. This has no effect if
%              bidsDir not specified {default: 1}
%
%
% Examples:
%  %Specify mgrid file and do NOT print
%  cfg=[];
%  cfg.mgridFname='/Applications/freesurfer/subjects/PT001/elec_recon/PT001.mgrid';
%  plotAllDepthsOnSlices('PT001','mgrid',cfg);
%
%  %Use FreeSurfer file structure and print
%  cfg=[];
%  cfg.printFigs=1;
%  plotAllDepthsOnSlices('PT001','mgrid',cfg);
%
%
% Author: David M. Groppe
% Feb. 2015
% Feinstein Institute for Medical Research/Univ. of Toronto





% ---- default cfg handling ----
if ~exist('cfg','var'), cfg = struct(); end
if ~isfield(cfg,'mgridFname'),    mgridFname=[];    else mgridFname=cfg.mgridFname; end
if ~isfield(cfg,'fullTitle'),     fullTitle=0;      else fullTitle=cfg.fullTitle; end
if ~isfield(cfg,'markerSize'),    markerSize=30;    else markerSize=cfg.markerSize; end
if ~isfield(cfg,'cntrst'),        cntrst=.5;       else cntrst=cfg.cntrst; end

% correct anatOverlay handling
if ~isfield(cfg,'anatOverlay')
    anatOverlay = 0;
else
    anatOverlay = cfg.anatOverlay;
end

if ~isfield(cfg,'colorLUT'),    colorLUT=0;          else colorLUT=cfg.colorLUT; end
if ~isfield(cfg,'pauseOn'),    pauseOn=0;          else pauseOn=cfg.pauseOn; end
if ~isfield(cfg,'printFigs'),    printFigs=0;          else printFigs=cfg.printFigs; end
if ~isfield(cfg, 'bidsDir'),      bidsDir=[];         else bidsDir=cfg.bidsDir; end
if ~isfield(cfg, 'bidsSes'),      bidsSes=1;         else bidsSes=cfg.bidsSes; end
if ~isfield(cfg, 'figOverwrite'),  figOverwrite=1;  else figOverwrite=cfg.figOverwrite; end

% region defaults (user can set cfg.region and cfg.regionColor)
region_spec = [];
if isfield(cfg,'region') && ~isempty(cfg.region)
    region_spec = cfg.region;
end
if isfield(cfg,'regionColor') && ~isempty(cfg.regionColor)
    region_color = cfg.regionColor;
else
    region_color = [0.4940, 0.1840, 0.5560];
end

% FreeSurfer Subject Directory
if exist('cfg','var') && isfield(cfg,'fsurfSubDir') && ~isempty(cfg.fsurfSubDir)
    fsdir = cfg.fsurfSubDir;
else
    fsdir = getFsurfSubDir();
end

% minimal cfg cleaning for downstream functions
if isfield(cfg,'fsurfSubDir'), cfg = rmfield(cfg,'fsurfSubDir'); end

% sanity checks
if (nargin<3),
   error('plotAllDepthsOnSlices requires 2 arguments'); 
end
validFormats={'mgrid', 'BIDS-iEEG', 'MNI'};
if isempty(findStrInCell(elecInfoType,validFormats))
    error('%s is not a valid value for elecInfoType parameter',elecInfoType);
end
if strcmpi(elecInfoType,'BIDS-iEEG') && isempty(bidsDir)
    error('You need to specifiy cfg.bidsDir when reading electrode information in BIDS-iEEG format.'); 
end

% Define path to neuroimaging files
if isempty(bidsDir)
    mriDir=fullfile(fsdir,fsSub,'mri');
else
    mriDir=fullfile(bidsDir,'derivatives','iELVis',['sub-' fsSub],'mri');
end

% Load MRI (brainmask)
mriFname=fullfile(mriDir,'brainmask.mgz');
if ~exist(mriFname,'file')
    error('File %s not found.',mriFname);
end
mri=MRIread(mriFname);
mx=max(max(max(mri.vol)))*cntrst;
mn=min(min(min(mri.vol)));
sVol=size(mri.vol);

% Get electrode info
if strcmpi(elecInfoType,'BIDS-iEEG'),
    coordFname=fullfile(bidsDir,['sub-' fsSub],'ieeg',sprintf('sub-%s_ses-%.2d_space-postimplant_electrodes.tsv',fsSub,bidsSes));
    fprintf('Taking electrode info from %s.\n',coordFname);
    elecCoordCsv=csv2Cell(coordFname,9,0); %9=tab
    nElecTotal=size(elecCoordCsv,1)-1;
    RAS_coor=zeros(nElecTotal,3);
    coordHdrs={'x','y','z'};
    for csvLoopB=1:3,
        colId=findStrInCell(coordHdrs{csvLoopB},elecCoordCsv(1,:),1);
        for csvLoopA=1:nElecTotal,
            RAS_coor(csvLoopA,csvLoopB)=str2double(elecCoordCsv{csvLoopA+1,colId});
        end
    end
    % convert RAS coordinates to LIP
    elecMatrix=zeros(nElecTotal,3);
    elecMatrix(:,1)=-RAS_coor(:,1)+129;
    elecMatrix(:,2)=-RAS_coor(:,3)+129;
    elecMatrix(:,3)=-RAS_coor(:,2)+129;
    
    elecLabels=cell(nElecTotal,3);
    nameId=findStrInCell('name',elecCoordCsv(1,:),1);
    typeId=findStrInCell('type',elecCoordCsv(1,:),1);
    hemId=findStrInCell('hemisphere',elecCoordCsv(1,:),1);
    elecRgb=zeros(nElecTotal,3); elecRgb(:,1)=1;
    for csvLoopA=1:nElecTotal,
        tempName=elecCoordCsv{csvLoopA+1,nameId};
        switch elecCoordCsv{csvLoopA+1,typeId}
            case 'grid', tempType='G';
            case 'strip', tempType='S';
            case 'depth', tempType='D';
            otherwise, error('Unrecognized electrode type: %s', elecCoordCsv{csvLoopA+1,typeId});
        end
        tempHem=elecCoordCsv{csvLoopA+1,hemId};
        elecLabels{csvLoopA}=sprintf('%s%s_%s',tempHem,tempType,tempName);
    end
elseif strcmpi(elecInfoType,'MNI'),
    elecReconDir=fullfile(fsdir,fsSub,'elec_recon');
    mniInfoFname=fullfile(elecReconDir,'persystElecInfo.tsv');
    mniPairsFname=fullfile(elecReconDir,'persystElecPairs.tsv');
    if ~exist(mniInfoFname,'file') || ~exist(mniPairsFname,'file')
        error('Missing MNI electrode files in elec_recon.');
    end
    [elecMatrix, elecLabels, elecRgb]=mni2Matlab(fsSub);
else
    if isempty(mgridFname)
        [elecMatrix, elecLabels, elecRgb]=ut_mgrid2matlab_mod(fsSub,0,fsdir);
    else
        [elecMatrix, elecLabels, elecRgb]=ut_mgrid2matlab_mod(fsSub,mgridFname,fsdir);
    end
end
            
nElec=length(elecLabels);
elecMatrix=round(elecMatrix);
xyz=zeros(size(elecMatrix));
xyz(:,1)=elecMatrix(:,2);
xyz(:,2)=elecMatrix(:,1);
xyz(:,3)=sVol(3)-elecMatrix(:,3);

% initialize marker for region membership (kept for potential debug)
elec_in_region = false(nElec,1);

% detect depth electrodes
depthElecs=zeros(nElec,1);
for a=1:nElec,
    if length(elecLabels{a})>=2 && strcmpi(elecLabels{a}(2),'D')
        depthElecs(a)=1;
    end
end

% If anatomical overlay requested, load segmentation and determine region IDs
region_segids = []; % default: empty -> draw everything
if universalYes(anatOverlay)
    if ischar(cfg.anatOverlay)
        segFname = cfg.anatOverlay;
    else
        segFname=fullfile(mriDir,'aparc+aseg.mgz');
    end
    if ~exist(segFname,'file')
        error('File %s not found.',segFname);
    end
    seg=MRIread(segFname);

    % Load segmentation color table
    if universalNo(colorLUT)
        pathstr = fileparts(which('mgrid2matlab'));
        inFile=fullfile(pathstr,'FreeSurferColorLUTnoFormat.txt');
        if ~exist(inFile,'file')
            error('Could not find file %s',inFile);
        end
    elseif exist(colorLUT,'file')
        inFile = colorLUT;
    else
        error('The defined color lookup table was not found');
    end
    fid=fopen(inFile,'r');
    tbl=textscan(fid,'%d%s%d%d%d%d');
    fclose(fid);
    
    % Map cfg.region (string or numeric) to region_segids
    if ~isempty(region_spec)
        if ischar(region_spec) || isStringScalar(region_spec)
            rs = lower(char(region_spec));
            switch rs
                case {'hippocampus','hip','hippocampi'} % INSERT MORE CASES AS YOU GO! 
                    region_segids = [17 53];
                otherwise
                    % try to match LUT names (case-insensitive)
                    lutNames = lower(tbl{2});
                    idx = find(contains(lutNames, rs));
                    if ~isempty(idx)
                        region_segids = tbl{1}(idx)';
                    else
                        region_segids = []; % unknown region string
                    end
            end
        elseif isnumeric(region_spec)
            region_segids = region_spec(:)';
        end
    end
    
    % compute elec_in_region boolean vector (for debugging / optional use)
    if ~isempty(region_segids)
        for ii = 1:nElec
            vx = elecMatrix(ii,1); vy = elecMatrix(ii,2); vz = elecMatrix(ii,3);
            if vx>=1 && vx<=size(seg.vol,1) && vy>=1 && vy<=size(seg.vol,2) && vz>=1 && vz<=size(seg.vol,3)
                elec_in_region(ii) = any(region_segids == seg.vol(vx,vy,vz));
            end
        end
    end
end

% Debug: show how many electrodes in the region
if ~isempty(region_spec) && ~isempty(region_segids)
    fprintf('Requested region: %s  -> segIDs: %s. Electrodes inside region: %d/%d\n', ...
        mat2str(region_spec), mat2str(region_segids), sum(elec_in_region), nElec);
elseif ~isempty(region_spec)
    fprintf('Requested region: %s  -> no matching segIDs found in LUT. Drawing all overlays.\n', mat2str(region_spec));
end

% Plotting loop (minimal changes: overlay drawing restricted to region_segids)
figId=0;
% --- minimalne, headless figure handling (TYLKO Visible = 'off') ---
figHandle = [];  % używamy uchwytu obiektowego (bez gcf)
for elecId=1:nElec
    if depthElecs(elecId)
        % utwórz lub użyj istniejącej figury (w zależności od figOverwrite)
        if isempty(figHandle) || ~isvalid(figHandle) || ~universalYes(figOverwrite)
            figHandle = figure('Visible','off');   % TYLKO to jedno ustawienie
        else
            figure(figHandle);                     % aktywuj istniejącą
            set(figHandle,'Visible','off');        % upewnij się, że jest niewidoczna
        end

        % czyść i ustaw rozmiar/pozycję przy użyciu uchwytu - BEZ gcf
        clf(figHandle);
        set(figHandle,'Position',[78 551 960 346],'PaperPositionMode','auto');

        hm=zeros(1,3);
        colormap gray;
        wdth=.35;
        wDelt=.33;
        xStart=-.005;
        yStart=.03;
        ht=.9;
        axes('Parent', figHandle,'position',[xStart yStart wdth ht]);
        imagesc(squeeze(mri.vol(:,xyz(elecId,2),:)),[mn mx]);
        axis square;
        set(gca,'xdir','reverse');
        hold on;
        
        if universalYes(anatOverlay)
            % Plot segmentation: ONLY draw patches for voxels in region_segids (if region specified)
            for a=1:sVol(1),
                for b=1:sVol(3),
                    curVal = seg.vol(a,xyz(elecId,2),b);
                    if curVal
                        if isempty(region_segids) || any(curVal == region_segids)
                            segId=find(tbl{1}==curVal);
                            tempRgb=double([tbl{3}(segId) tbl{4}(segId) tbl{5}(segId)])/255;
                            hM=patch([-.5 .5 .5 -.5]+b,[-.5 -.5 .5 .5]+a,tempRgb);
                            set(hM,'LineStyle','none','FaceAlpha',0.3);
                        end
                    end
                end
            end
        end
        
        % Plot electrode (use original elecRgb unless you explicitly want region color for electrode markers)
        hm(1)=plot(xyz(elecId,3),xyz(elecId,1),'r.');
        set(hm(1),'color',elecRgb(elecId,:),'markersize',markerSize);
        %find image limits
        mxX=max(squeeze(mri.vol(:,xyz(elecId,2),:)),[],2);
        mxY=max(squeeze(mri.vol(:,xyz(elecId,2),:)),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(1)/2),find(mxY==0)));
        limYb=min(intersect((sVol(1)/2:sVol(1)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax,
            axis([tempMin tempMax tempMin tempMax]);
        end
        set(gca,'xtick',[],'ytick',[]);
        
        %subplot(132);
        axes('Parent', figHandle,'position',[xStart+wDelt yStart wdth ht]);
        imagesc(squeeze(mri.vol(xyz(elecId,1),:,:)),[mn mx]);
        axis square;
        hold on;
        
        if universalYes(anatOverlay)
            for a=1:sVol(2),
                for b=1:sVol(3),
                    curVal = seg.vol(xyz(elecId,1),a,b);
                    if curVal
                        if isempty(region_segids) || any(curVal == region_segids)
                            segId=find(tbl{1}==curVal);
                            tempRgb=double([tbl{3}(segId) tbl{4}(segId) tbl{5}(segId)])/255;
                            hM=patch([-.5 .5 .5 -.5]+b,[-.5 -.5 .5 .5]+a,tempRgb);
                            set(hM,'LineStyle','none','FaceAlpha',0.3);
                        end
                    end
                end
            end
        end
        
        hm(2)=plot(xyz(elecId,3),xyz(elecId,2),'r.');
        set(hm(2),'color',elecRgb(elecId,:),'markersize',markerSize);
        %find image limits
        mxX=max(squeeze(mri.vol(xyz(elecId,1),:,:)),[],2);
        mxY=max(squeeze(mri.vol(xyz(elecId,1),:,:)),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(2)/2),find(mxY==0)));
        limYb=min(intersect((sVol(2)/2:sVol(2)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax,
            axis([tempMin tempMax tempMin tempMax]);
        end
        set(gca,'xtick',[],'ytick',[],'xdir','reverse');
        
        
        %subplot(133);
        axes('Parent', figHandle,'position',[xStart+wDelt*2 yStart wdth ht]);
        imagesc(squeeze(mri.vol(:,:,xyz(elecId,3))),[mn mx]);
        axis square;
        hold on;
        
        if universalYes(anatOverlay)
            for a=1:sVol(1),
                for b=1:sVol(2),
                    curVal = seg.vol(a,b,xyz(elecId,3));
                    if curVal
                        if isempty(region_segids) || any(curVal == region_segids)
                            segId=find(tbl{1}==curVal);
                            tempRgb=double([tbl{3}(segId) tbl{4}(segId) tbl{5}(segId)])/255;
                            hM=patch([-.5 .5 .5 -.5]+b,[-.5 -.5 .5 .5]+a,tempRgb);
                            set(hM,'LineStyle','none','FaceAlpha',0.3);
                        end
                    end
                end
            end
        end
        
        hm(3)=plot(xyz(elecId,2),xyz(elecId,1),'r.');
        set(hm(3),'color',elecRgb(elecId,:),'markersize',markerSize);
        %find image limits
        mxX=max(squeeze(mri.vol(:,:,xyz(elecId,3))),[],2);
        mxY=max(squeeze(mri.vol(:,:,xyz(elecId,3))),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(2)/2),find(mxY==0)));
        limYb=min(intersect((sVol(2)/2:sVol(2)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax,
            axis([tempMin tempMax tempMin tempMax]);
        end
        set(gca,'xtick',[],'ytick',[]);
        
        % Get anatomical label if aparc-file exists
        if exist(fullfile(mriDir,'aparc+aseg.mgz'),'file')
            anatLabel=ut_vox2Seg_mod(xyz(elecId,:),fsSub,fsdir);
        else
            anatLabel = 'NA';
        end
        
        % Remove first 3 characters that indicate hemisphere and electrode
        % type
        formattedLabel=elecLabels{elecId}(4:end);
        formattedLabel=rmChar(formattedLabel,'_'); % remove underscore between electrode stem and #
        
        if universalYes(fullTitle)
            ht=textsc2014([formattedLabel '; mgrid coords(' num2str(elecMatrix(elecId,:)-1) '); fsurf coords(' num2str(xyz(elecId,:)) '); ' anatLabel], ...
                'title');
            set(ht,'fontsize',14,'fontweight','bold');
        else
            ht=textsc2014([formattedLabel '; Anatomical Location: ' anatLabel], ...
                'title');
            set(ht,'fontsize',16,'fontweight','bold');
        end
        set(ht,'position',[.5 .97 0]);
        
        if ~universalNo(printFigs)
            if ischar(printFigs)
                outPath=printFigs; % user specified directory
            else
                erPath=fullfile(fsdir,fsSub,'elec_recon');
                outPath=fullfile(erPath,'PICS');
            end
            if ~exist(outPath,'dir')
                dirSuccess=mkdir(outPath);
                if ~dirSuccess,
                    error('Could not create directory %s',dirSuccess);
                end
            end
            
            drawnow;
            figFname=fullfile(outPath,sprintf('%s_%sSlices',fsSub,elecLabels{elecId}));
            figFname = regexprep(figFname, '[<>"/|?*]', '');
            fprintf('Exporting figure to %s\n',figFname);
            print(figHandle,figFname,'-djpeg');
            %pause(1);
        end
        
        if universalYes(pauseOn)
            fprintf('Paused. Press any key for next electrode.\n');
            pause;
        end
    end
end
fprintf('Done showing all electrodes.\n');
end