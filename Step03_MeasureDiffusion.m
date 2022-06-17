
xbinsize=20; %um
minDiffDistance=200;%um
maxDiffDistance=1400;%um

DefaultRange=500:1400;

filelist=IDArdir([boundarypath filesep '*_boundarycoords.mat']);
NFiles=numel(filelist);
clear BoundaryCoords
for ifile=1:NFiles
    load(filelist(ifile).name,'Coords');
    if isempty(Coords)
        continue;
    end
    if ~exist('BoundaryCoords','var')
        BoundaryCoords=nan(size(Coords,1),NFiles);
    end
    BoundaryCoords(:,ifile)=Coords(:,2);
end
if ~exist('BoundaryCoords','var')
    return;
end

clear TargetRes CoordRange
load(metadatafile,'TargetRes');
load(metadatafile,'RecordingDate');
load(metadatafile,'CoordRange');
acqtime=(RecordingDate-min(RecordingDate))*24*60;

if ~exist('CoordRange','var')
    CoordRange=DefaultRange;
end


options = optimoptions('lsqcurvefit');
expansionmodel = @(B,t)  B(1)-B(2)*exp(B(3)*t)+B(4)*t; 

x0=[200,200,-0.0001,0];

if ~doExpansionCorrection
    densityfile=strrep(densityfile,'.mat','_noexpansioncorrection.mat');
end
targetpath=fileparts(densityfile);
if ~isfolder(targetpath)
    mkdir(targetpath);
end

if calc_diffusion && (~(exist(densityfile,'file') && ~dooverwrite))
    filename=[alignedTIFpath filesep num2str(NFiles,'%04.0f') '.tif'];

    avBoundary=nanmedian(BoundaryCoords(CoordRange,1:min(1+9,size(BoundaryCoords,2))),2)';

    im=imread(filename);

    NSamples=1;
    position=cell(NSamples,1);
    for isample=1:NSamples
        figure;
        imagesc(im); colormap(gca,'gray'); axis image; hold all;
        for ipos=1:NSamples
            if isempty(position{ipos})
                continue;
            end
            x1=position{ipos}(1,1);
            y1=position{ipos}(1,2);
            x2=position{ipos}(2,1);
            y2=position{ipos}(2,2);
            plot([x1,x2],[y1,y2],'k');
        end
        xdata=Coords(CoordRange,1);
        ydata=avBoundary;
        plot(xdata,ydata,'r');            
        xdata=Coords(CoordRange,1);
        ydata=BoundaryCoords(CoordRange,NFiles);
        plot(xdata,ydata,'b');
        h = imline(gca); %#ok<IMLINE>

        position{isample} = wait(h);
        close(gcf);
    end
    
    linewidth=50; %px
    linestart=-100;
    lineend=800;
    linelength=lineend-linestart+1; %px

    if doRadialProjectionCorrection
        whichFiles=2:5;
        NwhichFiles=numel(whichFiles);
        im=cell(1,NwhichFiles);
        for ifile=1:NwhichFiles
            filename=[alignedTIFpath filesep num2str(whichFiles(ifile),'%04.0f') '.tif'];
            im{ifile}=imread(filename);
        end
        im=cat(3,im{:});
        im=median(im,3);
        y=(1:size(im,1))';
        x=(1:size(im,2))';    
        [X,Y] = meshgrid(x,y);
        XY=cat(3,X,Y);
        % Create Objective Function: 
        surfit = @(B,XY) ...
            B(5)*10^-6*(XY(:,:,1).^2) + B(4)*10^-6*(XY(:,:,2).^2) + ...
            B(3)*10^-3*(XY(:,:,1))    + B(2)*10^-3*(XY(:,:,2))    + B(1); 
        % Do Regression
        B = lsqcurvefit(surfit,rand(1,5),XY(1:4:end,1:4:end,:),log(double(im(1:4:end,1:4:end,:))),[],[],options);
        % Calculate Fitted Surface:
        Z = exp(surfit(B,XY));
        Z = Z/max(Z(:));
    end

    mean_density=zeros(linelength,NSamples,NFiles);
    std_density=zeros(linelength,NSamples,NFiles);
    parfor_progress(NFiles);
    firstROIpos=cell(NSamples,1);
    lastROIpos=cell(NSamples,1);
    for ifile=1:NFiles
        parfor_progress();
        if doExpansionCorrection
            avBoundary=nanmedian(BoundaryCoords(:,ifile:min(ifile+9,size(BoundaryCoords,2))),2); %#ok<*NANMEDIAN>
        else
            avBoundary=nanmedian(BoundaryCoords(:,1:min(1+9,size(BoundaryCoords,2))),2); %#ok<UNRCH>
        end


        filename=[alignedTIFpath filesep num2str(ifile,'%04.0f') '.tif'];
        im=imread(filename);

        if doRadialProjectionCorrection
            im=double(im)./Z;      
        end            

        switch ifile
            case 1
                firstIm=im;
                firstTime=acqtime(ifile);
            case NFiles
                lastIm=im;
                lastTime=acqtime(ifile);
        end

        xSurface=Coords(CoordRange,1);
        ySurface=avBoundary(CoordRange);

        for ipos=1:NSamples
            x1=position{ipos}(1,1);
            y1=position{ipos}(1,2);
            x2=position{ipos}(2,1);
            y2=position{ipos}(2,2);
            n=[y2-y1,-(x2-x1)];
            n=n/norm(n);

            r=[x2-x1,y2-y1];
            r=r/norm(r);

            dist2surface2=sum([xSurface-x1,ySurface-y1].^2,2);
            [~,ind]=min(dist2surface2);            

            xPx=nan(linelength,linewidth);
            yPx=nan(linelength,linewidth);
            for iwidth=1:linewidth
                xPx(:,iwidth)=xSurface(ind)+n(1)*(-linewidth/2+iwidth)+r(1)*(linestart:lineend);
                yPx(:,iwidth)=ySurface(ind)+n(2)*(-linewidth/2+iwidth)+r(2)*(linestart:lineend);
            end
            CTIntensity=interp2(double(im),xPx,yPx);

            switch ifile
                case 1
                    firstROIpos{ipos}=[...
                        xSurface(ind)+r(1)*[linestart,lineend];...
                        ySurface(ind)+r(2)*[linestart,lineend]];
                case NFiles
                    lastROIpos{ipos}=[...
                        xSurface(ind)+r(1)*[linestart,lineend];...
                        ySurface(ind)+r(2)*[linestart,lineend]];
            end

            mean_density(:,ipos,ifile)=mean(CTIntensity,2);
            std_density(:,ipos,ifile)=std(CTIntensity,[],2);
        end
        drawnow;
    end

    [xdata,tdata]=ndgrid((linestart:lineend)*TargetRes,acqtime);

    if dosave
        save(densityfile,'xdata','tdata','mean_density','std_density','position','B',...
            'firstIm','firstTime','lastIm','lastTime','firstROIpos','lastROIpos');
    end
end
