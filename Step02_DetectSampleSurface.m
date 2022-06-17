ROIY=[1 1];
ROIX=[1 1];

Thres=-5*10^-4;

ds_factor=1;
XROI=[];
YROI=[];

filelist=IDArdir([alignedTIFpath '*.tif']);
NFiles=numel(filelist);

ref_file=filelist(1).name;
fixed_im=imread(ref_file);
if size(fixed_im,3)>1
    fixed_im=rgb2gray(fixed_im);
end
if ds_factor<1
    fixed_im=imresize(fixed_im,ds_factor);
end
fixed_im_orig=fixed_im;
if ~isempty(XROI)
    fixed_im=fixed_im(:,XROI,:);
end
if ~isempty(YROI)
    fixed_im=fixed_im(YROI,:,:);
end

options = optimoptions('lsqcurvefit');
winsize=128;

B=zeros(5,NFiles);
if doRadialProjectionCorrection
    whichFiles=41:60;
    NwhichFiles=numel(whichFiles);
    im=cell(1,NwhichFiles);
    for ifile=1:NwhichFiles
        filename=[alignedTIFpath num2str(whichFiles(ifile),'%04.0f') '.tif'];
        im{ifile}=imread(filename);
    end
    im=cat(3,im{:});
    im=median(im,3);
    y=(1:size(im,1))';
    x=(1:size(im,2))';    
    [X,Y] = meshgrid(x,y);
    XY=cat(3,X,Y);
    % Create Objective Function: 
    surfit = @(B,XY)  B(1)*(XY(:,:,1).^2) + B(2)*(XY(:,:,2).^2) + ...
        B(3)*(XY(:,:,1)) + B(4)*(XY(:,:,2)) + B(5); 
    % Do Regression
    B = lsqcurvefit(surfit,randn(1,5),XY(1:4:end,1:4:end,:),double(im(1:4:end,1:4:end,:)),[],[],options);
    % Calculate Fitted Surface:
    Z = surfit(B,XY);
    Z = Z/max(Z(:));
end

I1=imread(filelist(1).name);
for ifile=1:NFiles
    if dosave
        [filepath,filename]=fileparts(filelist(ifile).name);
        savepath=[boundarypath filesep ];
        if ~isfolder(savepath)
            mkdir(savepath);
        end
        savefile=[savepath filename '_boundarycoords.mat'];
        if exist(savefile,'file') && ~dooverwrite
            continue;
        end
    end
    I1=imread(filelist(ifile).name);
    nI1=double(I1)./Z;
    med_anI1=movmedian(nI1,32*2);
    med_anI1=med_anI1-mean(med_anI1,1);

    I2=convn(nI1,ones(1,winsize)/winsize,'same');
    I2=(I2-min(I2(:)))/(max(I2(:))-min(I2(:)));

    dI2dy=diff(convn(I2,ones(30,1)/30,'same'),1,1);
    dI2dy=movmedian(dI2dy,32);
    ddI2ddy=diff(dI2dy,1,1);

    whichPoints=(1:1:size(I2,2))';
    NPoints=numel(whichPoints);
    candidates=nan(1,NPoints);
    for ix=1:NPoints
        values=-1*med_anI1(YStart:YEnd,whichPoints(ix));
        [~,peakinds] = findpeaks(values);
        dvalues=diff(values(peakinds),1,1);
        if numel(peakinds)<2
            continue;
        elseif numel(peakinds)==2
            continue;
        else
            dvalues=(dvalues-mean(dvalues))/std(dvalues);

            [~,candidateinds]=sort(dvalues,'descend');
            ind=[];
            if dvalues(candidateinds(1))-dvalues(candidateinds(2))>2
                ind=candidateinds(1);
            end

            if isempty(ind)
                ind=find(dvalues(2:end)>2 & dvalues(1:end-1)<0.1);
                if ~isempty(ind)
                    ind=ind+1;
                end
            end
            
            if isempty(ind)
                ind=find(dvalues>2);
            end
            if isempty(ind)
                continue;
            end
        end
        
        if isempty(ind)
            continue;
        end
        ind=ind(1);

        height=(values(peakinds(ind+1))-values(peakinds(ind)));
        target=values(peakinds(ind))+0.125*height;

        ind2=find(abs(values(peakinds(ind):min(numel(values),peakinds(ind)+300))-target)<0.03*height);
        if isempty(ind2)
            [~,ind2]=min(abs(values(peakinds(ind):min(numel(values),peakinds(ind)+300))-target));
        end
        ind2=ind2(1);
        ind2=peakinds(ind)-1+ind2;
        candidates(1,ix)=ind2+YStart-1;
    end   

    yPoints=candidates(:);

    invalid= isnan(yPoints);
    if sum(~invalid)>2
        [curve, goodness, output] =fit(whichPoints(~invalid),yPoints(~invalid),'linearinterp');
        curve=movmedian(curve(1:size(nI1,2)),64);
        Coords=cat(2,(1:size(nI1,2))',curve(1:size(nI1,2)));
    else
        Coords=[];
    end

    if dosave
        my_save(savefile,'Coords',Coords);
    end
    if doplot
        figure('position',get(0,'screensize')); 
        colormap gray;
        subplot(221);
        cla; imagesc(nI1); axis image; 
        hold on;
        plot(Coords(:,1),Coords(:,2));
        set(gca,'xlim',[ROIX(1),size(nI1,1)-ROIX(2)],'ylim',[ROIY(1),size(nI1,2)-ROIY(2)]);
        title([num2str(ifile) '/' num2str(NFiles)]);
        subplot(222); imagesc(I2(ROIY(1):end-ROIY(2)+1,ROIX(1):end-ROIX(2)+1)); axis image;
        hold on;
        plot(Coords(:,1),Coords(:,2));
        set(gca,'xlim',[ROIX(1),size(nI1,1)-ROIX(2)],'ylim',[ROIY(1),size(nI1,2)-ROIY(2)]);
        subplot(223); imagesc(dI2dy(ROIY(1):end-ROIY(2)+1,ROIX(1):end-ROIX(2)+1)); axis image; colorbar;
        title('dI2dy');
        set(gca,'clim',[Thres,-Thres]);
        hold on;
        plot(Coords(:,1),Coords(:,2));
        set(gca,'xlim',[ROIX(1),size(nI1,1)-ROIX(2)],'ylim',[ROIY(1),size(nI1,2)-ROIY(2)]);
        subplot(224); imagesc(ddI2ddy(ROIY(1):end-ROIY(2)+1,ROIX(1):end-ROIX(2)+1)); axis image; colorbar;
        set(gca,'clim',mean(reshape(ddI2ddy(ROIY(1):end-ROIY(2)+1,ROIX(1):end-ROIX(2)+1),[],1))+[-1,1]*std(reshape(ddI2ddy(ROIY(1):end-ROIY(2)+1,ROIX(1):end-ROIX(2)+1),[],1)))
        title('ddI2ddy');
        drawnow;
    end
end
