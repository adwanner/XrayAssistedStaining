maxOffset=128;
Neighbors2Consider=3;
RandPairs2Consider=0;
dovis=0;

if doextract_xrm
    [~,tempbasefilename]=fileparts(in_filename);
    tempbasepath=[TIFpath filesep];
    if ~isfolder(tempbasepath)
        mkdir(tempbasepath)
    end
    scriptpath= mfilename('fullpath');
    scriptpath=fileparts(scriptpath);
    systemCommand = ['cd "' scriptpath '"; python2.7 read_txrm_metadata.py' ...
        ' "' in_filename ...
        '" "' tempbasepath ...
        '" "' tempbasefilename '"'];
    [~, result] = system(systemCommand);
    if contains(result,{'error','Error'})
        disp(result);
        if contains(result, 'No module named olefile')
            disp('Try running "pip2.7 install --user olefile"');
        end
        if contains(result, 'No module named numpy')
            disp('Try running "pip2.7 install --user numpy"');
        end
        if contains(result, 'No module named libtiff')
            disp('Try running "pip2.7 install --user libtiff (maybe you need to run "sudo apt-get install python-dev")"');
        end
        return;
    end
end

filelist=IDArdir([TIFpath filesep '*.tif']);
NFiles=numel(filelist);

if read_metadata
    ImRes=zeros(NFiles,1);
    ImSize=zeros(2,NFiles);
    RecordingDate=zeros(NFiles,1);
    for ifile=1:NFiles
        filename=strrep(filelist(ifile).name,'.tif','.txt');
        fid=fopen(filename,'r');
        str=fgetl(fid);
        while ischar(str)
            if contains(str,'pixelsize')
                pos=strfind(str,' = ')+numel(' = ');
                value=str2double(str(pos:end));
                ImRes(ifile)=value;
            end
            if contains(str,'Date')
                pos=strfind(str,' = ')+numel(' = ');
                value=datenum(str(pos:end),'mm/dd/yy HH:MM:SS');
                RecordingDate(ifile)=value;
            end
            if contains(str,'ImageHeight')
                pos=strfind(str,' = ')+numel(' = ');
                value=str2double(str(pos:end));
                ImSize(1,ifile)=value;
            end
            if contains(str,'ImageWidth')
                pos=strfind(str,' = ')+numel(' = ');
                value=str2double(str(pos:end));
                ImSize(2,ifile)=value;
            end
            str=fgetl(fid);
        end
        fclose(fid);
    end
    [~,order]=sort(RecordingDate);
    filelist=filelist(order);
    ImRes=ImRes(order);
    ImSize=ImSize(:,order);

    RecordingDate=RecordingDate(order);
    TargetRes=max(ImRes);

    my_save(metadatafile,'ImRes',ImRes,'filelist',filelist,...
        'ImSize',ImSize,'RecordingDate',RecordingDate,...
        'TargetRes',TargetRes);
else
    load(metadatafile); %#ok<UNRCH>
end    

if calc_offsets
    NPairs=NFiles*(Neighbors2Consider+RandPairs2Consider);
    Pairs=zeros(2,NPairs);
    ipair=0;
    for ifile=1:(NFiles-1)
        for jfile=ifile+1:min(ifile+Neighbors2Consider,NFiles)
            ipair=ipair+1;
            Pairs(:,ipair)=[ifile,jfile];
        end
        if any(RandPairs2Consider)
            whichFiles=setdiff(1:NFiles,ifile:min(ifile+Neighbors2Consider,NFiles));
            whichFiles=whichFiles(randi(numel(whichFiles),RandPairs2Consider,1));
            for jfile=whichFiles
                ipair=ipair+1;
                Pairs(:,ipair)=[ifile,jfile];
            end
        end
    end
    NPairs=ipair;
    Pairs(:,NPairs+1:end)=[];

    NPairs=size(Pairs,2);
    dY=nan(NPairs,1);
    dX=nan(NPairs,1);
    cmax=nan(NPairs,1);
    parfor ipair=1:NPairs
        ifile=Pairs(1,ipair);
        jfile=Pairs(2,ipair); %#ok<PFBNS>

        iImage=imread(filelist(ifile).name); %#ok<PFBNS>
        iImage=imresize(iImage,ImRes(ifile)/TargetRes); %#ok<PFBNS>
        riImage=iImage;

        jImage=imread(filelist(jfile).name);
        jImage=imresize(jImage,ImRes(jfile)/TargetRes);
        rjImage=jImage;

        Peaks2Consider=2; NReps=10;
        cdY=nan(NReps,1);
        cdX=nan(NReps,1);
        tempCmax=nan(NReps,1);
        for irep=1:NReps
            [~,~,tempdY,tempdX,tempCmax(irep),C]= fcn_calc_relative_offset(...
                riImage,rjImage,'redxcorr2normimages',maxOffset,[],false,Peaks2Consider);
            disp([num2str(tempdY) ', ' num2str(tempdX) ', ' num2str(tempCmax(irep))]);
            cdX(irep)=tempdX;
            cdY(irep)=tempdY;
            if cdX(irep)==0 && cdY(irep)==0
                break;
            end
            rjImage=imtranslate(jImage,-[sum(cdX(1:irep)) sum(cdY(1:irep))],'FillValues',0);            
        end
        tempdY=nansum(cdY); %#ok<*NANSUM>
        tempdX=nansum(cdX);
        cmax(ipair)=nanmean(tempCmax); %#ok<*NANMEAN>

        dY(ipair)=tempdY;
        dX(ipair)=tempdX;
        if dovis
            iy=tempdY+maxOffset; %#ok<UNRCH>
            ix=tempdX+maxOffset;
            figure; 
            jtImage=imtranslate(jImage,-[tempdX tempdY],'FillValues',0);            
            ijImage = imfuse(iImage,jtImage,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
            imshow(ijImage);
            drawnow;

            Stack=cat(3,iImage,jtImage);
            fig=figure('position',get(0,'screensize')); 
            islice=1;
            ax=axes('parent',fig);
            handles=imagesc(Stack(:,:,islice),'parent',ax); axis image; 
            set(handles(1),'UserData',islice);
            set(gca,'clim',[0,6000]);
            colormap(gray(256));
            set(fig,'WindowScrollWheelFcn',@(src,evnt)next_slice_scrollwheel_cbk(...
                {Stack},handles,src,evnt));                
        end
    end
    parfor_progress(0);
    my_save(displacementfile,'ImRes',ImRes,'filelist',filelist,...
        'ImSize',ImSize,'RecordingDate',RecordingDate,...
        'TargetRes',TargetRes,...
        'dY',dY,'dX',dX,'Pairs',Pairs,'cmax',cmax);
else
    load(displacementfile); %#ok<UNRCH>
end

if writeAlignedImage
    NPairs=size(Pairs,2);
    validPairs=true(1,NPairs);
    NPairs=sum(validPairs);

    W=sparse([1:NPairs,1:NPairs],...
        [Pairs(1,validPairs),Pairs(2,validPairs)],...
        [ones(1,NPairs),-1*ones(1,NPairs)],NPairs,NFiles,2*NPairs);
    Offsets=[dY(validPairs),dX(validPairs)];
    Coords=W\Offsets;
    bp=1:100:size(Coords,1);
    Coords(:,1) = detrend(Coords(:,1),'linear',bp);
    Coords(:,2) = detrend(Coords(:,2),'linear',bp);    

    bp=51:100:size(Coords,1);
    Coords(:,1) = detrend(Coords(:,1),'linear',bp);
    Coords(:,2) = detrend(Coords(:,2),'linear',bp);    

    Coords(:,1)=Coords(:,1)-min(Coords(:,1));
    Coords(:,2)=Coords(:,2)-min(Coords(:,2));

    targetSize=floor(min(ImSize'-Coords));
    if ~isfolder(alignedTIFpath)
        mkdir(alignedTIFpath);
    end
    for ifile=1:NFiles
        filename=num2str(ifile,'%04.0f');
        alignedfile=[alignedTIFpath filename '.tif'];
        if exist(alignedfile,'file')
            continue;
        end
        iImage=imread(filelist(ifile).name);
        iImage=imtranslate(iImage,Coords(ifile,[2,1]),'FillValues',0);
        iImage=iImage((end-targetSize(1)+1):end,(end-targetSize(2)+1):end);
        imwrite(iImage,alignedfile);
        disp(['wrote ' num2str(ifile)]);
    end
end