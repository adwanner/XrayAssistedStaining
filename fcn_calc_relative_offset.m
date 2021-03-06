function [dY1,dX1,dY2,dX2,cmax,C]= fcn_calc_relative_offset(varargin)
%inputs: fcn_calc_relative_offset(image1,image2,CorrMethod,MaxOffset,CorrThres,CenterAsRef,NPeaks)
dY1=[];
dY2=[];
dX1=[];
dX2=[];
cmax=[];
C=[];
NargIn=numel(varargin);

CorrMethod='redxcorr2normimages';
MaxOffset=[];
CorrThres=[];
CenterAsRef=false;

if NargIn>=1
    if ~isnumeric(varargin{1})
        disp('Error: Input variable 1 has to be numeric.')
        return;
    end
end

if NargIn>=2
    if ~isnumeric(varargin{2})
        disp('Error: Input variable 2 has to be numeric.')
        return;
    end
end

if NargIn>=3
    if ischar(varargin{3})
        CorrMethod=varargin{3};
    else
        disp('Error: Input variable 3 has to be a string,')
        disp('corresponding to one of the available correlation')
        disp('methods:')
        disp('''redxcorr2normimages''')
        disp('''redxcorr2normgradimages''')
        return;
    end
end

if NargIn>=4
    if isnumeric(varargin{4})
        MaxOffset=varargin{4};
    else
        disp('Input variable 4 has to be numeric,')
        disp('corresponding to the maximal offset to consider.')
        disp('It is either a scalar or a 4 elements vector:')
        disp('maxOffset=[Top,Bottom,Left,Right];')
        return;
    end
end

if NargIn>=5
    if isnumeric(varargin{5}) && numel(varargin{5})==1
        CorrThres=varargin{5};
    else
        if ~isempty(varargin{5})
            disp('Input variable 5 has to be numeric scalar,')
            disp('corresponding to the correlation threshold to consider.')
            return;
        end
    end
end

if NargIn>=6
    if isnumeric(varargin{6})
        varargin{6}=logical(varargin{6});
    end
    if islogical(varargin{6}) && numel(varargin{6})==1
        CenterAsRef=varargin{6};
    else
        if ~isempty(varargin{6})
            disp('Input variable 6 has to be a logical scalar value,')
            disp('corresponding to wheter the center of the images should')
            disp('be used as reference for maximal offset.')
            return;
        end
    end
end

if NargIn>=7
    NPeaks=varargin{7};
else
    NPeaks=1;
end
if isempty(NPeaks)
    NPeaks=1;
end

%we want size(A)>=size(T)
if size(varargin{1},1)>size(varargin{2},1) || size(varargin{1},2)>size(varargin{2},2)
    is1greater2=1;
    A=varargin{1};
    T=varargin{2};
else
    is1greater2=0;
    A=varargin{2};    
    T=varargin{1};
end

if any([std(double(A(:))),std(double(T(:)))]==0)
    disp('Error: The values of input variable 1 and/ or 2 must not be')
    disp('all the same.')
    return;
end

clear varargin;
[THeight,TWidth,TNSlices]=size(T);
[AHeight,AWidth,ANSlices]=size(A);

if CenterAsRef
    CmaxHeight=2*floor(THeight/2)-floor((THeight-1)/2)+...
            floor((AHeight+1)/2)-1;
    CmaxWidth=2*floor(TWidth/2)-floor((TWidth-1)/2)+...
            floor((AWidth+1)/2)-1;
    switch numel(MaxOffset)
        case 0
            MaxOffset=[CmaxHeight,CmaxHeight,...
                CmaxWidth,CmaxWidth];
        case 1
            MaxOffset=[...
                min(MaxOffset,CmaxHeight),...
                min(MaxOffset,CmaxHeight),...
                min(MaxOffset,CmaxWidth),...
                min(MaxOffset,CmaxWidth)];
        case 2
            MaxOffset=[...
                min(MaxOffset(1),CmaxHeight),...
                min(MaxOffset(1),CmaxHeight),...
                min(MaxOffset(2),CmaxWidth),...
                min(MaxOffset(2),CmaxWidth)];
        case 4
            MaxOffset=[...
                min(MaxOffset(1),CmaxHeight),...
                min(MaxOffset(2),CmaxHeight),...
                min(MaxOffset(3),CmaxWidth),...
                min(MaxOffset(4),CmaxWidth)];
        otherwise
            disp('Error: The input for ''MaxOffset'' is not valid.')
            return;
    end
else
    switch numel(MaxOffset)
        case 0
            MaxOffset=[THeight,AHeight,...
                TWidth,AWidth];
        case 1
            MaxOffset=[...
                min(MaxOffset,THeight),...
                min(MaxOffset,AHeight),...
                min(MaxOffset,TWidth),...
                min(MaxOffset,AWidth)];
        case 2
            MaxOffset=[...
                min(MaxOffset(1),THeight),...
                min(MaxOffset(1),AHeight),...
                min(MaxOffset(2),TWidth),...
                min(MaxOffset(2),AWidth)];
        case 4
            MaxOffset=[...
                min(MaxOffset(1),THeight),...
                min(MaxOffset(2),AHeight),...
                min(MaxOffset(3),TWidth),...
                min(MaxOffset(4),AWidth)];
        otherwise
            disp('Error: The input for ''MaxOffset'' is not valid.')
            return;
    end
end
% fftsize(1) = FindClosestValidDimension(THeight+AHeight-1);
% fftsize(2) = FindClosestValidDimension(TWidth+AWidth-1);
fftsize(1)=2^nextpow2(THeight+AHeight-1);
fftsize(2)=2^nextpow2(TWidth +AWidth -1);

if CenterAsRef
    centerHeight=CmaxHeight+[-(MaxOffset(1)-1),(MaxOffset(2)-1)];
    centerWidth= CmaxWidth+[-(MaxOffset(3)-1),(MaxOffset(4)-1)];
else
    centerHeight=THeight+[-(MaxOffset(1)-1),(MaxOffset(2)-1)];
    centerWidth= TWidth +[-(MaxOffset(3)-1),(MaxOffset(4)-1)];
end

switch CorrMethod
    case 'redxcorr2normimages' %reduced xcorr2 with normalized images
    case 'redxcorr2normgradimages' %reduced xcorr2 with normalized gradient of images
        A=gradient(double(A));
        T=gradient(double(T));
%         [Ax,Ay]=gradient(double(A));
%         [Tx,Ty]=gradient(double(T));
%         A=Ax+Ay;
%         T=Tx+Ty;
%         clear Ax Ay Tx Ty
    otherwise
        CorrMethod='redxcorr2normimages'; %#ok<NASGU>
end
dX=nan(TNSlices,ANSlices);
dY=nan(TNSlices,ANSlices);
cmax=nan(TNSlices,ANSlices);
if nargout>=6
    C=zeros(centerHeight(2)-centerHeight(1)+1,...
        centerWidth(2)-centerWidth(1)+1,TNSlices,ANSlices);
end
for iaslice=1:ANSlices
    nA=rowcol_normalization(A(:,:,iaslice));
    if std(nA(:))==0
        continue;
    end
    fA=fft(nA,fftsize(1),1);
%         clear nA
    fA=permute(fA,[2,1]);
    fA=fft(fA,fftsize(2),1);
    fA=permute(fA,[2,1]);
    for itslice=1:TNSlices
        calcOK=false;
        while ~calcOK
            nT=rowcol_normalization(T(:,:,itslice));
            if std(nT(:))==0
                dY(itslice,iaslice)=NaN;
                dX(itslice,iaslice)=NaN;
                cmax(itslice,iaslice)=NaN;
                break;
            end
            nT=nT(end:-1:1,end:-1:1);
            fT=fft(nT,fftsize(1),1);
%             clear nT
            fT=permute(fT,[2,1]);
            fT=fft(fT,fftsize(2),1);
            fT=permute(fT,[2,1]);
%             inplaceprod(fT,fA);
            fT=fT.*fA;
            fC=fT;
%             clear fA fT
            if (MaxOffset(1)+MaxOffset(2))*fftsize(2)<=...
                    (MaxOffset(3)+MaxOffset(4))*fftsize(1)
                fC=ifft(fC,fftsize(1),1);
                fC=permute(fC,[2,1]);
                fC=fC(:,centerHeight(1):centerHeight(2));
                fC=ifft(fC,fftsize(2),1);
                fC=real(fC);
                fC=permute(fC,[2,1]);
                tempC=fC(:,centerWidth(1):centerWidth(2));
            else
                fC=permute(fC,[2,1]);
                fC=ifft(fC,fftsize(2),1);
                fC=permute(fC,[2,1]);
                fC=fC(:,centerWidth(1):centerWidth(2));
                fC=ifft(fC,fftsize(1),1);
                fC=real(fC);
                fC=permute(fC,[2,1]);
                fC=fC(:,centerHeight(1):centerHeight(2));
                tempC=permute(fC,[2,1]);            
            end
            
%             p=FastPeakFind(tempC);
%             p=reshape(p,2,[]);
%             ind=reshape(sub2ind(size(tempC),p(2,:),p(1,:)),[],1);
%             [~,order]=sort(tempC(ind),'descend');
%             ind=ind(order);
%             p=p(:,order);
%             imagesc(tempC); axis image; hold all;
%             plot(p(1,NPeaks),p(2,NPeaks),'r+')
%             uiwait();
            if NPeaks==1
                [ph ,ind] = max(tempC(:)); %#ok<ASGLU>
            else
                [ph,ind]=sort(tempC(:),'descend'); %#ok<ASGLU>
            end
%             [iy2 ,ix2] = ind2sub(size(tempC),ind2(1));
%             plot(ix2,iy2,'k+')            
            
            [iy ,ix] = ind2sub(size(tempC),ind(NPeaks));

            weight=tempC(ind(NPeaks));
            weight=weight/sum(weight);
            if CenterAsRef
                dY(itslice,iaslice)=iy'*weight-MaxOffset(1);
                dX(itslice,iaslice)=ix'*weight-MaxOffset(3);
            else
                dY(itslice,iaslice)=iy'*weight-MaxOffset(1);
                dX(itslice,iaslice)=ix'*weight-MaxOffset(3);
            end
            
            
%             overlap_area=(min(THeight,AHeight)-dY(itslice,iaslice))*...
%                 (min(TWidth,AWidth)-dX(itslice,iaslice));
%             overlap_area=min(THeight,AHeight)*min(TWidth,AWidth);
%             cmax(itslice,iaslice)=cmax(itslice,iaslice)/overlap_area;
            [iOverlapX,iOverlapY,~,jOverlapX,jOverlapY]=...
            fcn_get_overlap(...
                0,0,0,TWidth,THeight,1,...
                -dX(itslice,iaslice),-dY(itslice,iaslice),0,AWidth,AHeight,1);
            iOverlapX=round(iOverlapX);
            iOverlapY=round(iOverlapY);
            jOverlapX=round(jOverlapX);
            jOverlapY=round(jOverlapY);
            overlapT=double(T(iOverlapY(1):iOverlapY(2),iOverlapX(1):iOverlapX(2),itslice));
            overlapA=double(A(jOverlapY(1):jOverlapY(2),jOverlapX(1):jOverlapX(2),iaslice));

            NPx=numel(overlapT);
            sumOp=1/NPx*ones(1,NPx);
            
            overlapT=overlapT(:);
            tempMean=sumOp*overlapT;
            overlapT=overlapT-tempMean;
            tempStd=overlapT'*overlapT;
            tempStd=sqrt(tempStd/(NPx-1));
            overlapT=overlapT/tempStd;

            overlapA=overlapA(:);
            tempMean=sumOp*overlapA;
            overlapA=overlapA-tempMean;
            tempStd=overlapA'*overlapA;
            tempStd=sqrt(tempStd/(NPx-1));
            overlapA=overlapA/tempStd;

            cmax(itslice,iaslice)=overlapT'*overlapA/(NPx-1);

            if nargout>=6
                C(:,:,itslice,iaslice)=tempC;%#ok<AGROW> %/overlap_area;
            end
 
%             tempC = tempC - mean(tempC(:));
%             tempC = tempC/sqrt(tempC(:)'*tempC(:)/length(tempC(:)));        
%             [cmax(itslice,iaslice) ilin] = max(tempC(:));
%             C(:,:,itslice,iaslice)=tempC;
            
%             
%             imagesc(tempC); colorbar;
%             impoint(gca,ix,iy);
%             drawnow;
            if isempty(CorrThres)
                calcOK=true;
            else
                if cmax(itslice,iaslice)<CorrThres && ~all(MaxOffset==[THeight,AHeight,...
                        TWidth,AWidth])
                    disp(['cmax= ' num2str(cmax(itslice,iaslice)) '. CorrThres= ' num2str(CorrThres) '.'])
                    disp('Re-run calculations of offsets with new MaxOffset boundaries:')
                    MaxOffset=[THeight,AHeight,...
                        TWidth,AWidth];            
                    disp(num2str(MaxOffset));
                    calcOK=false;
                else
                    calcOK=true;
                end
            end
        end
    end        
end
if is1greater2
    dY1=dY;
    dX1=dX;
    dY2=-dY;
    dX2=-dX;
else
    dY1=-dY;
    dX1=-dX;
    dY2=dY;
    dX2=dX;
end
% close all;
% figure; imagesc(C); colorbar;
% impoint(gca,ix,iy);
% drawnow;
end

function A=rowcol_normalization(A)
    [AHeight,AWidth]=size(A);
    A=double(A);
    %% first normalize rows (in X direction) 
    if AWidth>1
        Mean=A*ones(AWidth,1)/AWidth; %faster than sum(...) or mean
        A=A-Mean*ones(1,AWidth);
        Norm=sqrt((A.*A)*ones(AWidth,1)/AWidth);
        Norm(Norm==0)=1;
        A=A./(Norm*ones(1,AWidth));
%         Norm=1./Norm;
%         inplaceprod(A,Norm*ones(1,AWidth));
    end
    if AHeight>1
        %% 2nd normalize columns (in Y direction) 
        A=A';
        Mean=A*ones(AHeight,1)/AHeight; %faster than sum(...) or mean
        A=A-Mean*ones(1,AHeight);
        Norm=sqrt((A.*A)*ones(AHeight,1)/AHeight);
        Norm(Norm==0)=1;
        A=A./(Norm*ones(1,AHeight));
%         Norm=1./Norm;
%         inplaceprod(A,Norm*ones(1,AHeight));
        A=A';
    end
end