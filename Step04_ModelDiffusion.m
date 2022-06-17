global xdata tdata diffusion tmax xmax optimizer xfittingRange tfittingRange
global exterior doplotprogress doweighted domasking
global dogrowth IntensityNorm


doweighted=0;


DefaultRange=500:1400;

xfittingRange=[50 1000]/1000;

depthrange=[-280,1500]/1000;
baselinerange=[-150,-50]/1000;

validation_xrange=[0,1000]/1000;
validation_trange=[0,1200];

tfittingRange=[0,inf];
tmax=1;
xmax=1;


if exist(modelfile,'file') && ~dooverwrite && dosavedata
    return;
end

filelist=IDArdir([boundarypath '*_boundarycoords.mat']);
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
clear TargetRes CoordRange
load(metadatafile,'TargetRes');
load(metadatafile,'RecordingDate');
load(metadatafile,'CoordRange');
acqtimes=(RecordingDate-min(RecordingDate))*24*60;

if ~exist('CoordRange','var')
    CoordRange=DefaultRange;
end

im_start=load(densityfile,'firstIm');

ydata=BoundaryCoords(CoordRange,:);
ydata=size(im_start,1)-ydata;
ydata=ydata./nanmean(ydata(:,1:min(10,size(ydata,2))),2); %#ok<*NANMEAN>
ydata=((nanmean(ydata-nanmean(ydata(:,1:min(10,size(ydata,2))),2))));
expansion=100*ydata;

options = optimoptions('lsqcurvefit');
expansionmodel = @(B,t)  B(1)-B(2)*exp(-B(3)*t); 
x0=[5.8697    5.6924    0.0017];    
global expansioncoeffs
expansioncoeffs = lsqcurvefit(expansionmodel,x0,...
    acqtimes,...
    expansion(:));            

[raw_diffusion,raw_xdata,raw_tdata]=microCT_extractDensity(densityfile,depthrange*1000);
raw_xdata=raw_xdata/1000;
raw_diffusion=raw_diffusion(:,:,1);

baseline=raw_xdata(:,1)'<=baselinerange(2) & raw_xdata(:,1)'>=baselinerange(1);
tempbaseline=mean(raw_diffusion(baseline,:),1);

if dobaselinecorr
    baselinecorrection=tempbaseline-median(tempbaseline);
    normdiffusion=raw_diffusion-baselinecorrection;
else   
    normdiffusion=raw_diffusion; %#ok<UNRCH>
end
%     
cdiffusion=normdiffusion;    

if donorm
    maxValue=6420;
    valid=raw_xdata>=-200/1000 & raw_xdata<=0 & ...
        raw_tdata>=0 & raw_tdata<=1200;
    baseline=mean(cdiffusion(valid));
    minValue=baseline-(maxValue-baseline)/0.9*0.1;    

    a=0.61749;
    IntensityRef=minValue-a*(maxValue-minValue);
    IntensityNorm=(maxValue-minValue)*(1+a);
    cdiffusion=(cdiffusion-IntensityRef)/IntensityNorm;
else
    cdiffusion=cdiffusion/10000; %#ok<UNRCH>
end

tdata=0:5:min(1200,max(raw_tdata(:)));
xdata=(baselinerange(1):0.005:depthrange(2));
tdata=ones(numel(xdata),1)*tdata;
xdata=xdata'*ones(1,size(tdata,2));

interp_fcn=scatteredInterpolant(raw_tdata(:),raw_xdata(:),cdiffusion(:));
interp_cdiffusion=interp_fcn(tdata,xdata);
interp_maxDiffusion=max(reshape(interp_cdiffusion,[],1));
interp_minDiffusion=min(reshape(interp_cdiffusion,[],1));

valid_files(ifile)=true;

diffusion=interp_cdiffusion;


baseline=xdata(:,1)'<=baselinerange(2) & xdata(:,1)'>=baselinerange(1);

maxDiffusion=max(diffusion(:));
minDiffusion=min(diffusion(:));

exterior=mean(reshape(diffusion(baseline,:),[],1));

xdata=xdata/xmax;
tdata=tdata/tmax;


if dosavedata && ~exist(modelfile,'file')
    save(modelfile,...
        'depthrange','baselinerange',...
        'xfittingRange','tfittingRange','tmax','xmax','diffusion','tdata','xdata',...
        'minDiffusion','maxDiffusion','exterior','expansioncoeffs',...
        'validation_xrange','validation_trange',...
        'IntensityRef','IntensityNorm');
end


if ~domodel
    return;
end


p0=[0.175680160795783   0.333314089809591   0.177002719038096   0.388441923147603   0.124182201469362   0.226634860359746];
NParameters=numel(p0);

optimizer='lsqnonlin';
global fig
fig=nan(1,2);
switch optimizer
    case 'fminunc'
        options = optimoptions('fminunc','Display','iter');
        p = fminunc(@fitfcn,p0,options);
    case 'lsqnonlin'
        lb=zeros(1,NParameters);
        ub=inf(1,NParameters);
        
        
        options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','Display','iter','FiniteDifferenceStepSize',1e-3);
        [p,resnorm,residual,exitflag,output] = lsqnonlin(@fitfcn,p0,lb,ub,options);
end

tic

tspan=[0:1:9,10:5:(60*30)]/tmax;
xspan=[0:1:9,10:5:3000]/1000/xmax;

[T,C,B,L]=diffreactmodel(p,xspan,tspan);

T=reshape(T,[numel(xspan),numel(tspan)]);
C=reshape(C,[numel(xspan),numel(tspan)]);
B=reshape(B,[numel(xspan),numel(tspan)]);
L=reshape(L,[numel(xspan),numel(tspan)]);

[xgrid,tgrid]=ndgrid(xspan,tspan);
interp_fcn=scatteredInterpolant(tgrid(:),xgrid(:),T(:));
interpT=interp_fcn(tdata,xdata);
interp_fcn=scatteredInterpolant(tgrid(:),xgrid(:),C(:));
interpC=interp_fcn(tdata,xdata);
interp_fcn=scatteredInterpolant(tgrid(:),xgrid(:),B(:));
interpB=interp_fcn(tdata,xdata);
interp_fcn=scatteredInterpolant(tgrid(:),xgrid(:),L(:));
interpL=interp_fcn(tdata,xdata);


valid=xmax*xdata(:)>=validation_xrange(1) & xmax*xdata(:)<=validation_xrange(2) & ...
    tmax*tdata(:)>=validation_trange(1) & tmax*tdata(:)<=validation_trange(2);
SSres=sum((diffusion(valid)-interpT(valid)).^2);
%degrees of freedom of residuals
k=numel(p);
n=numel(diffusion(valid));
df=n-k-1;
%Residual Standard Error
S=sqrt(SSres/df);
if dosavedata
    save(modelfile,...
        'p','xspan','tspan',...
        'T','C','B','L','tdata','xdata','S','-append');
end

whichX=([200,400,600,800,1000,1200]/1000)/xmax;
fig(1)=figure('position',get(0,'screensize')); 
for iplot=1:numel(whichX)
    subplot(2,3,iplot);
    [~,whichDataInd]=min(abs(xdata(:,1)-whichX(iplot)));
    plot(tdata(1,:)*tmax,interpT(whichDataInd,:)); hold all;
    plot(tdata(1,:)*tmax,interpC(whichDataInd,:)); hold all;
    plot(tdata(1,:)*tmax,interpB(whichDataInd,:)); hold all;
    plot(tdata(1,:)*tmax,interpL(whichDataInd,:)); hold all;
    plot(tdata(1,:)*tmax,interpL(whichDataInd,:)+interpB(whichDataInd,:),'k'); hold all;
    plot(tdata(1,:)*tmax,diffusion(whichDataInd,:)); hold all;
    set(gca,'xlim',[0,20*60]);
    set(gca,'ylim',[exterior 1]);
    title(num2str(whichX(iplot)*xmax));
end
legend({'T','C','B','L','data'});
ax=axes('position',[0 0 1 1],'Visible','off','HitTest','off');
text(ax,0.5,1,['SE_{res}=' num2str(S)],'VerticalAlignment','top');
drawnow;


whichT=[10,30,60,120,600,1200]/tmax;
fig(2)=figure('position',get(0,'screensize')); 
for iplot=1:numel(whichT)
    subplot(2,3,iplot);
    [~,whichDataInd]=min(abs(tdata(1,:)-whichT(iplot)));
    plot(xdata(xdata(:,1)>=0,1)*xmax,interpT(xdata(:,1)>=0,whichDataInd)); hold all;
    plot(xdata(xdata(:,1)>=0,1)*xmax,interpC(xdata(:,1)>=0,whichDataInd)); hold all;
    plot(xdata(xdata(:,1)>=0,1)*xmax,interpB(xdata(:,1)>=0,whichDataInd)); hold all;
    plot(xdata(xdata(:,1)>=0,1)*xmax,interpL(xdata(:,1)>=0,whichDataInd)); hold all;
    plot(xdata(xdata(:,1)>=0,1)*xmax,interpL(xdata(:,1)>=0,whichDataInd)+interpB(xdata(:,1)>=0,whichDataInd),'k'); hold all;
    plot(xdata(:,1)*xmax,diffusion(:,whichDataInd)); hold all;
    set(gca,'ylim',[exterior 1]);
    title([num2str(whichT(iplot)*tmax)]);
end
legend({'T','C','B','L','data'});
ax=axes('position',[0 0 1 1],'Visible','off','HitTest','off');
text(ax,0.5,1,['SE_{res}=' num2str(S)],'VerticalAlignment','top');
drawnow;
drawnow;

fig(3)=figure('position',get(0,'screensize')); 
subplot(121);
imagesc(tdata(1,:)*tmax,xdata(:,1)*xmax,diffusion); colorbar;
set(gca,'ydir','normal');%,'ylim',[0,1000]);
set(gca,'clim',[exterior 1]);
subplot(122);
imagesc(tdata(1,:)*tmax,xdata(:,1)*xmax,interpT); colorbar;
set(gca,'ydir','normal');%,'ylim',[0,1000]);
set(gca,'clim',[exterior 1]);
ax=axes('position',[0 0 1 1],'Visible','off','HitTest','off');
text(ax,0.5,1,['SE_{res}=' num2str(S)],'VerticalAlignment','top');
drawnow;
   

function difference=fitfcn(p)    
    global xdata tdata diffusion tmax xmax optimizer xfittingRange tfittingRange beta doweighted
    global var_diffusion

    Length=3000;
    tspan=[0:1:9,10:5:20*60,1260:60:(60*30)]/tmax;
    xspan=[0:1:9,10:5:1200,1205:5:Length]/1000/xmax;

    [T,C,B,L,Baseline]=diffreactmodel(p,xspan,tspan);    
    [xgrid,tgrid]=ndgrid(xspan,tspan);
    valid=xdata>=xfittingRange(1)/xmax & xdata<=xfittingRange(2)/xmax & ...
        tdata>=tfittingRange(1)/tmax & tdata<=tfittingRange(2)/tmax;
    
    interp_fcn=scatteredInterpolant(tgrid(:),xgrid(:),T(:));
    interpT=interp_fcn(tdata,xdata);
    
    if doweighted
        weights=1;
    end        
    r=corrcoef((interpT(valid)),diffusion(valid));
    switch optimizer
        case 'fminunc'
            difference= 1-r(2);
            disp([num2str(p) ', ' num2str(sum(abs(tempdiffusion(valid)-interpT(valid)))) ', ' num2str(r(2))]);
        case 'lsqnonlin'
           difference= (diffusion(valid)-interpT(valid));
            disp([num2str(p) ', ' num2str(sum(abs(difference))) ', ' num2str(r(2))]);
    end
    if doweighted
        difference=difference.*weights(valid);
    end
    global fig doplotprogress
    if doplotprogress
        interp_fcn=scatteredInterpolant(tgrid(:),xgrid(:),C(:));
        interpC=interp_fcn(tdata,xdata);
        interp_fcn=scatteredInterpolant(tgrid(:),xgrid(:),B(:));
        interpB=interp_fcn(tdata,xdata);
        interp_fcn=scatteredInterpolant(tgrid(:),xgrid(:),Baseline(:));
        interpBaseline=interp_fcn(tdata,xdata);
        interp_fcn=scatteredInterpolant(tgrid(:),xgrid(:),L(:));
        interpL=interp_fcn(tdata,xdata);
        
        if ishandle(fig(1))
            close(fig(1));
        end
        if ishandle(fig(2))
            close(fig(2));
        end
        whichX=[100,200,400,800,1000,1200]/1000/xmax;
        fig(1)=figure('position',get(0,'screensize')); 
        for iplot=1:numel(whichX)
            subplot(2,3,iplot);
            [~,whichDataInd]=min(abs(xdata(:,1)-whichX(iplot)));
            plot(tdata(1,:)*tmax,interpT(whichDataInd,:)); hold all;
            plot(tdata(1,:)*tmax,diffusion(whichDataInd,:)); hold all;
            plot(tdata(1,:)*tmax,interpB(whichDataInd,:)); hold all;
            plot(tdata(1,:)*tmax,interpC(whichDataInd,:)); hold all;
            plot(tdata(1,:)*tmax,interpL(whichDataInd,:)); hold all;
            plot(tdata(1,:)*tmax,interpBaseline(whichDataInd,:)); hold all;
            set(gca,'xlim',[0,20*60]);
            set(gca,'ylim',[0 1]);
            title([num2str(whichX(iplot)*xmax)]);
        end
        drawnow;

        whichT=[60,120,600,1200]/tmax;
        fig(2)=figure('position',get(0,'screensize')); 
        for iplot=1:numel(whichT)
            subplot(2,3,iplot);
            [~,whichDataInd]=min(abs(tdata(1,:)-whichT(iplot)));
            plot(xdata(:,1)*xmax,interpT(:,whichDataInd)); hold all;
            plot(xdata(:,1)*xmax,diffusion(:,whichDataInd)); hold all;
            plot(xdata(:,1)*xmax,interpB(:,whichDataInd)); hold all;
            plot(xdata(:,1)*xmax,interpC(:,whichDataInd)); hold all;
            plot(xdata(:,1)*xmax,interpL(:,whichDataInd)); hold all;
            plot(xdata(:,1)*xmax,interpBaseline(:,whichDataInd)); hold all;
            set(gca,'ylim',[0 1]);
            title([num2str(whichT(iplot)*tmax)]);
        end    
        drawnow;
    end
end
            
function [T,C,B,L,Baseline]=diffreactmodel(p,x,t)
    global exterior C0 Length xmax Lipid0 thres tau tmax dogrowth expansioncoeffs IntensityNorm domasking
    
    VialDiameter=25000/1000;
    c2I=2.146307905734463e+02;
    
    rescaledExterior=exterior/VialDiameter/c2I;
    
    B0=0;
    D1= abs(p(1))*10^-1 *(tmax/(xmax)^2);
    kon1= abs(p(2))*tmax   *10^3;
    kunmask=abs(p(3))*tmax *10^2;

    Lipid0 =abs(p(4))*10^-3;    
    if domasking
        Masked0=abs(p(5))*10^-3;
    else
        Masked0=0;
    end
    
    C0=rescaledExterior;
    
    tau=1;    
    H=p(6);
    R1=4000/2/1000;
    theta=0;
    t0=0;

    R = @(B,t)  (100+B(1)-B(2)*exp(-B(3)*t*tmax))/100; 
    dRdt=@(B,t) B(2)*B(3)*tmax*exp(-B(3)*t*tmax)/100;
        
    density=@(x,t) 2*geometry(x,H,R1,theta);
        
    options=odeset('RelTol',1e-4,'AbsTol',1e-4,'NormControl','on','InitialStep',1e-7);
    m = 0;
    sol = pdepe(m,@oscpde,@icfun,@bcfun,x(:)',t0+t(:)',options);
        
    C=sol(:,(1:numel(x)),1)'*c2I;
    C=reshape(C,[],1);
    B=sol(:,(1:numel(x)),2)'*c2I;
    B=reshape(B,[],1);

    L=sol(:,(1:numel(x)),3)'*c2I;
    
    Baseline=sol(:,(1:numel(x)),5)'*c2I;
    Baseline=reshape(Baseline,[],1);

    T=(C+B+Baseline);
    T=T(:);

    function [c,f,s] = oscpde(x,t,u,dudx)
        D = D1;
        c = [1;1;1;1;1];
        f = [D*dudx(1);0;0;0;0];
        s = [-1;1;-1;0 ;0]*u(3)*u(1)*kon1;
        if domasking
            s = s+[0;0;1;-1;0]*kunmask*u(1)*u(4);
        end
                            
        if dogrowth
            if R(expansioncoeffs,t)>=1
                V   = -1*x*xmax*dRdt(expansioncoeffs,t);
                dVdx= -1  *xmax*dRdt(expansioncoeffs,t);
            else
                V=0;
                dVdx=0;
            end

            growth=tau*(dVdx*u+V*dudx)/(xmax);      
            s=s-growth;
        end
    end

    function u0 = icfun(x)
        if x<0
            u0 = [C0*2*R1;0;0;0;rescaledExterior*VialDiameter-C0*2*R1];
        elseif x==0
            u0 = [C0*2*R1;B0*density(x);Lipid0*density(x);Masked0*density(x);...
                rescaledExterior*VialDiameter-C0*2*R1];
        else
            u0 = [0;B0*density(x);Lipid0*density(x);Masked0*density(x);...
                rescaledExterior*VialDiameter-C0*density(x)];
        end
    end

    function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t) %#ok<INUSD>
        pL = [uL(1)-C0*2*R1;uL(2)-B0*density(xL);uL(3)-Lipid0*density(xL);uL(4)-Masked0*density(xL);...
            uL(5)-(rescaledExterior*VialDiameter-C0*2*R1)];
        qL = [0;0;0;0;0];
        pR = [uR(1)-0;      uR(2)-B0*density(xR);uR(3)-Lipid0*density(xR);uR(4)-Masked0*density(xR);...
            uR(5)-(rescaledExterior*VialDiameter-C0*density(xR))];
        qR = [1;1;1;1;1];
    end

    function l=geometry(x,H,R,theta)
        l=zeros(numel(x),1);
        for ix=1:numel(x)
            if xmax*x(ix)<=0
                l(ix)=0;
                continue;
            elseif xmax*x(ix)<=H
                l(ix)=sqrt(xmax*x(ix)./H.*(H.^2+R.^2)-(xmax*x(ix)).^2);
                continue;
            else
                l(ix)=(R+tan(theta).*H-tan(theta).*xmax*x(ix));
                continue;
            end
        end
    end
end    