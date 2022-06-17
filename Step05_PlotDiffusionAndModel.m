global xdata tdata tmax xmax
global dounmask dogrowth 

depthrange=[-280,1500]/1000;
baselinerange=[-150,-50]/1000;

xvalidationRange=[0,1000]/1000;

tmax=1;%60*20;
xmax=1;

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

interp_tdata=0:5:min(1200,max(raw_tdata(:)));
interp_xdata=(baselinerange(1):0.005:depthrange(2));
tdata=ones(numel(interp_xdata),1)*interp_tdata;
xdata=interp_xdata'*ones(1,numel(interp_tdata));

interp_fcn=scatteredInterpolant(raw_tdata(:),raw_xdata(:),cdiffusion(:));
diffusion=interp_fcn(tdata,xdata);

diffusion(:,tdata(1,:)>floor(max(raw_tdata(:))))=nan;

baseline=xdata(:,1)'<=baselinerange(2) & xdata(:,1)'>=baselinerange(1);
exterior=nanmean(reshape(diffusion(baseline,:),[],1)); %#ok<*NANMEAN>


Imax=max(max(medfilt2(diffusion,[5,5])));

xdata=xdata/xmax;
tdata=tdata/tmax;
tspan=[0:1:9,10:5:(60*30)]/tmax;
xspan=[0:1:9,10:5:3000]/1000/xmax;
global factor
factor=1;

domodel=false;
if exist(modelfile,'file')
    domodel=true;
    clear p exterior expansioncoeffs;
    load(modelfile,'p','exterior','expansioncoeffs');    

    [T,C,B,L]=diffreactmodel(p,xspan,tspan,exterior,expansioncoeffs);
    [xgrid,tgrid]=ndgrid(xspan,tspan);
    interp_fcn=scatteredInterpolant(tgrid(:),xgrid(:),T(:));
    interpT=interp_fcn(tdata,xdata);

    valid=xmax*xdata>=xvalidationRange(1) & xmax*xdata<=xvalidationRange(2);
    a=interpT; b=diffusion;
    r=corrcoef((a(valid)),b(valid));
    SSres=sum((a(valid)-b(valid)).^2);

    n=numel(b(valid));
    k=5;
    Radj=1-((1-r(2)^2)*(n-1)/(n-k-1));
    df=n-k-1;
    %Residual Standard Error
    S=sqrt(SSres/df);
end

fig(1)=figure('position',get(0,'screensize')); 
whichX=([200,400,600,800,1000,1200])/xmax/1000;
for iplot=1:numel(whichX)
    subplot(3,6,iplot);
    [~,whichDataInd]=min(abs(xdata(:,1)-whichX(iplot)));
    plot(tdata(1,:)*tmax,squeeze(diffusion(whichDataInd,:)),'b'); hold all;
    if domodel
        plot(tdata(1,:)*tmax,squeeze(interpT(whichDataInd,:)),'k'); hold all;
    end
    set(gca,'xlim',[0,20*60]);
    set(gca,'ylim',[exterior Imax]*factor);
    ylabel('Intensity (a.u.)')
    xlabel('Time t (min)');
    title(['x=' num2str(whichX(iplot)*xmax*1000) ' ' char(181) 'm']);
    drawnow;
end
whichT=[30,60,120,360,600,1200]/tmax;
for iplot=1:numel(whichT)
    subplot(3,6,6+iplot);
    [~,whichDataInd]=min(abs(tdata(1,:)-whichT(iplot)));
    valid=xdata(:,1)>=0;
    plot(xdata(:,1)*xmax*1000,squeeze(diffusion(:,whichDataInd)),'b'); hold all;

    if dodetectpeaks
        peakRange=[0.7/xmax,max(xdata(:))]; %#ok<UNRCH>
        validpeakrange=xdata(:,1)>=peakRange(1) & xdata(:,1)<=peakRange(2);
        startpos=find(validpeakrange,1,'first');
        [maxval,maxind]=max(squeeze(diffusion(validpeakrange,whichDataInd)));
        [minval]=min(squeeze(diffusion(startpos-1+(1:maxind),whichDataInd)));
        (maxval-minval)/minval
        if (maxval-minval)/minval>=0.0225
            plot(xdata(maxind-1+startpos,1)*xmax*1000,...
                squeeze(diffusion(maxind-1+startpos,whichDataInd)),...
                'ko','markersize',12);
        end
    end

    if domodel
        plot(xdata(valid,1)*xmax*1000,squeeze(interpT(valid,whichDataInd)),'k'); hold all;
    end
    set(gca,'ylim',[exterior,Imax]*factor);
    set(gca,'xlim',[0,1200]);
    ylabel('Intensity (a.u.)')
    xlabel(['Depth x (' char(181) 'm)']);
    title(['t=' num2str(whichT(iplot)*tmax) ' min']);
    drawnow;
end

subplot(3,6,12+[1 2]);
imagesc(tdata(1,:)*tmax,xdata(:,1)*xmax*1000,diffusion); hc=colorbar; set(get(hc,'label'),'string','Intensity (a.u.)');
set(gca,'ydir','normal','ylim',[0,1200],'xlim',[0,1200]); axis square;
set(gca,'xtick',[0,400,800,1200],'ytick',[0,400,800,1200]);
ylabel(['Depth x (' char(181) 'm)']);
xlabel('Time t (min)');
set(gca,'clim',[exterior Imax]*factor);
title('data');
drawnow;
if domodel
    subplot(3,6,12+[3 4]);
    imagesc(tdata(1,:)*tmax,xdata(:,1)*xmax*1000,interpT); hc=colorbar; set(get(hc,'label'),'string','Intensity (a.u.)');
    set(gca,'ydir','normal','ylim',[0,1200],'xlim',[0,1200]); axis square;
    set(gca,'xtick',[0,400,800,1200],'ytick',[0,400,800,1200]);
    ylabel(['Depth x (' char(181) 'm)']);
    xlabel('Time t (min)');
    set(gca,'clim',[exterior Imax]*factor);
    title(['model: R^2_{adj}=' num2str(Radj) ', SE_{res}=' num2str(S)]);
    drawnow;
end

temp=load(densityfile);
ReferenceValue=10000; 

subplot(3,6,12+5);
im=double(temp.firstIm);
im=ReferenceValue-im;
im=im-baselinecorrection(1);
im=(im-IntensityRef)/IntensityNorm;
imagesc(im); axis image;
axis off; axis tight;
set(gca,'clim',[exterior Imax]*factor);
title(['t=' num2str(temp.firstTime) 'min']);
subplot(3,6,12+6);
im=double(temp.lastIm);
im=ReferenceValue-im;
im=im-baselinecorrection(end);
im=(im-IntensityRef)/IntensityNorm;
imagesc(im); axis image;
axis off; %axis tight;
set(gca,'clim',[exterior Imax]*factor);
title(['t=' num2str(temp.lastTime) 'min']);
drawnow;

function [T,C,B,L,Baseline]=diffreactmodel(p,x,t,exterior,expansioncoeffs)
    global C0 xmax Lipid0 tau tmax dogrowth dounmask
    
    VialDiameter=25000/1000;
    c2I=2.146307905734463e+02;
    
    rescaledExterior=exterior/VialDiameter/c2I;
    
    B0=0;
    D1= abs(p(1))*10^-1 *(tmax/(xmax)^2);
    kon1= abs(p(2))*tmax   *10^3;
    kunmask=abs(p(3))*tmax *10^2;

    Lipid0 =abs(p(4))*10^-3;    
    if dounmask
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
        if dounmask
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