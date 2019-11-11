function [R,t,E,pmat,RR_ind] = psth(data,bounds,sig,plt,T,err,t,type,N,p,comp)
%  function to plot trial averaged rate smoothed by 
%  a Gaussian kernel - visual check on stationarity 
%  Usage: [R,t,E] = psth(data,sig,plt,T,err,t)
%                                                   
%  Inputs:                                  
% Note that all times have to be consistent. If data
% is in seconds, so must be sig and t. If data is in 
% samples, so must sig and t. The default is seconds.
% data            structural array of spike times   
% sig             std dev of Gaussian (default 50ms)
%                 (minus indicates adaptive with    
%                  approx width equal to mod sig)   
% plt = 'n'|'r' etc      (default 'r')              
% T is the time interval (default all)              
% err - 0 = none                                    
%       1 = Poisson                                 
%       2 = Bootstrap over trials (default)         
% (both are based on 2* std err rather than 95%)    
% t   = times to evaluate psth at                   
%                                                   
% The adaptive estimate works by first estimating   
% the psth with a fixed kernel width (-sig) it      
% then alters the kernel width so that the number   
% of spikes falling under the kernel is the same on 
% average but is time dependent.  Reagions rich     
% in data therefore have their kernel width reduced 
%                                                    
% Outputs:                                 
%                                                   
% R = rate                                          
% t = times                                         
% E = errors (standard error)                       

if nargin <1 ; error('I need data!');end
[data]=padNaN(data); % create a zero padded data matrix from input structural array
sz=size(data);
% transposes data so that the longer dimension is assumed to be time
% Procedure used to bring input data into form compatible with that required
% for Murray's routine
if sz(1)>sz(2); data=data';end; 
if length(type)==size(data,2)
    data=data';
end
if nargin < 2; sig = 0.05;end
if nargin < 3; plt = 'r';end
if nargin < 5; err = 2; end

if isempty(sig); sig = 0.05;end
if isempty(plt); plt = 'r';end
if isempty(err); err = 2; end

adapt = 0;
if sig < 0; adapt = 1; sig = -sig; end

%  to avoid end effects increase time interval by the width
%  of the kernel (otherwise rate is too low near edges)

if nargin < 4; 
  T(1) = min(data(:,1));
  T(2) = max(max(data));
else
  T(1) = T(1)-4*sig;
  T(2) = T(2)+4*sig;
end

% calculate NT and ND and filter out times of interest

NT = length(data(:,1));
if NT < 4 && err == 2
  disp('Using Poisson errorbars as number of trials is too small for bootstrap')
  err = 1;
end    
    
m = 1;
D = zeros(size(data));
ND=zeros(1,NT);
for n=1:NT
  indx = find(~isnan(data(n,:)) & data(n,:)>=T(1) & data(n,:)<=T(2));
  ND(n) = length(indx);
  D(n,1:ND(n)) = data(n,indx);
  m = m + ND(n); 
end
N_tot = m;
N_max = max(ND);
D = D(:,1:N_max);

% if the kernel density is such that there are on average 
% one or less spikes under the kernel then the units are probably wrong

L = N_tot/(NT*(T(2)-T(1)));
if 2*L*NT*sig < 1 || L < 0.1 
  disp('Spikes very low density: are the units right? is the kernel width sensible?')
  disp(['Total events: ' num2str(fix(100*N_tot)/100) ' sig: ' ...
        num2str(fix(1000*sig)) 'ms T: ' num2str(fix(100*T)/100) ' events/sig: ' ...
        num2str(fix(100*N_tot*sig/(T(2)-T(1)))/100)])
end

%    Smear each spike out  
%    std dev is sqrt(rate*(integral over kernal^2)/trials)     
%    for Gaussian integral over Kernal^2 is 1/(2*sig*srqt(pi))

if nargin < 6
  N_pts =  fix(5*(T(2)-T(1))/sig);
  t = linspace(T(1),T(2),N_pts);
else
  N_pts = length(t);
end
  
bound_disc=discretize(bounds,t);
bound_disc(bounds<t(1))=1;
bound_disc(bounds>t(end))=length(t);

RR = NaN*ones(NT,N_pts);
f = 1/(2*sig^2);
for n=1:NT
  for m=1:ND(n)
    RR(n,bound_disc(n,1):bound_disc(n,2)) = sum([RR(n,bound_disc(n,1):bound_disc(n,2));exp(-f.*(bsxfun(@minus,t(bound_disc(n,1):bound_disc(n,2)),D(n,m))).^2)],'omitnan');
  end
end
RR = RR*(1/sqrt(2*pi*sig^2));

types=unique(type);
R=zeros(length(types),size(RR,2));
if NT > 1
    for ii=1:length(types)
        R(ii,:) = mean(RR(type==types(ii),:),'omitnan');
    end
else
    for ii=1:length(types)
        R(ii,:) = RR;
    end
end

for ii=1:length(types)
    RR_ind{ii}=RR(type==types(ii),:);
end

E = zeros(length(types),length(R));
for ii=1:length(types)
    NTvar=sum(~isnan(RR(type==types(ii),:)));
    err1=NTvar<4;
    
    if err == 1
        E(ii,:) = sqrt(R(ii,:)./(2*sig*sqrt(pi).*NTvar));
    elseif err == 2
        Nboot = 10;
        mE = 0;
        sE = 0;
        for b=1:Nboot
            indx = floor(NT*rand(1,NT)) + 1;
            mtmp = mean(RR(indx,~err1),'omitnan');
            mE = mE + mtmp;
            sE = sE + mtmp.^2;
        end
        E(ii,~err1) = sqrt((sE/Nboot - mE.^2/Nboot^2));
        E(ii,err1) = sqrt(R(ii,err1)./(2*sig*sqrt(pi).*NTvar(err1)));
    end
end

sigbounds=round(median(bound_disc));

fr1=NaN*ones(1,sigbounds(2)-sigbounds(1)+1);
fr2=fr1;

fr1=R(comp(1),sigbounds(1):sigbounds(2));
fr2=R(comp(2),sigbounds(1):sigbounds(2));

n1=nnz(type==types(comp(1)));
n2=nnz(type==types(comp(2)));
n=n1+n2;

alltrials=1:n;

N2=N;

randfr1=zeros(N,length(fr1));
randfr2=randfr1;

RR=RR(type==types(comp(1)) | type==types(comp(2)),:);

for ii=1:N
    randtrials1=randperm(n,n1);
    randtrials2=find(~ismember(alltrials,randtrials1));
    
    randfr1(ii,:)=mean(RR(randtrials1,sigbounds(1):sigbounds(2)),'omitnan');
    randfr2(ii,:)=mean(RR(randtrials2,sigbounds(1):sigbounds(2)),'omitnan');
end

iter=randperm(N,N2);

frdiff=fr2-fr1;

randfrdiff=randfr2-randfr1;
randfrdiffsort=sort(randfrdiff,'ascend');

poscomp=randfrdiffsort(round((1-p)*N),:);
negcomp=randfrdiffsort(round(p*N),:);

sig_ind_pos=frdiff>poscomp & frdiff>0;
sig_ind_neg=frdiff<negcomp & frdiff<0;

clusts_pos=regionprops(sig_ind_pos,'PixelList');
clusts_neg=regionprops(sig_ind_neg,'PixelList');

clusts=[clusts_pos;clusts_neg];

z=(frdiff-mean(randfrdiff))./std(randfrdiff);
zrand=bsxfun(@times, bsxfun(@minus, randfrdiff, mean(randfrdiff)), 1./std(randfrdiff));

%frall=mean(RR(:,sigbounds(1):sigbounds(2)),'omitnan');
%stdall=std(RR(:,sigbounds(1):sigbounds(2)),'omitnan');
%z=(frdiff-frall)./stdall;
%zrand=bsxfun(@times, bsxfun(@minus, randfrdiff, frall), 1./stdall);

zclust=zeros(1,length(clusts));
clustscell=cell(length(clusts),1);
for ii=1:length(clusts)
    zclust(ii)=sum(arrayfun(@(x) z(x),clusts(ii).PixelList(:,1)));
    clustscell{ii}=clusts(ii).PixelList(:,1);
end

sig_ind_pos_rand=bsxfun(@gt,randfrdiff,poscomp) & randfrdiff>0;
sig_ind_neg_rand=bsxfun(@lt,randfrdiff,negcomp) & randfrdiff<0;

zclustrandmax=zeros(1,N2);

for j=1:N2
    randfrdifftemp=randfrdiff(iter(j),:);
    
    clusts_pos_rand=regionprops(sig_ind_pos_rand(iter(j),:),'PixelList');
    clusts_neg_rand=regionprops(sig_ind_neg_rand(iter(j),:),'PixelList');
    
    clusts_rand=[clusts_pos_rand;clusts_neg_rand];
    
    if ~isempty(clusts_rand)
        zclustrand=zeros(1,length(clusts_rand));
        for ii=1:length(clusts_rand)
            zclustrand(ii)=sum(arrayfun(@(x) zrand(x),clusts_rand(ii).PixelList(:,1)));
        end
        zclustrandmax(j)=max(abs(zclustrand));
    else
        zclustrandmax(j)=0;
    end
    
end

pclust=1-sum(bsxfun(@gt,abs(zclust)',zclustrandmax),2)/N2;
pmat=NaN.*ones(1,size(RR,2));

for ii=1:length(pclust)
    pmat(clusts(ii).PixelList(:,1)+sigbounds(1)-1)=pclust(ii);
end


% if adaptive warp sig so that on average the number of spikes
% under the kernel is the same but regions where there is 
% more data have a smaller kernel

if adapt 
  sigt = mean(R,'omitnan')*sig./R;
  RR = zeros(NT,N_pts);
  f = 1./(2*sigt.^2);
  for n=1:NT
    for m=1:ND(n)
      RR(n,bound_disc(n,1):bound_disc(n,2)) = sum([RR(n,bound_disc(n,1):bound_disc(n,2));exp(-f(bound_disc(n,1):bound_disc(n,2)).*(bsxfun(@minus,t(bound_disc(n,1):bound_disc(n,2)),D(n,m))).^2)],'omitnan');
    end
    RR(n,:) = RR(n,:).*(1./sqrt(2*pi*sigt.^2));
  end
  
  if NT > 1
      R = mean(RR,'omitnan');
  else
      R = RR;
  end

  NTvar=sum(~isnan(RR));
  
  if err == 1
      E = sqrt(R./(2*sigt*sqrt(pi).*NTvar));
  elseif err == 2
    Nboot = 10;
    mE = 0;
    sE = 0;
    for b=1:Nboot
      indx = floor(NT*rand(1,NT)) + 1;
      mtmp = mean(RR(indx,~err1),'omitnan');
      mE = mE + mtmp;
      sE = sE + mtmp.^2;
    end
    E(~err1) = sqrt((sE/Nboot - mE.^2/Nboot^2));
    E(err1) = sqrt(R(err1)./(2*sig*sqrt(pi).*NTvar(err1)));
  end
end  


if plt == 'n';return;end
plot(t,R,plt)
hold on
if err > 0
  plot(t,R+2*E,'g')
  plot(t,R-2*E,'g')
end
%axis([T(1)+(4*sig) T(2)-(4*sig) 0 1.25*max(R)])
axis([T(1)+(4*sig) T(2)-(4*sig) 0 max(R+2*E)+10])
xlabel('time (s)')
ylabel('rate (Hz)')
title(['Trial averaged rate : Gaussian Kernel :'  ...
	    ' sigma = ' num2str(1000*sig) 'ms'])
hold off











