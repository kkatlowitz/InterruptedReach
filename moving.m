function [mov_time,mov_ind]=moving(coords,t,start_t,lpfilt,time_int)

pad=floor(length(lpfilt)/2);
t=t-start_t;
tinterp=0:0.01:max(t);
xinterp=interp1(t,coords(:,1),tinterp);
xinterp=filter(lpfilt,1,xinterp);
xinterp(end+1:end+pad)=NaN;
xinterp=xinterp(pad+1:end);
xvel=diff(xinterp)./diff(tinterp);
xacc=diff(xvel)./diff(tinterp(1:end-1));
xacc=filter(lpfilt,1,xacc);
xacc(end+1:end+pad)=NaN;
xacc=xacc(pad+1:end);

go_ind=discretize(time_int,tinterp);

xvel(1:go_ind(1)-1)=NaN;

if ~isnan(go_ind(2))
    xvel(go_ind(2)+1:end)=NaN;
end
[peaks,peaktime]=findpeaks(abs(xvel));
[peaks(end+1),peaktime(end+1)]=max(abs(xvel));
peaktime=peaktime(find(peaks>300,1));


[peakacc,peakacctime]=findpeaks(abs(xacc));
[peakacc(end+1),peakacctime(end+1)]=max(abs(xacc));
%[mins,mintime]=findpeaks(-abs(xvel));

%[mins(end+1),mintime(end+1)]=max(-abs(xvel));



if ~isempty(peaktime)
    %peakdiff=peaktime-mintime;
    peakdiff=peaktime-peakacctime;
    peakdiff(peakdiff<0)=NaN;
    
    [~,mindiffind]=min(peakdiff);
    %mov_ind=mintime(mindiffind);
    mov_ind=peakacctime(mindiffind);
    mov_time=tinterp(mov_ind);
    [~,mov_ind]=min(abs(t-mov_time));
    mov_time=t(mov_ind);
else
    mov_ind=NaN;
    mov_time=NaN;
end

end
