function [spikes,index] = spike_adapt(signal,par)

indlist=1:length(signal); %all the indices of the signal
subs=ceil(indlist/176000); %segments of 4 seconds
mins=accumarray(subs',signal',[],@min); %find the minimim of each segment
breakpoints=find(abs(diff(mins)./mins(1:end-1))>0.25) * 176000;
if ~isempty(breakpoints)
    split=[breakpoints(1);diff(breakpoints);length(signal)-breakpoints(end)];
    signalbreak=mat2cell(signal,split);
    indexadd=cumsum(split);
else
    signalbreak{1}=signal;
end
for j=1:length(signalbreak)
    signal=signalbreak{j};
    i=1;
    [spikestemp,thr{j}(i),indextemp]=spike_detect(signal',par,[]);
    while ~isempty(spikestemp)
        i=i+1;
        indmin=indextemp-19;
        deletion=reshape(bsxfun(@(x,y) x+y,indmin,(0:63)'),1,[]);
        deletion(deletion>length(signal))=[];
        signal(deletion)=[];
        [spikestemp,thr{j}(i),indextemp]=spike_detect(signal',par,[]);
    end
    
    [spikes{j},~,index{j}]=spike_detect(signalbreak{j}',par,thr{j}(end));
    if j>1
        index{j}=index{j}+breakpoints(j-1);
    end
end

spikes=spikes';
spikes=vertcat(spikes{:});

index=horzcat(index{:});

end
