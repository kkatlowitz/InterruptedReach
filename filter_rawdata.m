function ao=filter_rawdata(ao)
fs=ao.fs;
ao.fs_low=1375;%should be 1000

bpFilt = designfilt('bandpassfir','FilterOrder',500, ...
    'CutoffFrequency1',250,'CutoffFrequency2',10000, ...
    'SampleRate',fs);
lfpFilt = designfilt('bandpassfir','FilterOrder',500, ... %% for LFP
    'CutoffFrequency1',0.1,'CutoffFrequency2',50, ...
    'SampleRate',ao.fs_low);
ao.bp=filtfilt(bpFilt,ao.micro')';
for i=1:size(ao.bp,1)
    ao.lfp(i,:)=filtfilt(lfpFilt,resample(ao.micro(i,:),1,32));
end
if fs~=44000
    error('fix sampling rates')
end
%dont do
%lpFilt = designfilt('bandstopfir','FilterOrder',500, ... %% use to get rid of slow oscillations
%    'CutoffFrequency1',0.5,'CutoffFrequency2',3, ...
%    'SampleRate',ao.fs_low);
%ao.lfp=filtfilt(lpFilt,ao.lfp);  %% use to get rid of slow oscillations
