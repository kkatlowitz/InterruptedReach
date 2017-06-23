function ao=filter_rawdata(ao)
fs=ao.fs;
ao.fs_low=ao.fs/44;
bpFilt = designfilt('bandpassfir','FilterOrder',500, ...
    'CutoffFrequency1',250,'CutoffFrequency2',10000, ...
    'SampleRate',fs);
lfpFilt = designfilt('bandpassfir','FilterOrder',500, ... %% for LFP
    'CutoffFrequency1',0.1,'CutoffFrequency2',50, ...
    'SampleRate',ao.fs_low);

%lpFilt = designfilt('bandstopfir','FilterOrder',500, ... %% use to get rid of slow oscillations
%    'CutoffFrequency1',0.5,'CutoffFrequency2',3, ...
%    'SampleRate',ao.fs_low);

ao.bp=filtfilt(bpFilt,ao.dat);
ao.dat=resample(ao.dat,1,44); %% for syncing to LFP
ao.lfp=filtfilt(lfpFilt,ao.dat);
%ao.lfp=filtfilt(lpFilt,ao.lfp);  %% use to get rid of slow oscillations
