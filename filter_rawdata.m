function ao=filter_rawdata(ao)
fs=ao.fs;
bpFilt = designfilt('bandpassfir','FilterOrder',1000, ...
    'CutoffFrequency1',250,'CutoffFrequency2',10000, ...
    'SampleRate',fs);
ao.bp=filtfilt(bpFilt,ao.dat);
ao.lfp=resample(ao.dat,1,46);%switch this to a median filter?
ao.fs_low=ao.fs/44;