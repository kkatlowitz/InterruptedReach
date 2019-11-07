function ao=filter_rawdata(ao,pathbase,s,blocknum,varargin)

ao.fs_low=2000;
params.Fs=ao.fs_low;  

lfpall=resample(ao.dat,1,ao.fs/ao.fs_low);
lfpall=locdetrend(lfpall,ao.fs_low,[10,1]);
lfpall=rmlinesc(lfpall,params,0.05,'n');

if nargin>4
    SOS=varargin{1};
    G=varargin{2};
    lfpall=filtfilt(SOS,G,lfpall);
end

for ii=1:size(lfpall,2)
    lfp=lfpall(:,ii);
    lead=ii;
    save([pathbase,'/LFP/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(ii)],'lfp','s','blocknum','lead','-v7.3')
end
end
