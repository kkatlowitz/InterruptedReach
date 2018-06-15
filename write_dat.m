function write_dat(loc,sub,block,ao,tRange)
%this function will take in your neural data and turn it into a dat file
%for viewing and analysis. tRange is in seconds, use especially to cut off
%the end of the recording (like when we dont stop saving the file before
%the end of the experiment)
filePrefix=[loc,sub,'\datFiles\block' num2str(block,'%02d')];
datFileLoc=[loc,sub,'\datFiles\'];
if ~exist(datFileLoc,'dir')
    mkdir(datFileLoc)
end
%raw
fid=fopen([filePrefix '_micro.dat'],'w');
sig=ao.micro;
pad=max(0,3-size(sig,1));
sig=[sig;randn(pad,size(sig,2))/1000];
if ~isempty(tRange)
    t=(1:size(sig,2))/ao.fs;
    inds=t>tRange(1)&t<tRange(2);
    sig=sig(inds,:);
end
fwrite(fid,int16(sig*1e4),'int16');
fclose(fid);

%lfp
fid=fopen([filePrefix '_micro.lfp'],'w');
sig=ao.lfp;
pad=max(0,3-size(sig,1));
sig=[sig;randn(pad,size(sig,2))/1000];
if ~isempty(tRange)
    t=(1:size(sig,2))/ao.fs_low;
    inds=t>tRange(1)&t<tRange(2);
    sig=sig(inds,:);
end
fwrite(fid,int16(sig*1e4),'int16');
fclose(fid);

%high pass
fid=fopen([filePrefix '_micro.fil'],'w');
sig=ao.bp;
pad=max(0,3-size(sig,1));
sig=[sig;randn(pad,size(sig,2))/1000];
if ~isempty(tRange)
    t=(1:size(sig,2))/ao.fs;
    inds=t>tRange(1)&t<tRange(2);
    sig=sig(inds,:);
end
fwrite(fid,int16(sig*1e5),'int16');
fclose(fid);

%macro
fid=fopen([filePrefix '_macro.dat'],'w');
sig=ao.macro;
pad=max(0,3-size(sig,1));
sig=[sig;randn(pad,size(sig,2))/1000];
if ~isempty(tRange)
    t=(1:size(sig,2))/ao.lfp_fs;
    inds=t>tRange(1)&t<tRange(2);
    sig=sig(inds,:);
end
fwrite(fid,int16(sig*1e5),'int16');
fclose(fid);



% d = designfilt('highpassfir','StopbandFrequency',50, ...
%   'PassbandFrequency',200,'StopbandAttenuation',65, ...
%   'PassbandRipple',.1,'SampleRate',ao.fs,'DesignMethod','kaiserwin');
% sig=filtfilt(d,double(sig)')';

