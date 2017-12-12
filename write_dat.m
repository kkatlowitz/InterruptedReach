function write_dat(ao,tRange)
%this function will take in your neural data and turn it into a dat file
%for viewing and analysis. tRange is in seconds, use especially to cut off
%the end of the recording (like when we dont stop saving the file before
%the end of the experiment)
file=strrep(ao.file,'mat','dat');
file=strrep(file,'MatFiles','dat');
fid=fopen(file,'w');
sig=ao.dat;
pad=max(0,3-size(sig,1));
sig=[sig;randn(pad,size(sig,2))/1000];
if ~isempty(tRange)
    t=(1:size(sig,2))/ao.fs;
    inds=t>tRange(1)&t<tRange(2);
    sig=sig(inds,:);
end

d = designfilt('highpassfir','StopbandFrequency',50, ...
  'PassbandFrequency',200,'StopbandAttenuation',65, ...
  'PassbandRipple',.1,'SampleRate',ao.fs,'DesignMethod','kaiserwin');
% sig=filtfilt(d,double(sig)')';

fwrite(fid,int16(sig*1e5),'int16');
fclose(fid);
