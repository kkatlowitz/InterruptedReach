function write_dat(ao,tRange)
%this function will take in your neural data and turn it into a dat file
%for viewing and analysis. tRange is in seconds, use especially to cut off
%the end of the recording (like when we dont stop saving the file before
%the end of the experiment)
file=strrep(ao.file,'mat','dat');
fid=fopen(file,'w');
sig=ao.dat;
pad=3-size(sig,2);
sig=[sig,randn(size(sig,1),pad)/50];
if ~isempty(tRange)
    t=(1:size(sig,1))/ao.fs;
    inds=t>tRange(1)&t<tRange(2);
    sig=sig(inds,:);
end
fwrite(fid,int16(sig'*1e4),'int16');
fclose(fid);
