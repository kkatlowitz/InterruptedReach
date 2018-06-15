reboot;
addpath InterruptedReach
addpath(genpath('fcns\chronux'))
loc='C:\Users\akatl\Box Sync\Interrupted_Reach\Data\';
[num,txt,raw]=xlsread([loc,'interrupted_reach_data_summary.xlsx']);
% sub='SUB01';block=2;
% sub='SUB01';block=3;
% sub='SUB01';block=4;
% sub='SUB02';block=1;
% sub='SUB02';block=2;
% sub='SUB04';block=1;
% sub='SUB07';block=2;
% sub='SUB08';block=1;
% sub='SUB08';block=2;
sub='SUB08';block=3;
row=find(strcmp(txt(:,1),sub))-1+block;
F=[loc,sub,'\Block' num2str(block,'%02d') '\'];
ao_file=[raw{row,8},'.mat'];
behavioral_file=[raw{row,9},'.mat'];
ao=loadNeuralData(F,ao_file);
ao=filter_rawdata(ao);%this is more for visualization if you want to look at the data
write_dat(loc,sub,block,ao,[])%write to dat file. use http://neurosuite.sourceforge.net/ to view data easily
%%
% Load Data
ao=loadNeuralData(F,ao_file);
b=loadBehavior(F,ao_file,behavioral_file);

%%
%% Analyze Behavior
b=calculate_knee(b);%this need to be made
%% Spike Sort
% This part is still not solidified yet, so expect it to change alot
%klustakwik/kilosort section

% wav_clus section
% [spikes,thr,index] = spike_detect(x, par ,thr);%if you want single
[spk,index] = spike_adapt(signal,par);%if you want adaptive
%WAV CLUS HERE
%% no matter what we use above, the last steo will be the same: loading in
% the spike times and clusters
fpos=[];
spikes=load_spikes([F,'POS/' fpos],'POS');%give the file, and then say what algorithm you used to get the spikes.
%output will be a 2d mat: first dimension is spike times, second is spike
%cluser.
%% PSTH
b=align_spikes(b,spikes);%this loads the wav_clus data and aligns it to events
