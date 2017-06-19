reboot;
addpath InterruptedReach
F='Data/Sub01/';
ao_file='RT1D2.682F0002.mat';
behavioral_file='Sub01_DD_04192017_OR_HandL_STNR_Block02.mat';
fpos='Sub01_DD_04192017_OR_HandL_STNR_Block02_sorted.csv';
%% Load Data
[ao,b]=extractData(F,ao_file,behavioral_file,[1,2],[]);
% ao=filter_rawdata(ao);%not needed, spike detection does filtering for you.this is more for visualization if you want to look at the data
%% Analyze Behavior
b=calculate_knee(b);%this need to be made
%% Spike Sort
% This part is still not solidified yet, so expect it to change alot

%klustakwik/kilosort section
write_dat(ao,[0,535])%write to dat file. use http://neurosuite.sourceforge.net/ to view data easily

% wav_clus section
% [spikes,thr,index] = spike_detect(x, par ,thr);%if you want single
[spk,index] = spike_adapt(signal,par);%if you want adaptive
%WAV CLUS HERE

% no matter what we use above, the last steo will be the same: loading in
% the spike times and clusters
spikes=load_spikes([F,'POS/' fpos],'POS');%give the file, and then say what algorithm you used to get the spikes.
%output will be a 2d mat: first dimension is spike times, second is spike
%cluser.
%% PSTH
b=align_spikes(b,spikes);%this loads the wav_clus data and aligns it to events
