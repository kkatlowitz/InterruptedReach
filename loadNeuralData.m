function ao=loadNeuralData(F,ao_file)
%% Load alpha omega neural data
ao=struct();
ao.file=[F,ao_file];%store the data file for ref
if exist(ao.file,'file')
    variableInfo = who('-file', ao.file);
else
    error('File Not Found. Check yourself before you wreck yourself')
end
load([F,ao_file],'CRAW_01_KHz','CRAW_01_TimeBegin');%you need these for both
if ismember('CRAW_01', variableInfo) 
    load([F,ao_file],'CRAW_01','CRAW_01_BitResolution','CRAW_01_Gain');
    lead1=double(CRAW_01)./(CRAW_01_BitResolution*CRAW_01_Gain);
else
    lead1=[];
    disp('No Lead 1')
end
if ismember('CRAW_02', variableInfo) 
    load([F,ao_file],'CRAW_02','CRAW_02_BitResolution','CRAW_02_Gain');
    lead2=double(CRAW_02)./(CRAW_02_BitResolution*CRAW_02_Gain);
else 
    disp('No Lead 2')
    lead2=[];
end
%convert raw data on micro electrodes
ao.fs=CRAW_01_KHz*1000;%sampling rate
ao.micro=[lead1;lead2];

if ismember('CMacro_LFP_01', variableInfo)
    load([F,ao_file],'CMacro_LFP_01','CMacro_LFP_01_KHz','CMacro_LFP_01_BitResolution','CMacro_LFP_01_Gain');
    macro1=double(CMacro_LFP_01)./(CMacro_LFP_01_BitResolution*CMacro_LFP_01_Gain);%these numbers should be the same
    ao.macro_fs=CMacro_LFP_01_KHz*1e3;
else
    macro1=[];
    disp('No Macro 1')
end
if ismember('CMacro_LFP_02', variableInfo) 
    load([F,ao_file],'CMacro_LFP_02','CMacro_LFP_02_KHz','CMacro_LFP_02_BitResolution','CMacro_LFP_02_Gain');
    macro2=double(CMacro_LFP_02)./(CMacro_LFP_02_BitResolution*CMacro_LFP_02_Gain);
    ao.macro_fs=CMacro_LFP_02_KHz*1e3;%should overwrite with same thing if we have both
else
    macro2=[];
    disp('No Macro 2')
end
ao.macro=[macro1;macro2];