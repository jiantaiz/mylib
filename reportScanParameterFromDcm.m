function res = reportScanParameterFromDcm(dinfo)
% REPORTSCANPARAMETERFROMDCM ...
%    res = REPORTSCANPARAMETERFROMDCM(dinfo) ...
%
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 04/16/2019
%%
if ~isstruct(dinfo)
    dinfo = dicominfo(dinfo);
    
end
res.pat = dinfo.PatientName;
res.exno = dinfo.StudyID;
res.seno = dinfo.SeriesNumber;
res.sedesc = dinfo.SeriesDescription;
res.sedate = dinfo.SeriesDate;
res.TR = dinfo.RepetitionTime; %TR in msec
res.TE = dinfo.EchoTime;
res.ESP = dinfo.Private_0043_102c; %echo spaceing in usec
res.FOV = dinfo.ReconstructionDiameter; %in mm
ACQMTX = dinfo.AcquisitionMatrix; %acqusition matrix;
res.ACQMTX = ACQMTX(ACQMTX>0)';
res.PixSize = res.FOV./res.ACQMTX;
res.ScanTime = dinfo.Private_0019_105a; %AcquisitionDuration in usec
res.ScanTimeInMinute = res.ScanTime/1e6/60; %in minutes
res.PhaseOffsets = dinfo.NumberOfTemporalPositions;
res.PulseSeq = dinfo.Private_0019_109c;

group = '0019';
elem = hex2dec('10A7');
for k=1:22;
    fld = dicomlookup(group,elem+k-1);
    res.userCV(k) = dinfo.(fld);
end
res.CoilElements = res.userCV(3);
res.DriverAmp = res.userCV(15);
res.MEGFreq = res.userCV(14);
res.VibrationFreq = res.userCV(22);






