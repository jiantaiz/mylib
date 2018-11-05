function PlotPulse(TimeHistFile)
[~,basename,~] = fileparts(TimeHistFile); 
basename = matlab.lang.makeValidName(['plot_', basename]);
mfile = [basename, '.m'];
result = perl('readTimeHist.pl',TimeHistFile,mfile);
run(mfile);