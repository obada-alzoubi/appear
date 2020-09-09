% import parameters
function [alphacal,dirMATLAB,outform] = cmdin(filename)
	% The function will read settings file for the code 
	% Input is text file with attribute and value pairs 
	% Please use the demo for more information about the strcuture 

	if(exist(filename,'file')~=2)
	    error(sprintf('%s does not exist',filename));
end
	fp = fopen(filename);
	formatSpec = 'AlphaPower = %s MATLABdir = %s OutputFormat = %s';
	C = textscan(fp,formatSpec,'Delimiter','\n','CollectOutput',true);
	fclose(fp);
	alphacal = C{1}{1};
	dirMATLAB = C{1}{2};
	outform = C{1}{3};
	if(exist(strcat(dirMATLAB,'/MATLAB'),'dir')~=7)
	    error(sprintf('%s/MATLAB does not exist',dirMATLAB));
	end
end

