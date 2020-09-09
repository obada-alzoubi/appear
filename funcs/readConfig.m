function [Config] = readConfig(file_in)
% Author: Obada Al Zoubi oabad.y.alzoubi@gmail.com
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
% Read config file
f      = fopen(file_in,'r');                % open file
Config = []; 
while ~feof(f)  
    if isempty(f)
        continue;
    end;
    lin         = strtrim(fgetl(f)) ;   
    lin         = strtrim(lin); % Remove un-needed spaces
    prop_val    = strsplit(lin, '=');
    prop        = strtrim(prop_val{1});   
    val         = strtrim(prop_val{2});
    Config.(prop) = val;    
end
end

