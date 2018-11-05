function str = readString(fp)
% READSTRING reads in a null-terminated string from file
%    str = READSTRING(fp) returns a null-terminated string str from
%    file handler fp
%
%    See also: readPdFile

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017
%%
%% 
str =fread(fp,1,'uint8=>char');

while str(end) ~= 0
    char = fread(fp,1,'uint8=>char');
    str = [str char];
end
str = str(1:end-1);
