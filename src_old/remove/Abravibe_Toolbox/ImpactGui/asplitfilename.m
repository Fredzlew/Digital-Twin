function [Path,Prefix,Number,Ext] = asplitfilename(FileName)
% ASPLITFILENAME    Split file name in the form Path/Prefix<Number>.Ext
%                   into its parts
%
%       [Path,Prefix,Number,Ext] = asplitfilename(FileName)
%
%       Path        Full path to directory where file is located
%       Prefix      String beginning the file name. MUST NOT include . (periods)
%       Number      Number ending the file name, except for the extension
%       Ext         File extension, typically 'mat' or 'imptime' if valid
%                   ABRAVIVE file
%
%       FileName    String with one file name

% Start by splitting into path, file name, and extension
[Path,FullName,Ext] = fileparts(FileName);
% Next, split file name into Prefix and Number
m=length(FullName);
a=str2num(FullName(m));
while ~isempty(a) & m > 1 & isreal(a)   % if letter 'i' or 'j' in filename, imaginary number is returned!
    m=m-1;
    a=str2num(FullName(m));
end
Prefix=FullName(1:m);
Number=str2num(FullName(m+1:end));

