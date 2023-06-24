function rootPath = simseqRootPath()
%
% Return the path to root simSeq project folder
%
% This function must reside in the directory base of this code repository.
% It is used to determine the location of various subdirectories
%
% Example:
%   fullfile(simseqRootPath, 'stimulus')
%
% By Eline Kupers 2020


rootPath = which('simseqRootPath');
rootPath = fileparts(rootPath);

return

