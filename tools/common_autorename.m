function[SaveName] = common_autorename(PrevName, outputfolder,rootfolder)
% This function makes sure no run is overwritten.
% Check whether filename already exists, rename if necessary.
%
% Code by S. Delle Chiaie and F. Kurcz
%

[fPath, fName, fExt] = fileparts(PrevName);

if isempty(fExt)  % append .mat extension if there is none
  fExt     = '.mat';
  PrevName = fullfile(fPath, [fName, fExt]);
end

if exist([outputfolder,PrevName],'file')
    cd(outputfolder)                % move to output directory
    a = dir('*.mat');               % N x 1 structure of names of .mat files in current directory
    cd(rootfolder)                           % move back
    s = sprintf('%s',a.name);       % concatonate all names into one string
    b = count(s,[fPath,fName]);     % count occurences of PrevName
    SaveName = fullfile(fPath, [fName,'_', sprintf('%d', b + 1), fExt]);  % append number to filename
else
    SaveName = PrevName;
end