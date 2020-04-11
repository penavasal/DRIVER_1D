function plot_driver(varargin)

%clear all
clearvars -except varargin
format long
restoredefaultpath
warning('off', 'MATLAB:nearlySingularMatrix')

% Main file
% Version Dic 2019

VERSION='Code';
RUTA_unix='/';
RUTA_win ='\';
 
% PATH to files
if ismac || isunix  % Code to run on Mac or Unix plaform 
    s=strcat(pwd,RUTA_unix,VERSION);
    s1=strcat(s,'/Model');
elseif ispc         % Code to run on Windows platform
    s=strcat(pwd,RUTA_win,VERSION);
    s1=strcat(s,'\Model');
else
    disp('Platform not supported')
    stop
end 

addpath(path,s);
addpath(path,s1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE PARAMETERS, GEOMETRY, MATERIAL and BOUNDARY CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0
    str='FILE.mat';
else
    str=strcat(varargin{1},'.mat');
end
plot_constitutive(1,str);




