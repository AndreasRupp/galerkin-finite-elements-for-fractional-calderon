clear all; clc; close all;
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
result=fcp_run_experiment(1,'standard');
