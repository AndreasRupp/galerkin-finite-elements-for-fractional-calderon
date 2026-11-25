clear all; clc; close all;
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'src'));
profile='standard'; % This profile uses the meshes and noise levels of the paper.
results=cell(3,1);
for example=1:3
    results{example}=fcp_run_experiment(example,profile);
end
