function main(repo_basedir)

clearvars;close all;clc; % clear

repo_basedir = pwd;

if ~exist(fullfile(repo_basedir,'figs'))
    mkdir(fullfile(repo_basedir,'figs'))
end

% adding path for the repository
addpath(genpath(repo_basedir))

CGC_figParameters(repo_basedir) % define parameters for figures
disp('Parameter for figure generation was saved in the data folder..')
CGC_fig_Exp1 % save figures for Exp1
CGC_fig_Exp2 % save figures for Exp2
end
