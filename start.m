close all force;
clear all;

matlabrc;

ones(10)*ones(10);

warning ('off','all');

dirmat=[pwd,'/../'];
addpath([dirmat,'m_map'])
disp(['AMEDA path = ',pwd])
addpath([pwd])
addpath([pwd,'/tools'])
addpath([pwd,'/sources'])
addpath([pwd,'/plot'])

