clear all
close all
clc

file = "../data/bws_position.csv";
T = readtable(file);
Txy = tfestimate(T.Input,T.Output,[],[],[],1000)

% compare results with:
% .\etfe.exe -d ../data/bws_position.csv -n 10000 -f 1000