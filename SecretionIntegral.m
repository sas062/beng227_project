close all; clear all; clc;


%% Calculate sG
sG_peak = 30; % amol/min/islet
sG_peak = sG_peak/60; % amol/s/islet
N = 1000; % ~number of cells in an islet
Na = N*0.15; % Number of alpha-cells in the islet
sG_peak_cell = sG_peak/Na; % amol/s/cell


%% Calculate scaling factor
r = 5; % um
V1 = (4/3*pi*(r+0.05)^3-4/3*pi*r^3)*(1e-15); % Volume of realistic extracellular space
V2 = 5*5*10*4*(1e-15); % Volume of model extracellular space
C_per_t1 = sG_peak_cell/V1*1e-18*1e6; % dCdt for realistic volume
C_per_t2 = sG_peak_cell/V2*1e-18*1e6; % dCdt for model volume
ratio = C_per_t1/C_per_t2;