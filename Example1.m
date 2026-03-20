% Example1.m

% Environment
V = 5000000; % m3, Volume

% Sources
P = 0.015; % mg/h - rate of production of X
F_in = 0.02; % mg/h - inflow of X

% Sinks
F_out = 0.047; % mg/h - rate of outflow of X
L = 0.008; % mg/h - rate of loss of X

% Change in X with time
dc_dt = 0.001; % mg/m3-h

dm_dt = dc_dt * V; % mg/m3-h * (m3) = mg/h

E = dm_dt + L + F_out - P - F_in; % mg/h


