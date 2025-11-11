%% Hyperparameters for sensitivity analysis
%                        % --           -           ok         +         ++
betaval   = 0.0;         % Na           Na          0.0;       0.00001;  0.0001
Dcoreval  = 0.0000006;   % 0.000000006; 0.00000006; 0.0000006; 0.000006; 0.00006
Dcoat1val = 0.000005;    % 0.00000005;  0.0000005;  0.000005;  0.00005;  0.0005
Dcoat2val = 0.000001;    % 0.00000001;  0.0000001;  0.000001;  0.00001;  0.0001
acoreval  = 0.5;         % 
acoat1val = 0.2;         % 
acoat2val = 1.0;         % 
lambdaval = 0.05;        % 0.0 (N0);    0.01;       0.05;      0.25;     0.50;    1.0; 1000;



%% Global simulation configuration

nslice = 11;
nsubstrate = 2;
T = 4*60*60;

xscale = 'um';
tscale = 's';
mscale = 'ug';
substrateName = ["probiotics", "antioxydant"];


%% Erosion experimental data
exp_dX = [
    280  5.0000  0.6500  0.0180  0.0180  0.0180  0.0180  0.0180  0.0180  0.0180  0.0190
    280  5.0000  0.6500  0.0180  0.0180  0.0180  0.0180  0.0180  0.0180  0.0180  0.0175
    280  5.0000  0.6500  0.0180  0.0180  0.0180  0.0180  0.0180  0.0180  0.0180  0.0159
    280  5.0000  0.6500  0.0180  0.0180  0.0180  0.0180  0.0180  0.0180  0.0180  0.0149
    280  5.0000  0.6500  0.0180  0.0180  0.0180  0.0180  0.0180  0.0180  0.0180  0.0141
    280  5.0000  0.6500  0.0180  0.0180  0.0180  0.0180  0.0180  0.0180  0.0135  0
    280  5.0000  0.6500  0.0180  0.0180  0.0180  0.0180  0.0180  0.0108  0       0
    280  5.0000  0.6500  0.0180  0.0180  0.0180  0.0180  0.0104  0       0       0
    280  5.0000  0.6500  0.0180  0.0180  0.0180  0.0140  0       0       0       0
    280  5.0000  0.6500  0.0180  0.0180  0.0154  0       0       0       0       0
    280  5.0000  0.6500  0.0180  0.0180  0.0018  0       0       0       0       0
    280  5.0000  0.6500  0.0180  0.0119  0       0       0       0       0       0
    280  5.0000  0.6500  0.0180  0.0066  0       0       0       0       0       0
];
exp_t  = [
    0
    1800
    3600
    5400
    7200
    8100
    9000
    9900
    10800
    11700
    12600
    13500
    14400
];

%% Mass experimental data

exp_ul = 100; % microliters of solution

% Average initial internal concentration of TA [ug]
exp_cmax = 4.7210;

% Released mass concentration at exp_t
exp_crel = [
    0
    0.210666666666667
    0.410
    0.614
    0.901666666666667
    1.22966666666667
    1.50866666666667
    1.792
    2.009
    2.31533333333333
    2.71633333333333
    3.018
    3.25066666666667
];

exp_mrel = exp_crel * exp_ul; % ug

%% Set deltas

dxs = [ ...
   35.000 ... Core 1
    0.500 ... Core 2
    0.005 ... Core 3
    0.001 ... Coating 1
    0.001 ... Coating 2
    0.001 ... Coating 3
    0.003 ... Coating 4
    0.003 ... Coating 5
    0.003 ... Coating 6
    0.003 ... Coating 7
    0.001 ... Coating 8
];

dts = [...
    1.00 ... Core 1
    0.05 ... Core 2
    0.01 ... Core 3
    0.01 ... Coating 1
    0.01 ... Coating 2
    0.01 ... Coating 3
    0.01 ... Coating 4
    0.01 ... Coating 5
    0.01 ... Coating 6
    0.01 ... Coating 7
    0.01 ... Coating 8
];

dt_mass = 300; %60; % [s]


%% Set substrates
core_layers = 3;
exp_mtot = exp_cmax * exp_ul; % ug
r0 = sum(exp_dX(1, 1:core_layers));
r1 = sum(exp_dX(1, :));
exp_V = 4/3*pi*(r1^3-r0^3);
exp_c0 = exp_mtot/exp_V;

slice_edges = [0 cumsum(exp_dX(1,:))];


Vc1_2 = 4/3*pi * (slice_edges(core_layers+3)^3 - slice_edges(core_layers+1)^3);
Vc3_8 = 4/3*pi * (slice_edges(end)^3 - slice_edges(core_layers+3)^3);
cc3_8 = exp_c0 * (Vc1_2 + Vc3_8) / (2 * Vc1_2 + Vc3_8);
cc1_2 = 2 * cc3_8;
u0 = [
    1 0        % Core 1
    1 0        % Core 2
    1 0        % Core 3
    0 cc1_2   % Coating 1
    0 cc1_2   % Coating 2
    0 cc3_8   % Coating 3
    0 cc3_8   % Coating 4
    0 cc3_8   % Coating 5
    0 cc3_8   % Coating 6
    0 cc3_8   % Coating 7
    0 cc3_8   % Coating 8
];


%% Set diffusion coefficients

D = cell(nslice, 1);
D{1}  = [0.1       Dcoreval];  % Core 1
D{2}  = [0.0000001 Dcoreval];  % Core 2
D{3}  = [0.0000001 Dcoreval];  % Core 3
D{4}  = [0.0000001 Dcoat1val]; % Coating 1
D{5}  = [0.0000001 Dcoat1val]; % Coating 2
D{6}  = [0.0000001 Dcoat2val];  % Coating 3
D{7}  = [0.0000001 Dcoat2val];  % Coating 4
D{8}  = [0.0000001 Dcoat2val];  % Coating 5
D{9}  = [0.0000001 Dcoat2val];  % Coating 6
D{10} = [0.0000001 Dcoat2val];  % Coating 7
D{11} = [0.0000001 Dcoat2val];  % Coating 8

% Inward diffusivity radial asymmetric factor
use_asymmetric_diffusivity = true;
alpha = cell(nslice, 1);
alpha{1}  = [1, acoreval];
alpha{2}  = [1, acoreval];
alpha{3}  = [1, acoreval];
alpha{4}  = [1, acoat1val];
alpha{5}  = [1, acoat1val];
alpha{6}  = [1, acoat2val];
alpha{7}  = [1, acoat2val];
alpha{8}  = [1, acoat2val];
alpha{9}  = [1, acoat2val];
alpha{10} = [1, acoat2val];
alpha{11} = [1, acoat2val];

lambda = [100, lambdaval];


%% Set decay rates
beta = cell(nslice, 1);
beta{1}  = [0.000001 betaval];
beta{2}  = [0.000001 betaval];
beta{3}  = [0.000001 betaval];
beta{4}  = [0.000001 betaval];
beta{5}  = [0.000001 betaval];
beta{6}  = [0.000001 betaval];
beta{7}  = [0.000001 betaval];
beta{8}  = [0.000001 betaval];
beta{9}  = [0.000001 betaval];
beta{10} = [0.000001 betaval];
beta{11} = [0.000001 betaval];


%% Create output struct

P = struct();
P.nslice = nslice;
P.nsubstrate = nsubstrate;
P.T = T;

P.xscale = xscale;
P.tscale = tscale;
P.mscale = mscale;
P.substrateName = substrateName;

P.exp_dX = exp_dX;
P.exp_t = exp_t;

P.dxs = dxs;
P.dts = dts;
P.dt_mass = dt_mass;

P.exp_mrel = exp_mrel;
P.u0 = u0;

P.D = D;
P.use_asymmetric_diffusivity = use_asymmetric_diffusivity;
P.alpha = alpha;

P.lambda = lambda;
P.beta = beta;