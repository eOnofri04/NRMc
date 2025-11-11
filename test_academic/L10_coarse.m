%% Global simulation configuration

nslice = 10;
nsubstrate = 1;
T = 4*60*60;

xscale = 'um';
tscale = 's';
mscale = 'ug';
substrateName = [ "test" ];


%% Erosion experimental data
exp_dX = [
    10 10 10 10 10 10 10 10 10 10
    10 10 10 10 10 10 10 10 10 10
];
exp_t  = [
    0
    1
];

%% Mass experimental data


exp_mrel = [ 0 ]; % ug Placeholder

%% Set deltas

dxs = [ ...
    1 ... 1
    1 ... 2
    1 ... 3
    1 ... 4
    1 ... 5
    1 ... 6
    1 ... 7
    1 ... 8
    1 ... 9
    1 ... 10
];

dts = [...
    0.1 ... 1
    0.1 ... 2
    0.1 ... 3
    0.1 ... 4
    0.1 ... 5
    0.1 ... 6
    0.1 ... 7
    0.1 ... 8
    0.1 ... 9
    0.1 ... 10
];

dt_mass = 60; %60; % [s]


%% Set substrates
u0 = [
    1 % 1
    1 % 2
    1 % 3
    1 % 4
    1 % 5
    1 % 6
    1 % 7
    1 % 8
    1 % 9
    1 % 10
];


%% Set diffusion coefficients

D = cell(nslice, 1);
D{1}  = [ 0.5 ];
D{2}  = [ 0.5 ];
D{3}  = [ 0.5 ];
D{4}  = [ 0.5 ];
D{5}  = [ 0.5 ];
D{6}  = [ 0.5 ];
D{7}  = [ 0.5 ];
D{8}  = [ 0.5 ];
D{9}  = [ 0.5 ];
D{10} = [ 0.5 ];

% Inward diffusivity radial asymmetric factor
use_asymmetric_diffusivity = false;
alpha = cell(nslice, 1);
alpha{1}  = [ 1 ];
alpha{2}  = [ 1 ];
alpha{3}  = [ 1 ];
alpha{4}  = [ 1 ];
alpha{5}  = [ 1 ];
alpha{6}  = [ 1 ];
alpha{7}  = [ 1 ];
alpha{8}  = [ 1 ];
alpha{9}  = [ 1 ];
alpha{10} = [ 1 ];

lambda = [ 1 ];


%% Set decay rates
beta = cell(nslice, 1);
beta{1}  = [ 0 ];
beta{2}  = [ 0 ];
beta{3}  = [ 0 ];
beta{4}  = [ 0 ];
beta{5}  = [ 0 ];
beta{6}  = [ 0 ];
beta{7}  = [ 0 ];
beta{8}  = [ 0 ];
beta{9}  = [ 0 ];
beta{10} = [ 0 ];


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