%% Global simulation configuration

nslice = 2;
nsubstrate = 1;
T = 4*60*60;

xscale = 'um';
tscale = 's';
mscale = 'ug';
substrateName = [ "test" ];


%% Erosion experimental data
exp_dX = [
    25 75
    25 75
];
exp_t  = [
    0
    1
];

%% Mass experimental data


exp_mrel = [ 0 ]; % ug Placeholder

%% Set deltas

dxs = [ ...
    0.1 ... 1
    1   ... 2
];

dts = [...
    0.001 ... 1
    0.02 ... 2
];

dt_mass = 60; %60; % [s]


%% Set substrates
u0 = [
    1 % 1
    1 % 2
];


%% Set diffusion coefficients

D = cell(nslice, 1);
D{1}  = [ 0.5 ];
D{2}  = [ 0.5 ];

% Inward diffusivity radial asymmetric factor
use_asymmetric_diffusivity = false;
alpha = cell(nslice, 1);
alpha{1}  = [ 1 ];
alpha{2}  = [ 1 ];

lambda = [ 1 ];


%% Set decay rates
beta = cell(nslice, 1);
beta{1}  = [ 0 ];
beta{2}  = [ 0 ];


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