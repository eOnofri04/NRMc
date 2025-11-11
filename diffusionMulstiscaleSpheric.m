%% Initial variable and parameter loading
% Examples include, eg
%  exp_config_fn = "test_real/exp_2412"; do_plot = true; do_use_exp = true; do_save_U = true;
%  exp_config_fn = "test_academic/L1_coarse"; do_plot = true;

if ~exist('do_plot', 'var')
    do_plot = false;
    do_mov = false;
end
if ~exist('do_mov', 'var')
    do_mov = false;
elseif ~do_plot
    do_mov = false;
end
if ~exist('do_pause', 'var')
    do_pause = false;
end
if ~exist('do_plot_finer', 'var')
    do_plot_finer = false;
end
if ~exist('do_export', 'var')
    do_export = false;
end
if ~exist('do_save_U', 'var')
    do_save_U = false;
end
if ~exist('do_use_exp', 'var')
    do_use_exp = false;
end
if ~exist('exp_config_fn', 'var')
    exp_config_fn = "test_real/exp_2412";
end

P = loadParams(exp_config_fn);

%% Load simulation configuration
nslice = getParam(P, 'nslice');
nsubstrate = getParam(P, 'nsubstrate');
T = getParam(P, 'T');

xscale = getParam(P, 'xscale');
tscale = getParam(P, 'tscale');
mscale = getParam(P, 'mscale');


%% Load experimental dataset
exp_dX = getParam(P, 'exp_dX');
exp_t  = getParam(P, 'exp_t');
assert(isequal(size(exp_dX), [length(exp_t), nslice]));

if do_use_exp
    exp_mrel = getParam(P, 'exp_mrel');
end

%% Load deltas, determine ordering, and evaluate time discretisation
dxs = getParam(P, 'dxs');
assert(length(dxs) == nslice);

dx_min = min(dxs);
dx_ratios = dxs ./ dx_min;

dts = getParam(P, 'dts');
assert(length(dts) == nslice);

dt_min = min(dts);
for i = 1 : nslice
    assert( mod(dts(i), dt_min) == 0 )
end
dt_ratios = round(dts ./ dt_min);

[~, slice_order] = sort(dts);

dt_mass = getParam(P, 'dt_mass');
mass_steps = dt_mass / dt_min;

Ns = T ./ dts;
N  = T / dt_min;
Nm = T / dt_mass;

%% Create space discretization
dX = exp_dX(1,:);
assert(length(dX) == nslice);
slice_edges = [0 cumsum(dX)];

gridsB     = cell(nslice,1); % Radii of the cell centers (baricentric)
gridsN     = cell(nslice,1); % Radii of the cell edges (nodal)
gridsBsq   = cell(nslice,1);
gridsNsq   = cell(nslice,1);
V_B        = cell(nslice,1); % Volumes of the cells (baricentric)
V_BB_left  = cell(nslice,1); % Volumes of left-boundary cells
V_BB_right = cell(nslice,1); % Volumes of right-boundary cells
Ms = zeros(nslice,1);
for i = 1 : nslice
    assert( mod(dX(i), dxs(i)) == 0 )
    gridsB{i} = (slice_edges(i)+dxs(i)/2 : dxs(i) : slice_edges(i+1)-dxs(i)/2)';
    gridsBsq{i} = gridsB{i}.^2;
    Ms(i) = length(gridsB{i});
    gridsN{i} = (slice_edges(i) : dxs(i) : slice_edges(i+1))';
    assert( length(gridsN{i}) == Ms(i) + 1 )
    gridsNsq{i} = gridsN{i}.^2;
    V_B{i} = 4/3*pi * (gridsN{i}(2:end).^3 - gridsN{i}(1:end-1).^3);
    V_BB_left{i} = 4/3*pi * (gridsN{i}(1).^3 - (gridsN{i}(1)-dxs(i)).^3);
    V_BB_right{i} = 4/3*pi * ((gridsN{i}(end)+dxs(i)).^3 - gridsN{i}(end).^3);
end
clear("i");


%% Load substrates configuration
substrateName = getParam(P, 'substrateName');
u0 = getParam(P, 'u0');
assert( isequal(size(u0), [nslice, nsubstrate]) );


%% Load diffusion and evaluate constants

D = getParam(P, 'D');
assert( isequal(size(D), [nslice, 1]) );
for i = 1 : nslice
    assert( isequal(size(D{i}), [1, nsubstrate]) );
end
clear("i");

% Evaluate the constant that evals the mass transfered through interface i.
% Must be applied to u(i)-u(i-1)
% dt/(rb^2 dr) D rn^2/dr Vb
cnstN = cell(nslice, 1);
for i = 1 : nslice
    r_tmp = gridsB{i}(end)+dxs(i);
    cnstN{i} = (gridsN{i}.^2) ./ ([gridsB{i}; r_tmp].^2) ...
        .* dts(i) ./ (dxs(i) .^ 2) .* D{i} ...
        .* ([V_B{i}; V_BB_right{i}]);
    if i > 1
        cnstN{i}(1,:) = cnstN{i}(1,:) .* D{i-1} * 2 ./ (D{i} + D{i-1});
    end
    if i < nslice
        cnstN{i}(end,:) = cnstN{i}(end,:) .* D{i+1} * 2 ./ (D{i} + D{i+1});
    end
    clear("r_tmp", "V_tmp");
end
clear("i");

%% Load and process diffusivity radial asymmetric factor reduction factor
use_asymmetric_diffusivity = getParam(P, 'use_asymmetric_diffusivity');
if use_asymmetric_diffusivity
    alpha = getParam(P, 'alpha');
    assert( isequal(size(alpha), [nslice, 1]) );
    for i = 1 : nslice
        assert( isequal(size(alpha{i}), [1, nsubstrate]) );
    end
    clear("i");

    % Check whether it is really needed
    use_asymmetric_diffusivity = false;
    for slice_idx = 1 : nslice
        if any(alpha{slice_idx} ~= 1)
            use_asymmetric_diffusivity = true;
        end
    end
end

%% Load Robin constants for outer layer

lambda = getParam(P, 'lambda');
assert(length(lambda) == nsubstrate)

%% Load decay rates
beta = getParam(P, 'beta');

assert( isequal(size(beta), [nslice, 1]) );
for i = 1 : nslice
    assert( isequal(size(beta{i}), [1, nsubstrate]) );
end
clear("i");

%% Check CFL
for i = 1 : nslice
    assert((dts(i) / dxs(i)^2 * max(D{i})) < .5);
    display( "CFL: " + dts(i) / dxs(i)^2 * max(D{i}));
end
clear("i");

%% Initialise solution vectors
U = cell(nslice, 1);
M_internal = zeros(Nm+1, nsubstrate);
M_released = zeros(Nm+1, nsubstrate);
M_decayed  = zeros(Nm+1, nsubstrate);
for i = 1 : nslice
    U{i}          = ones(Ms(i), nsubstrate) .* u0(i,:);
    % eval initial mass
    M_internal(1, :) = M_internal(1, :) + V_B{i}' * U{i};
end
clear("i");
if do_save_U
    U_all = cell(nslice, Nm+1);
    for i = 1 : nslice
        U_all{i, 1} = U{i};
    end
end

%% Setup boundary conditions
b_left_f  = cell(nslice, 1);
b_right_f = cell(nslice, 1);
update_left  = false(nslice, 1);
update_right = false(nslice, 1);

% left boundary: Symmetrical condition due to spheric symmetry
b_left_f{1} = @(U, t) (U{1}(1,:));
% right boundary: Substrates go immetiatly to infinite
b_right_f{end} = @(U, t) (zeros(1, nsubstrate));
% b_right_f{end} = @(U, t) (lambda .* U{end}(end,:)); CANNOT DO THIS AS
% IT IS LATER RELEASED!!!
update_right(end) = true;

% inner boundaries:
%  the smaller dt reads from the bigger dt
%  the bigger dt will be updated on the smaller updates, hence symmetric
%  in case of a tie (should not happen) right is considered smaller
for i = 2 : nslice
    if dts(i) <= dts(i-1)
        b_left_f{i} = @(U, t) (U{i-1}(end,:));
        update_left(i) = true;
    else
        b_left_f{i} = @(U, t) (U{i}(1,:));
    end
end
for i = 1 : nslice-1
    if dts(i) < dts(i+1)
        b_right_f{i} = @(U, t) (U{i+1}(1,:));
        update_right(i) = true;
    else
        b_right_f{i} = @(U, t) (U{i}(end,:));
    end
end
clear("i");

%% Plot related variables and functions
axs = cell(nsubstrate, nslice+1);
scs = cell(nsubstrate, nslice+1);
scs_mi = cell(nsubstrate, 1);
scs_mr = cell(nsubstrate, 1);
scs_md = cell(nsubstrate, 1);
bgs = cell(nsubstrate, nslice);
if do_plot
    fig = figure(Position=[0 0 3840 2160]);% 1920 988
    if do_mov
        % Set up video writer
        videoHandle = VideoWriter('myAnimation.mp4', 'MPEG-4');
        videoHandle.FrameRate = 2;  % frames per second
        open(videoHandle);
    end
    for i = 1 : nsubstrate
        % mass panel is in proportion 3:1 with slices (x-axis)
        for j = 1 : nslice
            axs{i, j} = subplot(nsubstrate, nslice+3, (i-1)*(nslice+3)+j);
        end
        axs{i, nslice+1} = subplot(nsubstrate, nslice+3, ...
            (i-1)*(nslice+3)+(nslice+1:nslice+3) ...
        );
    end
    for j = 1 : nslice
        for i = 1 : nsubstrate
            bgs{i, j} = fill(axs{i, j}, ...
                [slice_edges(j) slice_edges(j) slice_edges(j+1) slice_edges(j+1)], ...
                1.1*[0, max(u0(:, i)), max(u0(:, i)), 0], ...
                [0.9290 0.6940 0.1250], ...
                FaceAlpha=.3 ...
            );
            hold(axs{i, j}, "on");
            set(axs{i, j}, "Layer", "top");
            scs{i, j} = bar(axs{i, j}, gridsB{j}, U{j}(:,i), .9);
            xlim(axs{i, j}, [slice_edges(j), slice_edges(j+1)]);
            ylim(axs{i, j}, [0, 1.1*max(u0(:, i))]);
            hold(axs{i, j}, "off");
            grid(axs{i, j}, "on")
        end
        title(axs{1, j}, "slice " + string(j));
        xlabel(axs{end, j}, "space [" + xscale + "]");
    end

    for i = 1 : nsubstrate
        ylabel(axs{i, 1}, substrateName(i));
        
        scs_mi{i} = scatter(axs{i, end}, ...
            0:dt_mass:T, [M_internal(1,i); nan(Nm, 1)]', '.' ...
        );
        hold(axs{i, end}, 'on');
        scs_mr{i} = scatter(axs{i, end}, ...
            0:dt_mass:T, [0; nan(Nm, 1)]', '.' ...
        );
        scs_md{i} = scatter(axs{i, end}, ...
            0:dt_mass:T, [0; nan(Nm, 1)]', '.' ...
        );
        if do_use_exp
            scatter(axs{i, end}, ...
                exp_t, exp_mrel', 'xr' ...
            );
        end
        hold(axs{i, end}, 'off');
        if i < nsubstrate
            lgnd = legend(axs{i, end}, ...
                [scs_mi{i}, scs_mr{i}, scs_md{i}], ...
                "Internal mass", "Released mass", "Decayed mass", ...
                Location='north', ...
                Orientation='horizontal' ...
            );
            lgnd.Position(2)=(nsubstrate-1)/nsubstrate-.01;
        end
    end
    title(axs{1, end}, "mass");
    xlabel(axs{end, end}, "time [" + tscale + "]");
    sgtitle(fig, 'Time = 00:00:00');

    if do_pause
        pause()
    end
    if do_mov
        frame = getframe(fig);
        writeVideo(videoHandle, frame);
    end
    if do_export
        iter = 0;
        saveas(fig, "iter_"+string(iter)+".png");
        iter = iter+1;
    end
    
end

%% Initialise values for the main loop
% U_left{slice_idx} and U_right{slice_idx} hold the boundary conditions for
% the `slice_idx` slice.
U_left   = arrayfun(@(i) b_left_f{i}(U, 1), 1:nslice, 'UniformOutput', false);
U_right  = arrayfun(@(i) b_right_f{i}(U, 1), 1:nslice, 'UniformOutput', false);
dM_left  = arrayfun(@(i) zeros(1, nsubstrate), 1:nslice, 'UniformOutput', false);
dM_right = arrayfun(@(i) zeros(1, nsubstrate), 1:nslice, 'UniformOutput', false);
prec_last_slice = nslice;
%% Main loop
wbar = waitbar(0, 'Starting');
for t_idx = 2 : N
    if mod(t_idx, 100) == 0
        waitbar(t_idx/N, wbar, sprintf('Progress: %d %%', floor(t_idx/N*100)));
    end
    % t_idx goes according dt_min
    % currtime states the effective time
    currtime = t_idx * dt_min;

    % t_idxs are the time indices for each single slice according the
    % following pattern:
    %   1 -> dt_r = 1
    %   dt_r+1 -> 2*dt_r = 2
    %   ...
    t_idxs = ceil(t_idx./dt_ratios);
    
    % tm_idx is the time index for the mass evaluation.
    % It follows analogous pattern to t_idxs but it is increased by one to
    % acknowledge time 0, ie.
    %   1 -> mass_steps = 2
    %   mass_steps+1 -> 2*mass_steps = 3
    %   ...
    % Released and decayed mass are stored incrementally on tm_idx.
    % On the last t_idx before updating, ie at mass_steps, 2*mass_steps ...
    % the total internal mass is evaluated before moving to the next t_idx.
    tm_idx = ceil(t_idx / mass_steps) + 1;

    % eval current outboundary
    % last_slice evaluates up to wich slice we should evaluate.
    % if a new last_slice is obtained, the prec_last_slice is flushed.
    % b_c is the index of the nodal point to the left w.r.t. r_curr.
    % Note. if b_c equals a nodal point, then that nodal point is chosen
    % OLD. r_curr = evalRcurr(currtime, Xtot);
    r_curr = evalRcurr(currtime, exp_t, exp_dX);
    last_slice = find(slice_edges > r_curr, 1, 'first') - 1;
    if isempty(last_slice)
        last_slice = nslice;
    end
    if prec_last_slice ~= last_slice
        assert(last_slice == prec_last_slice - 1)
        M_released(tm_idx, :) = M_released(tm_idx, :) + ...
            V_B{prec_last_slice}' * U{prec_last_slice};
        U{prec_last_slice}(:) = 0;
        update_right(last_slice) = true;
        b_right_f{last_slice} = @(U, t) (zeros(1, nsubstrate));
        cnstN{last_slice}(end, :) = (gridsN{last_slice}(end)^2) ...
            ./ ((gridsB{last_slice}(end)+dxs(last_slice))^2) ...
            .* dts(last_slice) ./ (dxs(last_slice)^2) ...
            .* D{last_slice} ...
            .* (V_BB_right{last_slice});
        if do_plot
            for sub_idx = 1 : nsubstrate
                bgs{sub_idx, prec_last_slice}.XData(3:4) = ...
                    slice_edges(prec_last_slice+1);
                bgs{sub_idx, prec_last_slice}.FaceColor = [0 0 0];
            end
        end
        prec_last_slice = last_slice;
    end
    b_c = find(gridsN{last_slice} >= r_curr, 1, 'first');
    
    %% Update each slice in dt order
    % It happens when t_idxs(slice_idx) = dt_r, 2*dt_r, ...
    for slice_idx_ord = 1 : nslice
        slice_idx = slice_order(slice_idx_ord);
        if slice_idx > last_slice
            continue;
        end
        if t_idxs(slice_idx) * dt_ratios(slice_idx) == t_idx

            % Eval the mass transfered at interface i, i.e. i -> i+1
            dM_N = cnstN{slice_idx} .* ( ...
                [U_left{slice_idx}; U{slice_idx}] ...
                - [U{slice_idx}; U_right{slice_idx}] ...
            );
            % diffusivity correction (negative elements are scaled)
            if logical(use_asymmetric_diffusivity) % logical prevents warning
                for sub_idx = 1 : nsubstrate
                    dM_N(dM_N(:, sub_idx) < 0, sub_idx) = ...
                        dM_N(dM_N(:, sub_idx) < 0, sub_idx) ...
                        .* alpha{slice_idx}(sub_idx);
                end
            end
            % Correct the mass movement on the last slice if different
            % conditions are chosen (eg. robin)
            if slice_idx == last_slice % && b_c <= Ms(last_slice)
                if b_c > 1
                    dM_N(b_c, :) = cnstN{slice_idx}(b_c) .* ( ...
                        U{slice_idx}(b_c-1, :) .* lambda ...
                    );
                else
                    dM_N(b_c, :) = cnstN{slice_idx}(b_c) .* ( ...
                        U_left{slice_idx} .* lambda ...
                    );
                end
            end

            % Check that boundary dM = 0 if it should not be updated
            assert( update_left(slice_idx) || all(dM_N(1, :) == 0) );
            assert( update_right(slice_idx) || all(dM_N(end, :) == 0) );

            % update the total deltas of the current slice
            % (from the last update of neighbour slices)
            % if update_left(slice_idx)
                dM_left{slice_idx} = dM_left{slice_idx} + dM_N(1, :);
            % end
            % if update_right(slice_idx)
                dM_right{slice_idx} = dM_right{slice_idx} + dM_N(end, :);
            % end

            % Update U{slice_idx}
            U{slice_idx} = U{slice_idx} - ( ...
                dM_N(2:end, :) - dM_N(1:end-1, :) ...
                ) ./ V_B{slice_idx} ...
            ;

            % Update the slice boundary conditions
            U_left{slice_idx} = U_left{slice_idx} ...
                + dM_N(1, :) / V_BB_left{slice_idx};
            U_right{slice_idx} = U_right{slice_idx} ...
                - dM_N(end, :) / V_BB_right{slice_idx};
            
            % Move the mass from the adjacent slices through the bounsaries
            % if it is required and set their total dM to 0
            if slice_idx > 1 && update_right(slice_idx-1)
                U{slice_idx}(1, :) = ...
                    U{slice_idx}(1, :) ...
                    + dM_right{slice_idx-1} ./ V_B{slice_idx}(1);
                dM_right{slice_idx-1} = zeros(1, nsubstrate);
                % dM_right{slice_idx-1}(:) = 0;
            end
            if slice_idx < last_slice && update_left(slice_idx+1)
                U{slice_idx}(end, :) = ...
                    U{slice_idx}(end, :) ...
                    - dM_left{slice_idx+1} ./ V_B{slice_idx}(end);
                dM_left{slice_idx+1} = zeros(1, nsubstrate);
                % dM_left{slice_idx+1}(:) = 0;
            end
            % If this is the last slice, move the mass outside of the
            % simulation
            if slice_idx == last_slice
                M_released(tm_idx, :) = M_released(tm_idx, :) ...
                    + dM_right{slice_idx} ...
                    + V_B{slice_idx}(b_c:end)' * U{slice_idx}(b_c:end, :);
                dM_right{slice_idx} = zeros(1, nsubstrate);
                % dM_right{slice_idx}(:) = 0;
                U{slice_idx}(b_c:end, :) = 0;
            end

            % Decay and store the decayed
            decayed = U{slice_idx} .* beta{slice_idx} .* dts(slice_idx);
            M_decayed(tm_idx, :) = ...
                M_decayed(tm_idx, :) ...
                + V_B{slice_idx}' * decayed;
            U{slice_idx} = U{slice_idx} - decayed;

            % Update boundary conditions of this and the neighbour slices
            if slice_idx == 1
                U_left{slice_idx} = b_left_f{slice_idx}(U, currtime);
            end
            if slice_idx > 1 && update_right(slice_idx-1)
                U_left{slice_idx}    = b_left_f{slice_idx}(U, currtime);
                U_right{slice_idx-1} = b_right_f{slice_idx-1}(U, currtime);
            end
            if slice_idx < nslice && update_left(slice_idx+1)
                U_right{slice_idx}  = b_right_f{slice_idx}(U, currtime);
                U_left{slice_idx+1} = b_left_f{slice_idx+1}(U, currtime);
            end
            if slice_idx == nslice
                U_right{slice_idx} = b_right_f{slice_idx}(U, currtime);
            end

            % Update plots
            if do_plot_finer && mod(t_idx, 10) == 1
                sgtitle(fig, sprintf( ...
                    'Time = %02d:%02d:%02d', ...
                    floorDiv(currtime,3600), ...
                    floorDiv(mod(currtime, 3600), 60), ...
                    floor(mod(currtime,60)) ...
                ));
                for sub_idx = 1 : nsubstrate
                    scs{sub_idx, slice_idx}.YData = U{slice_idx}(:,sub_idx);
                    if slice_idx == last_slice
                        bgs{sub_idx, slice_idx}.XData(3:4) = r_curr;
                    end
                end
                drawnow
                if do_pause
                    pause()
                else
                    pause(.01);
                end

                if do_mov
                    frame = getframe(fig);
                    writeVideo(videoHandle, frame);
                end
            end % update plot if
        end % slice if
    end % slice for

    %% Update internal mass and plot it if required
    % it happens when t_idx = mass_steps, 2*mass_steps, ...
    if (tm_idx - 1) * mass_steps == t_idx
        for slice_idx = 1 : nslice
            if do_save_U
                U_all{slice_idx, tm_idx} = U{slice_idx};
            end
            M_internal(tm_idx, :) = M_internal(tm_idx, :) ...
                + V_B{slice_idx}' * U{slice_idx};
        end

        if do_plot
            sgtitle(fig, sprintf( ...
                'Time = %02d:%02d:%02d', ...
                floorDiv(currtime,3600), ...
                floorDiv(mod(currtime, 3600), 60), ...
                floor(mod(currtime,60)) ...
            ));
            for sub_idx = 1 : nsubstrate
                for slice_idx = 1 : nslice
                    scs{sub_idx, slice_idx}.YData = U{slice_idx}(:,sub_idx);
                end
                bgs{sub_idx, last_slice}.XData(3:4) = r_curr;
                scs_mi{sub_idx}.YData(tm_idx) = M_internal(tm_idx, sub_idx);
                scs_mr{sub_idx}.YData(tm_idx) = sum(M_released(1:tm_idx, sub_idx));
                scs_md{sub_idx}.YData(tm_idx) = sum(M_decayed(1:tm_idx, sub_idx));
            end
            drawnow();
            if do_pause()
                pause()
            end
            if do_mov
                frame = getframe(fig);
                writeVideo(videoHandle, frame);
            end
        end
    end

    if do_plot && do_export && mod(currtime, 60*60) == 0
        saveas(fig, "iter_"+string(iter)+".png");
        iter = iter+1;
    end
end
clear("wbar");
if do_mov
    close(videoHandle);
end

%% Final inner boundary sync
arrayfun(@(i) isequal(dM_left{i},  zeros(1, nsubstrate)), 1:nslice);
arrayfun(@(i) isequal(dM_right{i}, zeros(1, nsubstrate)), 1:nslice);

%% Useless variable clean

clear("t_idx", "t_idxs", "currtime", "slice_idx", "dM_N", "dM_left", "dM_right", "sub_idx");

%% Check mass conservation

Minitial = M_internal(1, :);
Mfinal = M_internal(end, :) + sum(M_released) + sum(M_decayed);
Mloss = (Minitial-Mfinal)./Minitial*100;

for i = 1 : nsubstrate
    fprintf(substrateName(i)+" [sub %d]:\n", i);
    fprintf(" - Initial mass = %f \n", Minitial(i));
    fprintf(" - Final mass   = %f \n", Mfinal(i));
    fprintf("   - Internal --> %f \n", M_internal(end,i));
    fprintf("   - External --> %f \n", sum(M_released(:,i)));
    fprintf("   - Decayed ---> %f \n", sum(M_decayed(:,i)));
    fprintf(" - Mass loss    = %e [%%] \n\n", Mloss(i));
end
clear("i");

%% Dump workspace
if do_save_U
    out_name = exp_config_fn + '_U.mat';
else
    out_name = exp_config_fn + '.mat';
end
save(out_name)







%% Define domain disgregation function

function r = evalRcurr(currtime, exp_t, exp_dX)
    t0_idx = find(exp_t <= currtime, 1, 'last');
    if t0_idx == length(exp_t)
        t0_idx = t0_idx-1;
    end
    x0 = sum(exp_dX(t0_idx,:));
    x1 = sum(exp_dX(t0_idx+1,:));
    t0 = exp_t(t0_idx);
    t1 = exp_t(t0_idx+1);
    dx = (t1-currtime) / (t1-t0) * (x0-x1);
    assert(dx >= 0);
    r = x1 + dx;
end


%% Import parameters functions

function P = loadParams(configFile)
    run(configFile);

    % Validate 'P' existence
    if ~exist('P', 'var')
        error('Config file %s must define a struct ''P''.', configFile);
    end
end

function value = getParam(P, name)
    if isfield(P, name)
        value = P.(name);
    else
        error('Parameter "%s" not found in configuration.', name);
    end
end