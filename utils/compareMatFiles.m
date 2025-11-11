function compareMatFiles(file1, file2, tol)
%COMPAREMATFILES Compare two MATLAB .mat files variable by variable
%(GENERATED CODE)
%
%   compareMatFiles(file1, file2)
%   compareMatFiles(file1, file2, tol)
%
%   Compares the content of two .mat files.
%   Reports identical, differing, missing, or extra variables.
%   Numeric arrays are compared within a tolerance (default: 1e-12).
%
%   Example:
%       compareMatFiles('L1_coarse.mat', 'L1_coarse_.mat', 1e-10)

    if nargin < 3
        tol = 1e-12;
    end

    % Load both .mat files into structs
    A = load(file1);
    B = load(file2);

    varsA = fieldnames(A);
    varsB = fieldnames(B);

    fprintf('🔍 Comparing "%s" and "%s"\n\n', file1, file2);

    % Variables present only in one of the files
    onlyA = setdiff(varsA, varsB);
    onlyB = setdiff(varsB, varsA);

    if ~isempty(onlyA)
        fprintf('⚠️  Variables only in %s:\n', file1);
        disp(onlyA);
    end
    if ~isempty(onlyB)
        fprintf('⚠️  Variables only in %s:\n', file2);
        disp(onlyB);
    end

    % Compare common variables
    common = intersect(varsA, varsB);
    if isempty(common)
        fprintf('No common variables to compare.\n');
        return;
    end

    fprintf('Comparing %d common variables...\n\n', numel(common));

    for i = 1:numel(common)
        name = common{i};
        valA = A.(name);
        valB = B.(name);

        % Check class and size first
        if ~strcmp(class(valA), class(valB))
            fprintf('❌ %s: different types (%s vs %s)\n', ...
                name, class(valA), class(valB));
            continue;
        elseif ~isequal(size(valA), size(valB))
            fprintf('❌ %s: different sizes (%s vs %s)\n', ...
                name, mat2str(size(valA)), mat2str(size(valB)));
            continue;
        end

        % Compare values
        if isnumeric(valA)
            diffNorm = norm(valA(:) - valB(:), inf);
            if diffNorm < tol
                fprintf('✅ %s: numerically identical (max diff = %.3e)\n', ...
                    name, diffNorm);
            else
                fprintf('⚠️  %s: numeric difference (max diff = %.3e)\n', ...
                    name, diffNorm);
            end
        elseif isequaln(valA, valB)
            fprintf('✅ %s: identical (non-numeric)\n', name);
        else
            fprintf('⚠️  %s: different values (non-numeric)\n', name);
        end
    end

    fprintf('\n✅ Comparison complete.\n');
end