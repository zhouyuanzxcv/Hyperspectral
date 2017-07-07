function conv = test_convergence(errors, convergence_t, window_size, verbose)
if nargin < 3
    window_size = 5;
end
if nargin < 4
    verbose = 1;
end

conv = 0;
if length(errors) < window_size*2
    return;
end

iter = length(errors) - 1;
obj_mag = abs(errors(window_size*2));
prev_mean = mean(errors(end - (2*window_size-1) : end - window_size));
curr_mean = mean(errors(end - (window_size-1) : end));
if prev_mean - curr_mean <= convergence_t * obj_mag
    if verbose
        disp(['Finishes at iteration ',num2str(iter), '. The objective ',...
            'function has value ',num2str(errors(end))]);
    end
    conv = 1;
end
