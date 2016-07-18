function conv = test_convergence(errors, convergence_t)
conv = 0;
if length(errors) <= 10
    return;
end

iter = length(errors) - 1;
obj_mag = abs(errors(10));
prev_mean = mean(errors(end-9:end-5));
curr_mean = mean(errors(end-4:end));
if prev_mean - curr_mean <= convergence_t * obj_mag
    disp(['Finishes at iteration ',num2str(iter), '. The objective ',...
        'function has value ',num2str(errors(end))]);
    conv = 1;
end
