function val = eval_obj_fun(I, I1, T, s, sigma, options)
% I1 is the reference image
try
    [rows,cols,B] = size(I1);
    I2 = transform(I, T, s, sigma, cols, rows, options);
    val = calc_reg_metric(I1, I2, options);
catch err
    if strcmp(err.identifier, 'Registration:ExceedBound')
        val = Inf;
    else
        throw(err);
    end
end
