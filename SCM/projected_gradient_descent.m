function X = projected_gradient_descent(X, fcn_obj, fcn_der, fcn_project, delta_t0)
X2_old = X;
E0 = fcn_obj(X);
der = fcn_der(X);

delta_t = delta_t0;

while 1
    X2_new = X - delta_t*der;
    X2_new = fcn_project(X2_new);

    E0 = [E0; fcn_obj(X2_new)];
    if E0(end) >= E0(end-1)
        X = X2_old;
%         disp(['Tried times of delta_t is ',num2str(log10(delta_t/delta_t0))]);
        break;
    else
        X2_old = X2_new;
        delta_t = delta_t*10;
    end
end    
