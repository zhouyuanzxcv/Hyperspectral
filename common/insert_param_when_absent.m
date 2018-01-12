function options = insert_param_when_absent(options, fieldname, default_value)
if ~isfield(options, fieldname)
    options.(fieldname) = default_value;
end

end

