function options = insert_into_options(options, varargin)
for i = 1:length(varargin)/2
    options.(varargin{i*2-1}) = varargin{i*2};
end

