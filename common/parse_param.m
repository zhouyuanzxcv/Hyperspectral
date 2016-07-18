function value = parse_param(options, field_name, default_value)
%PARSE_PARAM Summary of this function goes here
%   Detailed explanation goes here
if isempty(options) || ~isstruct(options) || ~isfield(options,field_name)
    value = default_value;
else
    value = options.(field_name);
end