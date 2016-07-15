classdef ParameterAnalysis
    %PARAMETERANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        StepSizeMultiplier 
    end
    
    methods
        function obj = ParameterAnalysis()
            obj.StepSizeMultiplier = 10;
        end
        
        function options = autoParamSelection(obj, fcn_run, options, ...
                list_params, step_size, range)   
            disp('Start searching best parameters. This could take a while.');
%             if 0
            curr_val = fcn_run(options);
            list_search_options = options;
            list_search_value = curr_val;
            while 1
                nbhds = obj.findNeighboringOptions(options, list_params, ...
                    step_size, range);
                nbhds_vals = zeros(1,length(nbhds));
                for i = 1:length(nbhds)
                    [val, list_search_options, list_search_value] = ...
                        obj.run(fcn_run, nbhds(i), list_search_options, ...
                        list_search_value);
                    nbhds_vals(i) = val;
                end
                
                if obj.testLocalMinimum(curr_val, nbhds_vals)
                    break;
                else
                    [min_val,ind] = min(nbhds_vals);
                    options = nbhds(ind);
                    curr_val = min_val;
                end
            end
            save('param_analysis.mat','list_search_options',...
                'list_search_value','options');
%             end
            obj.printListSearch(list_search_options, list_search_value, ...
                list_params);
            disp('The best searched parameters are');
            options
        end
        
        function printListSearch(obj, list_search_options, list_search_value,...
                list_params)
            disp('--------------------- Search Results -----------------');
            list_params = [list_params,{'value'}];
            disp(list_params);
            
            N = length(list_search_options);
            M = length(list_params);
            X = zeros(N, M);
            for i = 1:N
                options = list_search_options(i);
                value = list_search_value(i);
                for j = 1:M-1
                    name = list_params{j};
                    X(i,j) = options.(name);
                end
                X(i,M) = value;
            end
            disp(X);
        end
        
        function [val, list_search_options, list_search_value] = run(obj, ...
                fcn_run, options, list_search_options, list_search_value)
            find = 0;
            for i = 1:length(list_search_options)
                if isequal(options, list_search_options(i))
                    find = 1;
                    val = list_search_value(i);
                    break;
                end
            end
            if ~find
                val = fcn_run(options);
                list_search_options = [list_search_options, options];
                list_search_value = [list_search_value, val];
            end
        end
        
        function T = testLocalMinimum(obj, curr_val, nbhds_vals)
            T = all(curr_val <= nbhds_vals);
        end
                
        function nbhds = findNeighboringOptions(obj, options, list_params, ...
                step_size, range)
            nbhds = [];
            options_start = options;
            for i = 1:length(list_params)
                options = options_start;
                name = list_params{i};
                forw_val = obj.findNeighboringValue(options_start.(name),...
                    step_size, i, 1);
                back_val = obj.findNeighboringValue(options_start.(name),...
                    step_size, i, 0);
                upper_bd = range(2,i);
                lower_bd = range(1,i);
                
                options.(name) = min(forw_val, upper_bd);
                nbhds = [nbhds, options];
                options.(name) = max(back_val, lower_bd);
                nbhds = [nbhds, options];
            end
        end
        
        function new_val = findNeighboringValue(obj, curr_val, ...
                step_size, i, direction)
            order_of_magnitude = @(x) floor(log10(x));
            
            if isnumeric(step_size)
                size = step_size(i);
                new_val = obj.calcNewValue(direction, curr_val, size);
            elseif isequal(step_size,'medium')
                possible_vals = [1,2,5];
                size = 2;
                new_val = obj.calcNewValue(direction, curr_val, size);
                order = order_of_magnitude(new_val);
                new_val = new_val/(10^order);
                [~,ind] = min(abs(new_val - possible_vals));
                new_val = possible_vals(ind);
                new_val = new_val * (10^order);
            end        
        end
        
        function new_val = calcNewValue(obj, direction, curr_val, size)
            if direction == 1
                new_val = curr_val * size; 
            else
                new_val = curr_val / size; 
            end
        end
    end
end

