classdef ParameterAnalysis
    %PARAMETERANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % 1 for multiplication, 2 for addition
        StepSizeMode
    end
    
    properties
        % 1 for 1-order.
        % 2 for 2-order (cartesian product).
        % A vector to indicate size for each parameter, e.g. (5,5) means
        % 5 neighbors (including current value) for each of the 2
        % parameters.
        NeighborhoodMode
    end
    
    properties
        PrintListSearch
    end
    
    properties
        % Use 'Greedy' or 'BrutalForce' or 'BlockCoordinateDescent'
        Algorithm 
    end
    
    %% parameters used by brutal force search
    properties
        BrutalForceDepth
    end
    
    %% parameters used by BCD
    properties
        BlockCoordinateDescentDepth
    end
    
    properties
        BlockCoordinateDescentGroup
    end
    
    properties
        MaxNumberOfIterations
    end
    
    %% parameters used by all
    properties
        CheckValueInListSearchOptions
    end
    
    methods
        function obj = ParameterAnalysis()
            obj.StepSizeMode = 1;
            obj.NeighborhoodMode = 1;
            obj.PrintListSearch = 0;
            obj.Algorithm = 'Greedy';
            obj.BrutalForceDepth = 3;
            obj.CheckValueInListSearchOptions = 1;
            
            obj.MaxNumberOfIterations = 50;
            obj.BlockCoordinateDescentDepth = 5;
        end
        
        function options = autoParamSelection(obj, fcn_run, options, ...
                list_params, step_size, range)
            disp('Start searching best parameters. This could take a while.');

            switch obj.Algorithm
                case 'Greedy'
                    options = autoParamSelectionByGreedyAlgo(obj, ...
                        fcn_run, options, list_params, step_size, range);
                case 'BrutalForce'
                    options = autoParamSelectionByBrutalForce(obj, ...
                        fcn_run, options, list_params, step_size, range);
                case 'BlockCoordinateDescent'
                    options = autoParamSelectionByBCD(obj, ...
                        fcn_run, options, list_params, step_size, range);
                otherwise
            end
        end
        
        function options = autoParamSelectionByGreedyAlgo(obj, fcn_run, ...
                options, list_params, step_size, range)   
%             if 0
            curr_val = fcn_run(options);
            list_search_options = options;
            list_search_value = curr_val;
            while 1
                nbhds = obj.findNeighboringOptions(options, list_params, ...
                    step_size, range, obj.NeighborhoodMode);
                nbhds_vals = zeros(1,length(nbhds));
                for i = 1:length(nbhds)
                    [val, list_search_options, list_search_value] = ...
                        obj.run(fcn_run, list_params, nbhds(i), ...
                        list_search_options, list_search_value);
                    nbhds_vals(i) = val;
                end
                
                if obj.testLocalMinimum(curr_val, nbhds_vals)
                    break;
                else
                    [min_val,ind] = min(nbhds_vals);
                    options = nbhds(ind);
                    curr_val = min_val;
                end
                
                if obj.PrintListSearch
                    obj.printListSearch(list_search_options, list_search_value, ...
                        list_params);
                end
            end
            save('param_analysis.mat','list_search_options',...
                'list_search_value','options');
%             end

            disp('The best searched parameters are');
            options
        end
        
        function options = autoParamSelectionByBrutalForce(obj, fcn_run, ...
                options, list_params, step_size, range)
            list_search_options = [];
            list_search_value = [];
            for k = 1:obj.BrutalForceDepth
                nbhds = obj.findNeighboringOptions(options, list_params, ...
                    step_size, range, obj.NeighborhoodMode);
                nbhds_vals = zeros(1,length(nbhds));
                for i = 1:length(nbhds)
                    [val, list_search_options, list_search_value] = ...
                        obj.run(fcn_run, list_params, nbhds(i), ...
                        list_search_options, list_search_value);
                    nbhds_vals(i) = val;
                end
                
                [min_val,ind] = min(nbhds_vals);
                options = nbhds(ind);
                step_size = step_size / 2;                
            end
            save('param_analysis.mat','list_search_options',...
                'list_search_value','options');

            if obj.PrintListSearch
                obj.printListSearch(list_search_options, list_search_value, ...
                    list_params);
            end
            disp('The best searched parameters are');
            options

        end
        
        function options = autoParamSelectionByBCD(obj, ...
                fcn_run, options, list_params, step_sizes, range)
            last_options = options;
            for iter = 1:obj.MaxNumberOfIterations
                for j = 1:length(obj.BlockCoordinateDescentGroup)
                    inds = obj.BlockCoordinateDescentGroup{j};
                    list_search_options = [];
                    list_search_value = [];
                    step_size = step_sizes(inds);
                    for k = 1:obj.BlockCoordinateDescentDepth
                        nbhds = obj.findNeighboringOptions(options, list_params(inds), ...
                            step_size, range(:,inds), obj.NeighborhoodMode(inds));
                        nbhds_vals = zeros(1,length(nbhds));
                        for i = 1:length(nbhds)
                            [val, list_search_options, list_search_value] = ...
                                obj.run(fcn_run, list_params, nbhds(i), ...
                                list_search_options, list_search_value);
                            nbhds_vals(i) = val;
                        end
                        
                        [min_val,ind] = min(nbhds_vals);
                        options = nbhds(ind);
                        step_size = step_size / 2;
                    end
                end
                if ParameterAnalysis.isOptionsEqual(options, last_options, list_params)
                    disp(['BCD stopped at iteration ',num2str(iter)]);
                    break;
                else
                    last_options = options;
                end
            end
            
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
        
        function [val, list_search_options, list_search_value] = run( ...
                obj, fcn_run, list_params, options, list_search_options,...
                list_search_value)
            find = 0;
            if obj.CheckValueInListSearchOptions
                for i = 1:length(list_search_options)
                    if ParameterAnalysis.isOptionsEqual(options, ...
                            list_search_options(i), list_params)
                        find = 1;
                        val = list_search_value(i);
                        break;
                    end
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
                step_size, range, nbhd_mode)
            if length(obj.StepSizeMode) == 1
                step_mode = repmat(obj.StepSizeMode, [length(step_size),1]);
            else
                step_mode = obj.StepSizeMode;
            end
            
            if isscalar(nbhd_mode) && nbhd_mode == 1
                % 1 order neighborhood will not consider the current point
                nbhds = find1OrderNeighborhood(obj, options, ...
                    list_params, step_size, range, step_mode);
            elseif isscalar(nbhd_mode) && nbhd_mode == 2
                nbhds = find2OrderNeighborhood(obj, options, ...
                    list_params, step_size, range, step_mode);
            else % it's a vector specifying the size of nbhd
                nbhds = findBigNeighborhood(obj, options, ...
                    list_params, step_size, range, step_mode, ...
                    nbhd_mode);
            end
        end
        
        function nbhds = find1OrderNeighborhood(obj, options, ...
                list_params, step_size, range, modes)
            nbhds = [];
            options_start = options;
            for i = 1:length(list_params)
                options = options_start;
                name = list_params{i};
                forw_val = obj.findNeighboringValue(options_start.(name),...
                    step_size, i, 1, modes, range);
                back_val = obj.findNeighboringValue(options_start.(name),...
                    step_size, i, 0, modes, range);
                
                options.(name) = forw_val;
                nbhds = [nbhds, options];
                options.(name) = back_val;
                nbhds = [nbhds, options];
            end
        end
        
        function nbhds = find2OrderNeighborhood(obj, options, ...
                list_params, step_size, range, modes)
            % modes are the step mode (multiplication or addition)
            nbhd_size = 3 * ones(size(list_params));
            nbhds = findBigNeighborhood(obj, options, ...
                list_params, step_size, range, modes, nbhd_size);
        end
        
        function nbhds = findBigNeighborhood(obj, options, ...
                list_params, step_size, range, modes, nbhd_size)
            % num_steps takes the form like (3,3) for 2 parameters where
            % 3 values including the current value are considered as nbhd
            num_steps = floor(nbhd_size/2);
            options_start = options;
            params = cell(1,length(list_params));
            for i = 1:length(list_params)
                name = list_params{i};
                curr_val = options_start.(name);
                vals = curr_val;
                for j = 1:num_steps(i)
                    back_val = obj.findNeighboringValue(vals(1), ...
                        step_size, i, 0, modes, range);
                    vals = cat(1,back_val,vals);
                end
                for j = 1:num_steps(i)
                    forw_val = obj.findNeighboringValue(vals(end), ...
                        step_size, i, 1, modes, range);
                    vals = cat(1,vals,forw_val);
                end
                params{i} = vals;
            end
            % cartprod can't handle duplicate entries
%             C = cartprod(params{:}); 
            C = ParameterAnalysis.cartesian(params{:});

            nbhds = repmat(options_start, [size(C,1),1]);
            for i = 1:size(C,1)
                for j = 1:length(list_params)
                    nbhds(i).(list_params{j}) = C(i,j);
                end
            end
        end
                
        
        function new_val = findNeighboringValue(obj, curr_val, ...
                step_size, i, direction, modes, range)
            order_of_magnitude = @(x) floor(log10(x));
            
            if isnumeric(step_size)
                size = step_size(i);
                mode = modes(i);
                new_val = obj.calcNewValue(direction, curr_val, size, mode);
            elseif isequal(step_size,'medium')
                possible_vals = [1,2,5];
                size = 2;
                mode = 1;
                new_val = obj.calcNewValue(direction, curr_val, size, mode);
                order = order_of_magnitude(new_val);
                new_val = new_val/(10^order);
                [~,ind] = min(abs(new_val - possible_vals));
                new_val = possible_vals(ind);
                new_val = new_val * (10^order);
            end
            
            if direction == 1
                upper_bd = range(2,i);
                new_val = min(new_val, upper_bd);
            else                
                lower_bd = range(1,i);
                new_val = max(new_val, lower_bd);
            end
        end
        
        function new_val = calcNewValue(obj, direction, curr_val, size, mode)
            if mode == 1 % multiplication
                if direction == 1
                    new_val = curr_val * size; 
                else
                    new_val = curr_val / size; 
                end
            elseif mode == 2 % addition
                if direction == 1
                    new_val = curr_val + size; 
                else
                    new_val = curr_val - size; 
                end
            end
        end
    end
    methods (Static)
        function C = cartesian(varargin)
            args = varargin;
            n = nargin;
            
            [F{1:n}] = ndgrid(args{:});
            
            for i=n:-1:1
                G(:,i) = F{i}(:);
            end
            
            C = unique(G , 'rows');
        end

        function t = isOptionsEqual(options1, options2, list_params)
            t = 1;
            for i = 1:length(list_params)
                if options1.(list_params{i}) ~= options2.(list_params{i})
                    t = 0;
                    break;
                end
            end
        end
    end
end

