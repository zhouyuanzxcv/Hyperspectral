classdef CachedCalculationResults
    %CACHEDCALCULATIONRESULTS Summary of this class goes here
    %   Detailed explanation goes here
    % Usage:
    %     cached_params = [];
    %     cached_params.Y = Y;
    %     cached_params.k = k;
    %     cached_file = '../../../../CachedResults/cached_nmf_adjacency_graph.mat';
    %     ccr = CachedCalculationResults(cached_file);
    %     A = ccr.GetMatchedResult(cached_params);
    %     if isempty(A)
    %         [A,~,~] = create_adjacency_graph(Y,'nn',k,1);
    %         if N*B > 1e6
    %             ccr.AddResult(cached_params,A);
    %         end
    %     end

    
    properties
        FileName
        MaxSignatureSize
        MaxEntries
    end
    
    methods
        function obj = CachedCalculationResults(file_name,max_param_size)
            if nargin < 2
                max_param_size = 1e5;
            end
            obj.FileName = file_name;
            obj.MaxSignatureSize = max_param_size;
            obj.MaxEntries = 30;
        end
        
        function AddResult(obj, params, result)
            [list_params,list_results] = LoadCachedResults(obj);
            if length(list_params) >= obj.MaxEntries
                list_params(1) = [];
                list_results(1) = [];
            end
            list_params = [list_params,{obj.GetSignature(params)}];
            list_results = [list_results,{result}];
            save(obj.FileName,'list_params','list_results');
        end
        
        function matched = GetMatchedResult(obj,params)
            [list_params,list_results] = LoadCachedResults(obj);
            signature = obj.GetSignature(params);

            matched = [];
            for i = 1:length(list_params)
                if isequal(list_params{i},signature)
                    matched = list_results{i};
                    break;
                end
            end
            
        end
        
        function [list_params,list_results] = LoadCachedResults(obj)
            if ~exist(obj.FileName,'file')
                list_params = {};
                list_results = {};
                save(obj.FileName,'list_params','list_results');
            end
            S = load(obj.FileName);
            list_params = S.list_params;
            list_results = S.list_results;
        end
        
        function new_params = GetSignature(obj,params)
            new_params = params;
            names = fieldnames(new_params);
            for i = 1:length(names)
                X = new_params.(names{i});
                x = X(:);
                if length(x) > obj.MaxSignatureSize
                    step = ceil(length(x) / obj.MaxSignatureSize);
                    x = x(1:step:end);
                end
                new_params.(names{i}) = x;
            end
        end
        
    end
    
end

