% WorkerObjWrapper - manage persistent state on MATLABPOOL workers
%
% The WorkerObjWrapper is designed for situations where a piece of data is
% needed multiple times inside the body of a PARFOR loop or an SPMD block, and
% this piece of data is both expensive to create, and does not need to be
% re-created multiple times. Examples might include: database connection
% handles, large arrays, and so on.

% Copyright 2011-2013 The MathWorks, Inc.

classdef WorkerObjWrapper < handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public interface:
    properties ( Transient, Dependent, SetAccess = private )
        Value % Underlying value, accessible only on the workers
    end
    methods
        function obj = WorkerObjWrapper( ctor, args, dtor )
        % WORKEROBJWRAPPER - client-side constructor for worker object wrapper
        %
        %    W = WORKEROBJWRAPPER( X ) evaluated on the client creates a
        %    WorkerObjWrapper using the specified value X in the client's
        %    workspace. When W.Value is subsequently used within SPMD or PARFOR,
        %    it will have the value X.
        %
        %    W = WORKEROBJWRAPPER( C ) evaluated on the client creates a
        %    WorkerObjWrapper using the specified Composite C.
        %
        %    W = WORKEROBJWRAPPER( FCN, ARGSCELL ) creates a WorkerObjWrapper by
        %    invoking FCN inside an SPMD block, supplying ARGSCELL as
        %    arguments. I.e., W.Value has the value obtained from running
        %    FCN(ARGSCELL{:}) on the workers.
        %
        %    W = WORKEROBJWRAPPER( FCN, ARGSCELL, CLEANUPFCN ) creates the
        %    WorkerObjWrapper as in the previous case. When the variable W is
        %    cleared or goes out of scope on the client, the cleanup function
        %    CLEANUPFCN is invoked on the value. I.e., CLEANUPFCN(W.Value) is
        %    evaluated on the workers.
        %
        %    Example 1, simply copy a value to the workers once and use it
        %    multiple times:
        %
        %      w = WorkerObjWrapper( magic(5) ); % copies "magic(5)" to each worker
        %      spmd, size(w.Value), end          % returns [5 5] on each worker
        %
        %    Example 2, build the value on the workers:
        %
        %      w = WorkerObjWrapper( @rand, {5} ); % invokes "rand(5)" on each worker
        %      spmd, size(w.Value), end            % returns [5 5] on each worker
        %      spmd, isreplicated(w.Value), end    % returns false
        %
        %    Example 3, build the value on the workers based on labindex:
        %
        %      % invokes "rand(labindex)" on each worker:
        %      w = WorkerObjWrapper( @() rand(labindex), {} );
        %      spmd, size(w.Value,1) == labindex, end % true on each worker
        %
        %    Example 4, build the value on the workers, and use a cleanup function:
        %
        %      % build a function handle to open a numbered text file:
        %      fcn = @() fopen( sprintf( 'worker_%d.txt', labindex ), 'wt' );
        %
        %      % opens the file handle on each worker, specifying that fclose
        %      % will be used later to "clean up" the file handle created.
        %      w = WorkerObjWrapper( fcn, {}, @fclose );
        %
        %      % Run a parfor loop, logging to disk which worker operated on which
        %      % loop iterates
        %      parfor ii=1:10
        %         fprintf( w.Value, '%d\n', ii );
        %      end
        %
        %      clear w; % causes "fclose(w.Value)" to be invoked on the workers
        %      type worker_1.txt % see which iterates worker 1 got
        %

            if ~isempty( getCurrentTask() )
                error( 'The WorkerObjWrapper should be constructed on the client' );
            end
            if nargin < 2
                if isa( ctor, 'Composite' )
                    if ~ all( exist( ctor ) ) %#ok<EXIST> Composite.exist
                        error( 'The supplied Composite must have a value on every worker.' );
                    end
                end
                % "ctor" is actually the value, just map things through so that
                % workerInit is none the wiser. Do not wrap args in a cell so
                % that Composites can be correctly transmitted.
                args = ctor;
                ctor = @iReturnInput;
            end
            if nargin < 3
                dtor = [];
            end
            tmpId = WorkerObjWrapper.getNextID();
            WorkerObjWrapper.workerInit( tmpId, ctor, args, dtor );
            obj.ID = tmpId;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % "Public" implementation details - these methods must have public
    % visibility, but they add nothing useful to the API.
    methods
        function v = get.Value( obj )
        % V = GET.VALUE(OBJ) - retrieve the value, only on the workers
            if ~isempty( getCurrentTask() )
                assert( WorkerObjWrapper.Map.isKey( obj.ID ) );
                valdtor = WorkerObjWrapper.Map( obj.ID );
                v = valdtor{1};
            else
                v = 'Only available on the workers';
            end
        end
        function delete( obj )
        % DELETE(OBJ) - when fired on the client, cleans up the value on the workers
            if isempty( getCurrentTask() ) && ~isempty( obj.ID )
                WorkerObjWrapper.workerDelete( obj.ID );
            else
                % on the workers, do nothing - we hold on to the value until the client releases
                % its handle.
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Private implementation
    properties ( Access = private )
        ID = [] % Unique key used to access the data, uint32
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constant / static data
    properties ( GetAccess = private, Constant )
        Map = containers.Map( 'KeyType', 'uint32', 'ValueType', 'any' );
    end

    methods ( Access = private, Static )
        function id = getNextID()
        % GETNEXTID - return the next uint32 ID key
            persistent NEXT_ID;
            if isempty( NEXT_ID )
                NEXT_ID = uint32(1);
            end
            id = NEXT_ID;
            if NEXT_ID == intmax( class( NEXT_ID ) )
                NEXT_ID = 1;
            else
                NEXT_ID = NEXT_ID + 1;
            end
        end
        function workerInit( id, ctor, args, dtor )
        % WORKERINIT - called to set the value in the map
            spmd, ?WorkerObjWrapper; end % workaround first-time access problems
            spmd
                if iscell( args )
                    obj = ctor( args{:} );
                else
                    obj = ctor( args );
                end
                m = WorkerObjWrapper.Map;
                m(id) = { obj, dtor }; %#ok<NASGU> handle object
            end
        end
        function workerDelete( id )
        % WORKERDELETE - called when the client object goes out of scope
        % Cleans up data on the workers
            spmd
                assert( WorkerObjWrapper.Map.isKey( id ) );
                valdtor = WorkerObjWrapper.Map( id );
                WorkerObjWrapper.Map.remove( id );
                [val, dtor] = deal( valdtor{:} );
                if ~isempty( dtor )
                    dtor( val );
                end
                % Clear these variables so Composites are not created
                val = []; dtor = []; valdtor = []; %#ok<NASGU>
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iReturnInput - no-op function
function x = iReturnInput( x )
end
