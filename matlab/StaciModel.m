classdef StaciModel < handle
    %STACIMODEL In-memory MATLAB interface to a STACI hydraulic network.
    % All public numerical properties use SI units.

    properties (SetAccess = private)
        Info struct
        NetworkFile string
    end

    properties (Access = private)
        Handle uint64 = uint64(0)
    end

    methods
        function object = StaciModel(networkFile)
            arguments
                networkFile {mustBeTextScalar}
            end
            if exist('staci_mex', 'file') ~= 3
                build_staci_mex();
            end
            object.NetworkFile = string(networkFile);
            [object.Handle, object.Info] = staci_mex( ...
                'load', char(object.NetworkFile));
        end

        function delete(object)
            object.release();
        end

        function release(object)
            if object.Handle ~= 0
                try
                    staci_mex('release', object.Handle);
                catch exception
                    warning('STACI:MatlabRelease', ...
                        'Could not release STACI model: %s', exception.message);
                end
                object.Handle = uint64(0);
            end
        end

        function ids = nodeIds(object)
            object.requireOpen();
            ids = string(object.Info.node_ids(:));
        end

        function ids = linkIds(object, type)
            object.requireOpen();
            ids = string(object.Info.link_ids(:));
            if nargin >= 2 && strlength(string(type)) > 0
                requested = object.normalizeLinkType(string(type));
                types = string(object.Info.link_types(:));
                ids = ids(strcmpi(types, requested));
            end
        end

        function types = linkTypes(object)
            object.requireOpen();
            types = string(object.Info.link_types(:));
        end

        function values = getNodeProperty(object, ids, property)
            object.requireOpen();
            values = staci_mex('get_node_property', object.Handle, ...
                object.asCellText(ids), char(string(property)));
        end

        function setNodeProperty(object, ids, property, values)
            object.requireOpen();
            staci_mex('set_node_property', object.Handle, ...
                object.asCellText(ids), char(string(property)), double(values));
        end

        function values = getLinkProperty(object, ids, property)
            object.requireOpen();
            values = staci_mex('get_link_property', object.Handle, ...
                object.asCellText(ids), char(string(property)));
        end

        function setLinkProperty(object, ids, property, values)
            object.requireOpen();
            staci_mex('set_link_property', object.Handle, ...
                object.asCellText(ids), char(string(property)), double(values));
        end

        function status = solveHydraulics(object)
            object.requireOpen();
            status = staci_mex('solve_hydraulics', object.Handle);
        end

        function resetHydraulicState(object)
            object.requireOpen();
            staci_mex('reset', object.Handle);
        end

        function results = nodeResults(object)
            object.requireOpen();
            results = staci_mex('node_results', object.Handle);
        end

        function results = linkResults(object)
            object.requireOpen();
            results = staci_mex('link_results', object.Handle);
        end

        function results = solveSteadyQuality(object, mode)
            arguments
                object
                mode {mustBeTextScalar} = "both"
            end
            object.requireOpen();
            results = staci_mex('solve_steady_quality', object.Handle, ...
                lower(char(string(mode))));
        end

        function results = solveSteadyWaterAge(object)
            results = object.solveSteadyQuality("age");
        end

        function results = hydraulicSensitivity(object, elementId, property)
            object.requireOpen();
            results = staci_mex('hydraulic_sensitivity', object.Handle, ...
                char(string(elementId)), char(string(property)));
        end

        function tableValue = nodeTable(object)
            values = object.nodeResults();
            tableValue = table(string(values.id(:)), values.elevation_m(:), ...
                values.pressure_head_m(:), values.total_head_m(:), ...
                values.demand_m3s(:), ...
                'VariableNames', {'id', 'elevation_m', 'pressure_head_m', ...
                                  'total_head_m', 'demand_m3s'});
        end

        function tableValue = linkTable(object)
            values = object.linkResults();
            tableValue = table(string(values.id(:)), string(values.type(:)), ...
                string(values.from_node(:)), string(values.to_node(:)), ...
                values.flow_rate_m3s(:), values.velocity_ms(:), ...
                values.headloss_m(:), logical(values.status(:)), ...
                'VariableNames', {'id', 'type', 'from_node', 'to_node', ...
                                  'flow_rate_m3s', 'velocity_ms', ...
                                  'headloss_m', 'status'});
        end
    end

    methods (Access = private)
        function requireOpen(object)
            if object.Handle == 0
                error('STACI:InvalidHandle', ...
                    'This StaciModel object has already been released.');
            end
        end
    end

    methods (Static, Access = private)
        function values = asCellText(ids)
            values = cellstr(string(ids(:)));
            if isempty(values)
                error('STACI:EmptySelection', ...
                    'At least one node or link id is required.');
            end
        end

        function type = normalizeLinkType(type)
            switch lower(type)
                case {"pipe", "cso"}
                    type = "Cso";
                case {"pump", "szivattyu"}
                    type = "Szivattyu";
                case {"powerpump", "epanetpowerpump"}
                    type = "EpanetPowerPump";
                case {"valve", "jelleggorbesfojtas"}
                    type = "JelleggorbesFojtas";
            end
        end
    end
end
