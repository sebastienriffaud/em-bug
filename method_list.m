function [methods,labels,names,markers] = method_list(problem)
%METHOD_LIST List of available low-rank methods and plotting metadata.
%   [methods,labels,names,markers] = METHOD_LIST(problem) returns function
%   handles and labels for the two BUG integrators.
%
%   Input
%   -----
%   problem : "ou" or "multiplicative".
%
%   Outputs
%   -------
%   methods : cell array of integrator handles.
%   labels  : string array of plot legend labels.
%   names   : string array of short method identifiers.
%   markers : string array of plotting markers.

    switch problem
        case {"ou","multiplicative"}
            methods = {@bug_standard,@bug_closed};
            labels = ["Standard BUG","Closed EM--BUG"];
            names = ["BUG_standard","BUG_closed"];
            markers = ["o","x"];

        otherwise
            error('method_list:UnknownTestcase','Unknown testcase "%s".',problem)
    end

end
