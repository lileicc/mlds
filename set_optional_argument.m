% This source code is (c) Copyright by Lei Li, Mark Rogers.
% All rights preserved.
%
% Permission is granted to use it for non-profit purposes,
% including research and teaching. For-profit use requires
% the express consent of the author (leili@cs.berkeley.edu).
%
% Details in the following paper:
%   Mark Rogers, Lei Li and Stuart J. Russell (2013),
%     "Multilinear Dynamical Systems for Tensor Time Series",    
%     In Advances in Neural Information Processing Systems 26. 
%
function argval = set_optional_argument(argname, default_argval, optional_args)
%
% set optional argument with name argname.  Default value is default_argval.
%
% @author: Mark Rogers
% @last modified date: 2013/12/13
%
a = find(strcmp(argname, optional_args), 1);
if isempty(a)
  argval = default_argval;
else
  argval = optional_args{a+1};
end
