function [dout, stat]= seisstack(din,method,par)
%This is a wrapper function calling other stacking functions.
%Usage: SEISSTACK(din,method,par)
% Input arguments:
%   din: seismic traces to stack. By default, each column contains one
%       trace. The default stacking dimension is 2 (columns). This can be
%       specified in "par", which is a struct containing all processing
%       parameters, as par.stackdim.
%
%   method: string specifying the stacking method. Currently implemented
%   methods includes 'linear', 'selective','robust','phase','autovariance'.
%       'linear' - sum all traces and get the average. The traces might be
%           normalized before stacking. Specify this as par.normalize=1.
%       'selective' - stack of selected traces with given correlation
%           coefficient threshold, as par.ccmin.
%       'robust' - Robust stacking that starts from median of all
%           traces and compute the weight based on similrity of the trace
%           to the stack. This is a iterative method. There is an option to
%           specify a time window (as the indices of the starting and
%           ending points) as par.win=[startindex,endindex].
%       'pws' - phase-weighted stacking method.
%       'acf' - Adaptive covirance filter.
%
%   par: This is a struct containing all processing parameters, customized
%   for each stacking method. Here is a list of the members:
%        par.normalize - specifies the normalization option for LINEAR 
%           stacking. If 1, the program would normalize each trace beforeing averaging. 
%        par.win - stacking window specified as [min_index max_index] for
%           Robust stacking. The stacking weight for each trace will be
%           computed within the stacking window.
%        par.maxit - maximum number of iteration in iterative stacking
%        par.stackdim - stacking dimension in DIN, default is the second
%           dimension (columns).
%        par.reference - Reference trace used to initiate the iterative
%           stacking or used in computing the correlation coefficients for 
%           selective stacking. This is optional. If not specified, the
%           program uses median as the starting stacking for Robust stacking 
%           and the linear stack of all traces for selective stacking.
%        par.stacktype - this applies for Robust stacking, which could be
%           'robust' or 'corrcoef' when computing the weights.
%        par.ccmin - minimum correlation coefficient as the cutoff
%           threshold.
%
% Output arguments:
%   dout: the stack of all traces. The length of 'dout' is the same as the
%       individual trace in 'din'
%   stat: this is a struct containing all statistical information for the
%       stacking operation. Different method will save different
%       statistical data. All stacking methods save stat.method, stat.t (the
%       cpu time spent for the stacking operation), and stat.w (stacking
%       weight). stat.par will be saved for reference. stat.method and stat.t
%       are set by this wrapper function. Additional parameters are saved for
%       the following methods:
%       'selective' stacking saves stat.ccstack as the correlation
%           coefficient between the final stack and the reference trace, 
%           stat.ar (acceptance ratio), and stat.ccall (correlation coefficient
%           between each trace and the reference.
%       'robust' stacking saves stat.nit (number of iteration used to generate the
%           final stack) 
%       'pws' stacking saves TBA
%       'acf' stacking saves TBA
%       
%
%Wrote by: Xiaotao Yang Jan 2020
%

par0=struct('stackdim',2); %set default parameters.

if nargin < 2
    error('Not enough input arguments');
end

if nargout > 3
    error('Too many output arguments');
end

if nargin == 2;par=par0;end

tic;
if strcmp(method,'linear')
    [dout, stat]=linstack(din, par);
elseif strcmp(method,'robust')
    [dout, stat]=robuststack(din, par);
elseif strcmp(method,'pws')
    [dout, stat]=pwstack(din, par);
elseif strcmp(method,'tf-pws')
    par.stacktype='tf-pws';
    [dout, stat]=pwstack(din, par);
elseif strcmp(method,'acf')
    [dout, stat]=acfstack(din, par);
elseif strcmp(method,'selective')
    [dout, stat]=selectivestack(din, par);
elseif strcmp(method,'nroot')
    [dout, stat]=nrootstack(din, par);
else
    error(['Stacking method NOT recoganized: ',method]);
end
stat.t = toc;
stat.method = method;

return;
end