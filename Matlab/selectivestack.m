function [dout, stat] = selectivestack(din, par)
%% Stacking with given threshold of correlation coeffcient between individual trace and the reference
% This function was modified from Kurama Okubo's function with the same name.
% The original function was downloaded here: https://github.com/kura-okubo/ccstack.m/blob/master/ccstack.m.
% Process flow:
%   1. compute Pearson's correlation coefficient between each trace and
%   ref.
%   2. thresholding out them and stack
%
%   din: seismic traces to stack. By default, each column contains one
%       trace. The default stacking dimension is 2 (column). This can be
%       specified in "par", which is a struct containing all processing
%       parameters, as par.stackdim.
%   par: This is a struct containing all processing parameters.
%       par.ccmin - minimum correlation coefficient as the cutoff
%           threshold.
%       par.stackdim - stacking dimension in DIN, default is the second
%           dimension (columns).
%       par.reference - Reference trace used to in computing the correlation 
%           coefficients for selective stacking. This is optional. If not 
%           specified, the program uses the linear stack of all traces.
%       par.maxit - maximum number of iteration in iterative stacking
%
% Output arguments:
%   dout: the stack of all traces. The length of 'dout' is the same as the
%       individual trace in 'din'
%   stat: this is a struct containing all statistical information for the
%       stacking operation.
%       stat.w - stacking weight.
%       stat.ar - acceptance ratio.
%       stat.ccall - correlation coefficient between each trace and the reference.
% 
% Adapted by Xiaotao Yang @ Harvard University March 2020
%
par0=struct('ccmin',0,'stackdim',2,'verbose',0,'reference',[],'maxit',100);
epsilon = 0.00001; %used to mark convergence when the stack beam changes less than this value.
largenumber = 100000; %used as the initial value of model norm.
messagebase='SELECTIVESTACK ';
if nargin < 1
    error('Not enough input arguments');
end

if nargout > 2
    error('Too many output arguments');
end

if nargin == 1;par=par0;end

if ~isfield(par,'maxit'); par.maxit = par0.maxit;end
if ~isfield(par,'stackdim'); par.stackdim = par0.stackdim;end
if ~isfield(par,'verbose'); par.verbose = par0.verbose;end
if ~isfield(par,'reference')
    if par.verbose
        warning([messagebase 'reference trace not set in par. use thelinear stack as the reference. Specify as par.reference.']);
    end
    par.reference = par0.reference;
end
if ~isfield(par,'ccmin')
    if par.verbose
        warning([messagebase 'ccmin not set in par. use 0.0 as default. Specify as par.ccmin for other values.']);
    end
    par.ccmin = par0.ccmin;
end
if par.stackdim ==1
    din=din';
end
% Nan check
if any(any(isnan(din)))
    warning([messagebase 'selectivestack(): Trace contain NaN value. Ignore the trace.']);
    nancol = any(isnan(din));
    din = din(:, ~nancol);
end

[nsamp,ndata]=size(din);
if isempty(par.reference)
    ref=nanmean(din,2);
%     ref=nanmedian(din,2);
else
    ref=par.reference;
end

%update iteratively
%get initial deltad: left hand side term in Pavlis and Vernon, equation (6)
deltad = largenumber;
itecount = 0;
slast = ref;
while deltad > epsilon && itecount < par.maxit
    w = zeros(ndata,1);

    cclist = nan(ndata, 1);

    for i = 1:ndata
        temp = corrcoef(din(:, i), slast);
        cclist(i) = temp(1,2);
    end
    ind = find(cclist >= par.ccmin);
    dout = mean(din(:, ind), 2);
    w(ind)=1;
    
    deltad = norm(slast - dout,1)/(norm(slast)*nsamp);
    slast = dout;
    itecount=itecount+1;
end
if par.stackdim ==1
    dout=dout';
end

stat.ar = 100*length(ind)/ndata; 
stat.ccall=cclist;
if par.verbose
    fprintf("Acceptance ratio: %4.2f [%%]\n", stat.ar);
end

stat.par=par;
stat.w=w;
stat.nit=itecount;
return;
end