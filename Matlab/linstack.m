function [dout,stat]=linstack(din,par)
%   din: seismic traces to stack. By default, each column contains one
%       trace. The default stacking dimension is 2 (columns). This can be
%       specified in "par", which is a struct containing all processing
%       parameters, as par.stackdim.
%   par: This is a struct containing all processing parameters.
%       par.normalize - specifies the normalization option for LINEAR 
%           stacking. If 1, the program would normalize each trace beforeing averaging.
%       par.stackdim - stacking dimension in DIN, default is the second
%           dimension (columns).
%
% Output arguments:
%   dout: the stack of all traces. The length of 'dout' is the same as the
%       individual trace in 'din'
%   stat: this is a struct containing all statistical information for the
%       stacking operation.
%       stat.w - stacking weight, which is 1/max(abs(din(:,j))).
% 
% Xiaotao Yang @ Harvard University March 2020
par0=struct('stackdim',2,'verbose',0,'normalize',0);
if nargin < 1
    error('Not enough input arguments');
end

if nargout > 2
    error('Too many output arguments');
end

if nargin == 1;par=par0;end

if ~isfield(par,'normalize'); par.normalize = par0.normalize;end
if ~isfield(par,'stackdim'); par.stackdim = par0.stackdim;end
if ~isfield(par,'verbose'); par.verbose = par0.verbose;end

if par.stackdim ==1
    din=din';
end

ndata=size(din,2);
if par.normalize    
    dintemp=nan(size(din));
    for i=1:ndata
        stat.w(i) = 1./nanmax(abs(din(:,i)));
        dintemp(:,i)=din(:,i)*stat.w(i);
    end
    dout = nanmean(dintemp,par.stackdim);
else
    dout = nanmean(din,par.stackdim);
    stat.w = ones(ndata,1);
end
if par.stackdim ==1
    dout=dout';
end
stat.par=par;
return;
end