function [dout, stat]=nrootstack(din,par)
%% N^th root stacking (default par.N=2).
%   din: seismic traces to stack. By default, each column contains one
%       trace. The default stacking dimension is 2 (column). This can be
%       specified in "par", which is a struct containing all processing
%       parameters, as par.stackdim.
%   par: This is a struct containing all processing parameters.
%       par.N - N^th root (default is 2).
%       par.stackdim - stacking dimension in DIN, default is the second
%           dimension (columns).
%
% Output arguments
%   dout: the stack of all traces. The length of 'dout' is the same as the
%       individual trace in 'din'
%   stat: this is a struct containing all statistical information for the
%       stacking operation.
% 
% Adapted by Xiaotao Yang @ Harvard University March 2020
% 
% REFERENCES:
% Millet, F., Bodin, T., & Rondenay, S. (2019). Multimode 3-D Kirchhoff 
%   migration of receiver functions at continental scale. Journal of Geophysical
%   Research: Solid Earth, 124, 8953?8980. https://doi.org/10.1029/ 2018JB017288
% Ruckemann (2012). Comparison of Stacking Methods Regarding Processing and Computing
%   of Geoscientific Depth Data
%
par0=struct('N',2,'stackdim',2,'verbose',0);

if nargin < 1
    error('Not enough input arguments');
end

if nargout > 2
    error('Too many output arguments');
end

if nargin == 1;par=par0;end

if ~isfield(par,'stackdim'); par.stackdim = par0.stackdim;end
if ~isfield(par,'verbose'); par.verbose = par0.verbose;end
if ~isfield(par,'N'); par.N = par0.N;end
messagebase='NROOTSTACK ';

if mod(log2(par.N),1)>0.0
    error([messagebase 'par.N should satisfy N=2^p, p = 0, 1, ... '])
end
if par.stackdim ==1
    din=din';
end
% Nan check
if any(any(isnan(din)))
    warning([messagebase 'Trace contain NaN value. Ignore the trace.']);
    nancol = any(isnan(din));
    din = din(:, ~nancol);
end

[nsamp,ndata]=size(din);
dout=zeros(nsamp,1);
for i=1:ndata
    dat=din(:,i);
    dout=dout+sign(dat).*abs(dat).^(1/par.N);
end
dout=dout/ndata;
dout=sign(dout).*abs(dout).^par.N;
stat.par=par;
return;
end