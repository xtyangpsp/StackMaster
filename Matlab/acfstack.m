function [dout, stat]= acfstack(din, par)
%% Stacking after applying adaptive covariance filter (ACF) following Nakata at el. (2015).
%   din: seismic traces to stack. By default, each column contains one
%       trace. The default stacking dimension is 2 (column). This can be
%       specified in "par", which is a struct containing all processing
%       parameters, as par.stackdim.
%   par: This is a struct containing all processing parameters. It has the
%       following fields:
%       verbose - how verbose the program is, 1 or 0;
%       stackdim - stacking dimension in DIN, default is the second
%           dimension (columns).
%       window - Window length to apply adaptive filter in number of samples.
%       harshness - value to adjust the filter harshness;
%       overlap - Percentage overlap between windows (0 to 1).
%
% Output arguments
%   dout: the stack of all traces. The length of 'dout' is the same as the
%       individual trace in 'din'
%   stat: this is a struct containing all statistical information for the
%       stacking operation.
% 
% Adapted from Tim Clements' implementation in Julia 
% (https://github.com/tclements/SeisNoise.jl/blob/master/src/stacking.jl).
% The performance is optimized for MATLAB through reducing for loops.
%
% by Xiaotao Yang @ Harvard University March 2020
%
% REFERENCE:
%Nakata, N., Chang, J. P., Lawrence, J. F., & Boué, P. (2015). Body wave 
% extraction and tomography at Long Beach, California, with ambient?noise 
% interferometry. Journal of Geophysical Research: Solid Earth, 120(2), 1159-1173.

%

%default parameters
par0=struct('stackdim',2,'verbose',0,'harshness',2.0,'overlap',0.9);

if nargin < 1
    error('Not enough input arguments');
end

if nargout > 2
    error('Too many output arguments');
end

if nargin == 1;par=par0;end
if ~isfield(par,'stackdim'); par.stackdim = par0.stackdim;end
if ~isfield(par,'verbose'); par.verbose = par0.verbose;end
if ~isfield(par,'harshness'); par.harshness = par0.harshness;end
if ~isfield(par,'overlap'); par.overlap = par0.overlap;end

if par.stackdim ==1
    din=din';
end

messagebase=['ACFSTACK '];
if ~isfield(par,'window')
    error([messagebase '[window] is not set in par!']);
end
% par.window = 2^(log2(par.window) - mod(log2(par.window),1) + 1);
% Nan check
if any(any(isnan(din)))
    warning([messagebase 'Trace contain NaN value. Ignore the trace.']);
    nancol = any(isnan(din));
    din = din(:, ~nancol);
end

[nsamp, ndata] = size(din);
if ndata <=1
    dout=din;
else
    overlap_factor=1 - par.overlap;
    window_step = round(par.window * (1 - par.overlap));
    minind = 1:window_step:nsamp - par.window;

    % allocate out array
    dout_temp = zeros(nsamp, ndata);
    W=repmat(tukeywin(par.window,overlap_factor/4),1,ndata);
    % loop through each window
    for ii=1:length(minind)
        dout_temp(minind(ii):minind(ii)+par.window-1,:) =  dout_temp(minind(ii):minind(ii)+par.window-1,:) + ...
            ACF_kernel(din(minind(ii):minind(ii)+par.window-1,:).*W,par.harshness);
    end
    % windows at right edge
    reverseind = nsamp:-window_step:nsamp-par.window;
    for ii=1:length(reverseind)
        dout_temp(reverseind(ii)-par.window+1:reverseind(ii),:) = dout_temp(reverseind(ii)-par.window+1:reverseind(ii),:) + ...
            ACF_kernel(din(reverseind(ii)-par.window+1:reverseind(ii),:).*W,par.harshness);
    end
    % get the linear stack of the data after applying ACF.
    dout=linstack(dout_temp);
end
if par.stackdim ==1
    dout=dout';
end
stat.par=par;
return;

end

function Aout=ACF_kernel(A, g)
    if nargin==1
        g=2.0;
    end
    [Nrows, Ncols] = size(A);
    % fft the 2D array
    spec = fft(A,Nrows,1);

    % create auto-covariance function

    S2=sum(spec.* conj(spec),2);
    S1=sum(spec.*sum(conj(spec),2),2);
%     S1=sum(spec,2).*sum(conj(spec),2) - S2;
%     Nspec = size(spec,1);
%     S2 = complex(zeros(Nspec,1));
%     for ii=1:Ncols
%        S2 = S2 + spec(:,ii).*conj(spec(:,ii)); 
%     end
%     S1 = complex(zeros(Nspec,1));
% 
%     for ii = 1:Ncols
%         S1 = S1 + spec(:,ii).*sum(conj(spec),2);
% %         for jj = 1:Ncols
% %             S1 = S1 + spec(:,ii).*conj(spec(:,jj));
% %         end
%     end
    
    % make ifft
    Aout = ifft(spec.*(abs((S1 - S2) ./ (S2 * (Ncols -1))).^g),Nrows,1);
    
    return;
end