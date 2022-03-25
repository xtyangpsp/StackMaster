function [dout, stat] = pwstack(din, par)
%% Phase-weighted stacking, with the option of time-frequency PWS.
%   din: seismic traces to stack. By default, each column contains one
%       trace. The default stacking dimension is 2 (column). This can be
%       specified in "par", which is a struct containing all processing
%       parameters, as par.stackdim.
%   par: This is a struct containing all processing parameters.
%       par.stacktype - This could be 'pws' [default] or 'tf-pws' (for 
%           time-frequency PWS, see Schimmel et al, GJI, 2007).
%       par.stackdim - stacking dimension in DIN, default is the second
%           dimension (columns).
%       par.pow - This is the power index in Equation (4) of Schimmel and 
%           Paulssen (1997) defining the transition between coherent and 
%           less coherent signal summation. The default is 2 used by examples
%           in the paper.
%
% Output arguments
%   dout: the stack of all traces. The length of 'dout' is the same as the
%       individual trace in 'din'
%   stat: this is a struct containing all statistical information for the
%       stacking operation.
%       stat.w - time series showing the phase weight.
% 
% Adapted by Xiaotao Yang @ Harvard University March 2020
% 
% REFERENCES:
% Original PWS:
%       Schimmel and Paulssen, 1997, Noise reductiion and detection of weak,
%           coherent signals through phase-weighted stacks, GJI
% tf-PWS: 
%       Schimmel and Gallart, 2007, Frequency-dependent phase coherence for
%           noise suppression in seismic array data, JGR
%       Schimmel, Stutzmann, and Gallart, 2011, Using instantaneous phase 
%           coherence for signal extraction from ambient noise data at a 
%           local to a global scale, GJI
%       This implementation mainly followed Baig et al. (2009) using DOST 
%           for S-transform. 
%
par0=struct('pow',2,'stackdim',2,'verbose',0,'stacktype','pws','waterlevel',0.01);

if nargin < 1
    error('Not enough input arguments');
end

if nargout > 2
    error('Too many output arguments');
end

if nargin == 1;par=par0;end

if ~isfield(par,'stackdim'); par.stackdim = par0.stackdim;end
if ~isfield(par,'verbose'); par.verbose = par0.verbose;end
if ~isfield(par,'stacktype'); par.stacktype = par0.stacktype;end
if ~isfield(par,'waterlevel'); par.waterlevel = par0.waterlevel;end
messagebase=['PWSTACK (' par.stacktype ') '];
if ~isfield(par,'pow')
    if par.verbose
        warning([messagebase 'pow not set in par. use ' num2str(par0.pow) ' as default. Specify as par.pow for other values.']);
    end
    par.pow = par0.pow;
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

stat.pw = zeros(nsamp,1);
dout = zeros(nsamp,1);

if strcmp(par.stacktype,'pws')
    %get the phase weight 
    pw=sum(exp(1i.*angle(din + 1i.*hilbert(din))),2);
    pw=abs(pw/ndata).^par.pow;
    dout = nanmean(din,2).*pw;
elseif strcmp(par.stacktype,'tf-pws')
    if mod(log2(nsamp),1)>0.0
        ns2=2^(log2(nsamp)-rem(log2(nsamp),1)+1);
        trimsamp=ns2 - nsamp;  
        trimstart=floor(trimsamp/2);
        trimend=trimsamp-trimstart;

        %append with zeros to form 2^N samples.
        if trimstart >0 && trimend >0
            din2=[zeros(trimstart,ndata);din;zeros(trimend,ndata)];
        elseif trimend >0
            din2=[din;zeros(trimend,ndata)];
        end
    else
        trimstart=0;
        din2 = din;
    end
    Sall=ST.dost(din2);%.*din2;
    pw=sum(Sall./(abs(Sall) + par.waterlevel*mean(abs(Sall),2)*ones(1,ndata)),2);
    pw=abs(pw/ndata).^par.pow;
    dout2=nanmean(din2,2).*pw;%
    dout=dout2(trimstart+1:trimstart+nsamp);
    pw=pw(trimstart+1:trimstart+nsamp);
end

if par.stackdim ==1
    dout=dout';
end
stat.w=pw;
stat.par=par;

return;
end