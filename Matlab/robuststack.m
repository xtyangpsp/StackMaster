function [dout, stat]= robuststack(din,par)
%Robust stacking scheme is based on Pavlis and Vernon,
%Computers & Geosciences, 2010.
%   din: seismic traces to stack. By default, each column contains one
%       trace. The default stacking dimension is 2 (columns). This can be
%       specified in "par", which is a struct containing all processing
%       parameters, as par.stackdim.
%   par: This is a struct containing all processing parameters.
%        par.win - stacking window specified as [min_index max_index] for
%           Robust stacking. The stacking weight for each trace will be
%           computed within the stacking window.
%        par.maxit - maximum number of iteration in iterative stacking
%        par.stackdim - stacking dimension in DIN, default is the second
%           dimension (columns).
%        par.reference - Reference trace uses to initiate the iterative
%           stacking. This is optional. If not specified, the program used
%           median as the starting stacking.
%        par.stacktype - For Robust stacking, this could be 'robust' or 
%           'corrcoef' when computing the weights.
%
% Output arguments:
%   dout: the stack of all traces. The length of 'dout' is the same as the
%       individual trace in 'din'
%   stat: this is a struct containing all statistical information for the
%       stacking operation.
%       stat.w - stacking weight for the final iteration.
%       stat.nit - number of iteration used to generate the final stack.
%
%Wrote by: Xiaotao Yang Jan 2020
%
par0=struct('maxit',100,'stackdim',2,'verbose',0,'stacktype','robust','win',[],'reference',[]);
if nargin < 1
    error('Not enough input arguments');
end

if nargout > 2
    error('Too many output arguments');
end

if nargin == 1;par=par0;end

if ~isfield(par,'maxit'); par.maxit = par0.maxit;end
if ~isfield(par,'stacktype'); par.stacktype = par0.stacktype;end
if ~isfield(par,'stackdim'); par.stackdim = par0.stackdim;end
if ~isfield(par,'win'); par.win = par0.win;end
if ~isfield(par,'verbose'); par.verbose = par0.verbose;end
if ~isfield(par,'reference'); par.reference = par0.reference;end

epsilon = 0.00001; %used to mark convergence when the stack beam changes less than this value.
%this number is used only in RobustSNR stack type. See Pavlis and Vernon,
%Computers & Geosciences, 2011 for details.
largenumber = 100000; %used as the initial value of model norm.

if par.stackdim ==1
    din=din';
end

[nsamp,ndata]=size(din);
if ndata<=1
    dout=din;
    itecount=1;
    w=1;
else
    if isempty(par.win) 
        if par.verbose;warning('par.win not specified, use the whole data length instead!');end
        par.win=[1,nsamp];
    end
    idmin=par.win(1);
    idmax=par.win(2);

    data=din(idmin:idmax,:);
    dout=zeros(nsamp,1);

    if isempty(par.reference)
        s = nanmedian(data,2); %start from median stack.
    else
        s = par.reference(idmin:idmax);
    end

    %get initial deltad: left hand side term in Pavlis and Vernon, equation (6)
    deltad = largenumber;
    itecount = 0;
    w = zeros(ndata,1);
    ampscale = zeros(ndata,1);
    % figure;hold on;

    if strcmp(par.stacktype,'robust')
        while deltad > epsilon && itecount < par.maxit
            for i = 1:ndata
                dtemp = data(:,i);
                if range(dtemp) == 0 || ~isempty(find(isnan(dtemp),1))
                    w(i) = 0; % for zero trace or NaN trace, use zero weight.
                else                
                    ampscale(i)=abs(dot(s,dtemp));
                    d = norm(dtemp);
                    r = norm(dtemp - s);%norm(dtemp - (ampscale(i)/abs(dot(dtemp,dtemp)))*s) + 1;
                    w(i) = ampscale(i)./(d.*r);
    %                 w(i) = abs(dot(s,dtemp))./(norm(dtemp).*norm(dtemp - s));
    %                 w(i) = w(i)^2;
                end
            end

            % normalize weight
            w=w./nansum(w);
            stemp = zeros(idmax - idmin +1,1);
            for i = 1:ndata
                dtemp = data(:,i);
                if abs(range(dtemp)) > 0 && isempty(find(isnan(dtemp),1))
                    stemp = stemp + w(i)*data(:,i);
                end
            end
    %         stemp=stemp/norm(stemp);
            slast = s;
            s = stemp;
    %         plot(s);
            deltad = norm(s - slast,1)/(norm(s)*nsamp);
            itecount = itecount + 1;
        end
    elseif strcmp(par.stacktype,'corrcoef')
        data=din(idmin:idmax,:);
        dout=zeros(nsamp,1);

        while deltad > epsilon && itecount < maxit
            for i = 1:ndata
                dtemp = data(:,i);
                if range(dtemp) == 0 || ~isempty(find(isnan(dtemp),1))
                    w(i) = 0; % for zero trace or NaN trace, use zero weight.
                else       
                    ce=corrcoef(s,dtemp);
                    w(i) = ce(1,2)*heaviside(ce(1,2));
                end
            end

            % normalize weight
            w=w./nansum(w);
            stemp = zeros(idmax - idmin +1,1);
            for i = 1:ndata
                dtemp = data(:,i);
                if abs(range(dtemp)) > 0 && isempty(find(isnan(dtemp),1))
                    stemp = stemp + w(i)*data(:,i);
                end
            end
            slast = s;
            s = stemp;
            deltad = norm(s - slast,1)/(norm(s)*nsamp);
            itecount = itecount + 1;
        end
    end

    for j = 1:ndata
        dtemp = din(:,j);
        if abs(range(dtemp)) > 0 && isempty(find(isnan(dtemp),1))
            dout = dout + w(j)*din(:,j);
        end
    end
end

if par.stackdim ==1
    dout=dout';
end
%assemble statistical parameters
stat.w = w;
stat.nit = itecount;
stat.par=par;
return;
end