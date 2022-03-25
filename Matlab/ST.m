%
% Code by: U. Battisti and L. Riba
% 7th November 2015
%
%         ST class which implements O(N log N) algorithms to compute 
%             DOST, DCST and FFT both in 1d and 2d.      
% 
% A tutorial, called ST_tutorial is distributed together with this file
% A graphical explanation of the code can be found in the included ST_algorithms.pdf 
%
%                                 HOW TO CITE
% If this code helps you in your research, please cite
% [1]   U. Battisti, L. Riba, "Window-dependent bases for efficient
%       representations of the Stockwell transform",
%       Applied and Computational Harmonic Analysis, 23 February 2015,
%       http://dx.doi.org/10.1016/j.acha.2015.02.002.
%
%                                 CONTACTS
% If you find a bug or you need more information regarding the code feel free
% to contact Luigi Riba at: ribaluigi (at) gmail (dot) com
%
%
%              THE S-TRANSFORM, A BRIEF INTRODUCTION
% The Stockwell transform (S-transform), introduced by R. G. Stockwell in 1996 [6],
% can be seen as an hybrid between the wavelet and Gabor transform
% (short time Fourier transform). The link between these transforms
% becomes evident in the multidimensional case as shown in [2] and [3].
% In practice, we can think of the S-transform as a short time Fourier transform
% in which the analyzing window gives better time localization for high frequencies
% and better time localization for low frequencies.
%
% The S-transform is computationally expensive. In [5], R. G. Stockwell himself 
% showed how to decompose a signal on a basis, called "DOST basis", which is linked 
% to the S-transform.
% 
% Wang and Orchard found a FFT-fast algorithm to compute the DOST coefficients, 
% see [7] and [8]. This made feasible the use of the DOST for 2 dimensional signals.
%
% The mathematical description of the fast DOST algorithmhs is based on the Plancharel 
% formula (see [1], Proposition 13 and Remark 4).
%
% In this code we also implemented the Discrete Cosine Stockwell Transform (DCST),
% which follows the DOST construction but uses the DCT instead of the FFT.
%
%                                REFERENCES
% [1]   U. Battisti, L. Riba, "Window-dependent bases for efficient
%       representations of the Stockwell transform",
%       Applied and Computational Harmonic Analysis, 23 February 2015,
%       http://dx.doi.org/10.1016/j.acha.2015.02.002.
% [2]   L. Riba, M. Wong, Continuous inversion formulas for multi-dimensional
%       modified Stockwell transforms, Integral Transforms Spec. Funct. 2015
% [3]   L. Riba, Multi-dimensional Stockwell transforms and applications,
%       PhD thesis, Università degli Studi di Torino, Italy, 2014.
% [4]   R.G. Stockwell, "Why use the S-Transform", Pseudo-differential
%       operators partial differential equations and time-frequency
%       analysis, vol. 52 Fields Inst. Commun., pages 279--309,
%       Amer. Math. Soc., Providence, RI 2007;
% [5]   R.G. Stockwell, "A basis for efficient representation of the
%       S-transform", Digital Signal Processing, 17: 371--393, 2007;
% [6]   R.G. Stockwell, L. Mansinha, R.P. Lowe, Localization of the complex spectrum:
%       the S transform, IEEE Trans. Signal Process. 44 (1996) 998–1001.
% [7]   Y. Wang and J. Orchard, "Fast-discrete orthonormal
%       Stockwell transform", SISC: 31:4000--4012, 2009;
% [8]   Y. Wang, "Efficient Stockwell transform with applications to
%       image processing", PhD thesis, University of Waterloo,
%       Ontario Canada, 2011;
%
%                             FUNCTION DESCRIPTIONS
%
%   ST.dost()
%     computes the DOST coefficients of a given signal
%     input: vector (or matrix) of size 2^k;
%     output: the DOST coefficients of the vector (of the columns of the matrix)
%
%   ST.rearrangeDost()
%     rearranges the DOST coefficients to make them more readable
%     input: a vector of DOST coefficients (length n=2^k)
%     output: a n x n matrix which describes the time frequency plane using the
%     DOST coefficients. See ST_algorithms.pdf for a graphical explanation.
%
%   ST.idost()
%     computes the inverse DOST
%     input: vector (or matrix) of size 2^k
%     output: the iDOST coefficients of the vector (of the columns of the matrix)
%
%   ST.dost2()
%     computes the 2 dimensional dost of a given matrix (ex: image) taking 
%	  the dost on the columns and then on the rows
%     input: matrix of size 2^k
%     output: dost2 coefficients of the matrix
%
%   ST.idost2()
%     computes the inverse 2 dimensional dost of a given matrix (ex: image)
%     input: matrix of size 2^k
%     output: idost2 coefficients of the matrix
%
%   ST.dcst()
%     computes the dcst coefficients of a given signal
%     input: vector (or matrix) of size 2^k
%     output: the dcst coefficients of the vector (of the columns of the matrix)
%
%   ST.rearrangeDcst()
%     rearranges the dcst coefficients to make them more readable
%     input: a vector of dost coefficients (length n=2^k)
%     output: a n x n matrix which describes the time frequency plane using the
%     dcst coefficients. See ST_algorithms.pdf for a graphical explanation
%
%   ST.idcst()
%     computes the inverse DCST, i.e. given a vector of DCST coefficients gives
%     back the original signal
%     input: vector (or matrix) of size 2^k
%     output: the idcst coefficients of the vector (of the columns of the matrix)
%
%   ST.dcst2()
%     computes the 2 dimensional dcst of a given matrix (ex: image)
%     taking the dcst on the columns and then on the rows
%     input: matrix of size 2^k
%     output: dcst2 coefficients of the matrix
%
%   ST.idcst2()
%     computes the inverse 2 dimensional dcst of a given matrix (ex: image)
%     input: matrix of size 2^k
%     output: dcst2 coefficients of the matrix
%
%   ST.fourier(), ST.ifourier(), ST.fourier2(), ST.ifourier2()
%     these functions compute the normalized and centered FFT the 1 and 2
%     dimensional cases and their inverses
%
%                             PRIVATE FUNCTIONS
%   ST.dostbw() and ST.dcstbw()
%     give the frequency bandwidth decomposition of the dost and dcst resepectively
%     input: number of samples (a power of 2)
%     output: the bandwidth decomposition.
%     examples: ST.dostbw(16) = [1 4 2 1 1 1 2 4], ST.dcstbw(16) = [1 1 2 4 8]
%
% Additional details:
% Copyright (c) by U. Battisti and L. Riba
% $Revision: 1.0 $
% $Date: 7th November 2015$




classdef ST
    methods(Static = true, Access = public)
        
        % normalized and centered fft
        function out = fourier(in)
            switch isvector(in)
                case true
                    out = (1 ./ sqrt(length(in))) .* fftshift(fft(ifftshift(in)));
                otherwise
                    out = (1 ./ sqrt(size(in, 1))) .* fftshift(fft(ifftshift(in, 1), [], 1), 1);
            end
        end
        
        % normalized and centered ifft
        function out = ifourier(in)
            switch isvector(in)
                case true
                    out = sqrt(length(in)) .* fftshift(ifft(ifftshift(in)));
                otherwise
                    out = sqrt(size(in, 1)) .* fftshift(ifft(ifftshift(in, 1), [], 1), 1);
            end
        end
        
        % fourier2 transform
        function out = fourier2(in)
            switch (isvector(in))
                case true
                    out = ST.fourier(in);
                otherwise
                    out = ST.fourier(ST.fourier(in).').';
            end
        end
        
        % ifourier2 transform
        function out = ifourier2(in)
            switch (isvector(in))
                case true
                    out = ST.ifourier(in);
                otherwise
                    out = ST.ifourier(ST.ifourier(in).').';
            end
        end
        
        % DOST - discrete orthonormal Stockwell transforms
        % dost transform
        function out = dost(in)
            switch(isvector(in))
                case true
                    out = ST.fourier(in);
                    D = length(in);
                    bw = ST.dostbw(D);
                    k = 1;
                    for ii = bw
                        if ii == 1
                            k = k + ii;
                        else
                            out(k : k + ii - 1) = ST.ifourier(out(k : k + ii - 1));
                            k = k + ii;
                        end
                    end
                otherwise
                    out = ST.fourier(in);
                    M = size(in, 1);
                    bw = ST.dostbw(M);
                    k = 1;
                    for ii = bw
                        if ii == 1
                            k = k + ii;
                        else
                            out(k : k + ii - 1, :) = ST.ifourier(out(k:k + ii - 1, :));
                            k = k + ii;
                        end
                    end
            end
        end
        
        %idost transform
        function out = idost(in)
            switch (isvector(in))
                case true
                    out = in;
                    D = length(in);
                    bw = ST.dostbw(D);
                    k = 1;
                    for ii = bw
                        if ii == 1
                            k = k + ii;
                        else
                            out(k : k + ii - 1) = ST.fourier(out(k : k + ii - 1));
                            k = k + ii;
                        end
                    end
                    out = ST.ifourier(out);
                otherwise
                    out = in;
                    M = size(in, 1);
                    bw = ST.dostbw(M);
                    k = 1;
                    for ii = bw
                        if ii == 1
                            k = k + ii;
                        else
                            out(k : k + ii - 1, :) = ST.fourier(out(k : k + ii - 1, :));
                            k = k + ii;
                        end
                    end
                    out = ST.ifourier(out);
            end
            
        end
        
        % dost2 transform
        function out = dost2(in)
            switch (isvector(in))
                case true
                    out = ST.dost(in);
                otherwise
                    out = ST.dost(ST.dost(in).').';
            end
        end
        
        % idost2 transform
        function out = idost2(in)
            switch (isvector(in))
                case true
                    out = ST.idost(in);
                otherwise
                    out = ST.idost(ST.idost(in).').';
            end
        end
        
        % DCST - discrete cosine Stockwell transforms
        % dcst transform
        function out = dcst(in)
            switch(isvector(in))
                case true
                    out = dct(in);
                    D = length(in);
                    bw = ST.dcstbw(D);
                    k = 1;
                    for ii = bw
                        if ii == 1
                            k = k + ii;
                        else
                            out(k : k + ii - 1) = idct(out(k : k + ii - 1));
                            k = k + ii;
                        end
                    end
                otherwise
                    out = dct(in);
                    M = size(in, 1);
                    bw = ST.dcstbw(M);
                    k = 1;
                    for ii = bw
                        if ii == 1
                            k = k + ii;
                        else
                            
                            out(k : k + ii - 1, :) = idct(out(k:k + ii - 1, :));
                            k = k + ii;
                        end
                    end
            end
        end
        
        %idcst transform
        function out = idcst(in)
            switch (isvector(in))
                case true
                    out = in;
                    D = length(in);
                    bw = ST.dcstbw(D);
                    k = 1;
                    for ii = bw
                        if ii == 1
                            k = k + ii;
                        else
                            out(k : k + ii - 1) = dct(out(k : k + ii - 1));
                            k = k + ii;
                        end
                    end
                    out = idct(out);
                otherwise
                    out = in;
                    M = size(in, 1);
                    bw = ST.dcstbw(M);
                    k = 1;
                    for ii = bw
                        if ii == 1
                            k = k + ii;
                        else
                            out(k : k + ii - 1, :) = dct(out(k : k + ii - 1, :));
                            k = k + ii;
                        end
                    end
                    out = idct(out);
            end
            
        end
        
        % dcst2 transform
        function out = dcst2(in)
            switch (isvector(in))
                case true
                    out = ST.dcst(in);
                otherwise
                    out = ST.dcst(ST.dcst(in).').';
            end
        end
        
        % idcst2 transform
        function out = idcst2(in)
            switch (isvector(in))
                case true
                    out = ST.idcst(in);
                otherwise
                    out = ST.idcst(ST.idcst(in).').';
            end
        end
        
        % rearrangeDost
        % the dost coefficients have a meaning as time-frequency energy
        % content. This function rearrange the coefficients in a matrix form
        % which encodes these information
        function out = rearrangeDost(in)
            switch(isvector(in))
                case true
                    D = length(in);
                    bw = ST.dostbw(D);
                    tbw = D ./ bw;
                    out = zeros(D, D);
                    count = 1;
                    for hh = 1:length(bw)
                        for kk = 1: bw(hh)
                            ii = D + 1 - sum(bw(1:hh));
                            jj = tbw(hh) * (kk - 1) + 1;
                            tmp = repmat(count, bw(hh), tbw(hh));
                            out(ii:ii + size(tmp, 1) - 1, jj:jj + size(tmp, 2) - 1) = tmp;
                            count = count + 1;
                        end
                    end
                    out = in(out);
                otherwise
                    error('ERROR: This method is intended to work with 1d vectors');
            end
        end
        
        % rearrangeDcst
        % the dcst coefficients have a meaning as time-frequency energy
        % content. This function rearrange the coefficients in a matrix form
        % which encodes these information
        function out = rearrangeDcst(in)
            switch (isvector(in))
                case true
                    D = length(in);
                    bw = ST.dcstbw(D);
                    tbw = D ./ bw;
                    out = zeros(D, D);
                    count = 1;
                    for hh = 1:length(bw)
                        for kk = 1: bw(hh)
                            ii = D + 1 - sum(bw(1:hh));
                            jj = tbw(hh) * (kk - 1) + 1;
                            tmp = repmat(count, bw(hh), tbw(hh));
                            out(ii:ii + size(tmp, 1) - 1, jj:jj + size(tmp, 2) - 1) = tmp;
                            count = count + 1;
                        end
                    end
                    out = in(out);
                otherwise
                    error('ERROR: This method is intended to work with 1d vectors');
            end
        end
    end
    
    methods(Static = true, Access = private)
        % bandwidth partitioning for the dost, i.e. dyadic negative to positive
        % it gives the size of the dost bandwidths
        % dostbw(16) = [1, 4, 2, 1, 1, 1, 2, 4]
        function out = dostbw(in)
            out = 2 .^ ([0, log2(in) - 2: - 1:0, 0, 0:log2(in) - 2]);
        end
        
        % bandwidth partitioning for the dost, i.e. dyadic zero to positive
        % it gives the size of the dcst bandwidths
        % dcstbw(16) = [1, 1, 2, 4, 8]
        function out = dcstbw(in)
            out = 2 .^ ([0, 0:log2(in) - 1]);
        end
        
    end
    
end
