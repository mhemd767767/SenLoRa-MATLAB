function [kappa, detail] = senlora_sample_entropy(signal, opts)
%SENLORA_SAMPLE_ENTROPY Sample entropy from SenLoRa Eq. (9), Appendix (18-22).
%   [KAPPA, DETAIL] = SENLORA_SAMPLE_ENTROPY(SIGNAL, OPTS)
%   counts distinct unordered pairs with Chebyshev distance < r, excluding
%   self-matches. The ratio of pair counts is identical to the ratio of
%   ordered-pair counts in the paper because both are multiplied by two.
%
%   OPTS.m (default 2) and OPTS.r_fraction (default 0.2) are researcher
%   choices: THE PAPER DOES NOT SPECIFY THESE PARAMETER VALUES.
%   Set OPTS.r_absolute to specify a fixed tolerance. In particular,
%   the paper's gamma=0.01496 CANNOT be transferred blindly to these
%   numerical defaults or to a different CSI time-series representation.
%
%   No match pairs => undefined (NaN) or infinite (Inf) entropy. Never
%   treat insufficient matches as an entropy indicating adequate sensing.
if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'m'), opts.m = 2; end
if ~isfield(opts,'r_fraction'), opts.r_fraction = 0.2; end
validateattributes(opts.m,{'numeric'},{'scalar','integer','>=',1});
validateattributes(opts.r_fraction,{'numeric'},{'scalar','real','positive'});
x = double(signal(:));
if ~isreal(x), error('Sample entropy input must be a real-valued CSI feature.'); end
if any(~isfinite(x)), error('Sample entropy input cannot contain NaN or Inf.'); end
n = numel(x); m = double(opts.m);
detail = struct('count_m',0,'count_m1',0,'n',n,'m',m, ...
    'r',NaN,'valid',false,'reason','insufficient_samples');
kappa = NaN;
if n < m+3, return; end
if isfield(opts,'r_absolute') && ~isempty(opts.r_absolute)
    r = double(opts.r_absolute);
    validateattributes(r,{'numeric'},{'scalar','real','positive','finite'});
else
    r = max(opts.r_fraction*std(x), ...
        100*eps(max(1,max(abs(x)))));
end
detail.r = r;
% The paper counts length-m matches for all N-m+1 templates, and length
% (m+1) matches for N-m templates; follow those exact ranges.
b = 0; a = 0;
for i = 1:(n-m+1)
    for j = (i+1):(n-m+1)
        if max(abs(x(i:i+m-1)-x(j:j+m-1))) < r
            b = b+1;
        end
    end
end
for i = 1:(n-m)
    for j = (i+1):(n-m)
        if max(abs(x(i:i+m)-x(j:j+m))) < r
            a = a+1;
        end
    end
end
detail.count_m = b;
detail.count_m1 = a;
if b == 0
    detail.reason = 'no_m_matches';
elseif a == 0
    kappa = Inf;
    detail.reason = 'no_m_plus_one_matches';
else
    kappa = -log(a/b);
    detail.valid = isfinite(kappa);
    detail.reason = 'ok';
end
end
