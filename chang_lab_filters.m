% chang lab frequencies
function get_filter_param_chang_lab(fq_min,fq_max)
if nargin<2
    fq_min = 4.0749286538265;
    fq_max = 200.;
end
scale = 7.;
cfs = 2.^((log2(fq_min) * scale:1: log2(fq_max) * scale) / scale);
sds = 10.^( log10(.39) + .5 * (log10(cfs)));
sds = (sds) * sqrt(2.);
end

% standard neuro bands 
%bands = ['theta', 'alpha', 'beta', 'gamma', 'high gamma']
%min_freqs = [4., 8., 15., 30., 70.]
%max_freqs = [7., 14., 29., 59., 150.]
%HG_freq = 200.
