function k= gaussian_filter(X, rate, center, sd)
    N = size(X,2);
    d=1./rate;
    freq=[0:ceil(N/2-1), ceil(-(N)/2):1:-1]/(d*N);
    k = exp((-(abs(freq) - center).^2)/(2 * (sd^2)));
    k=k/sum(k);
end 