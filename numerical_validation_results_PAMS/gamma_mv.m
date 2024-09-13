function g = gamma_mv(x, m)
    g = pi^(m*(m-1)/4) * prod(gamma(x - (0:(m-1))/2));
end