function exp = laplace_exponent(alpha, ttm, k, sigma)
% Compute the laplace exponent

switch alpha
    case 0
        exp = @(w) -ttm/k*log(1+k*w*sigma^2);
    otherwise
        if alpha<0 || alpha>1
            disp('Error')
            exp = -1;
        else
            exp = @(w) ttm*(1-alpha)/(k*alpha) * (1-(1+w*k*sigma^2/(1-alpha)).^alpha);
        end
end

end