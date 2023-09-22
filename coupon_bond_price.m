function [cbp,B_t_i]  = coupon_bond_price(set_date, t_alpha, CouponPaymentDates, CouponValue, discounts, discounts_dates, sigma, a, x)
% Computes coupon bond price

coupons_ttm= yearfrac(t_alpha, CouponPaymentDates, 6); %% t_alpha --- t_i with i belonging (alpha, w)
total_ttm = yearfrac(set_date, CouponPaymentDates, 6); %% from t_0 ---- t_i with i belonging (alpha, w)
ttm_alpha = yearfrac(set_date, t_alpha, 6); %% t0 ----- t_alpha

deltas = yearfrac([t_alpha; CouponPaymentDates(1:end-1)], CouponPaymentDates, 6); %% delta between coupon payment dates starting from t_alpha
zRates = zeroRates(discounts_dates, discounts);

zRate_alpha = interp1(discounts_dates(2:end), zRates, t_alpha, 'linear', 'extrap');
B_alpha = exp(- yearfrac(set_date, t_alpha, 3).* zRate_alpha / 100);

zRate_zero = interp1(discounts_dates(3:end), zRates(2:end), CouponPaymentDates, 'linear', 'extrap'); %% extracting wanted value for coupon payments dates
B_0_i =  exp(-yearfrac(set_date, CouponPaymentDates, 3).* zRate_zero / 100)/B_alpha; %% forward discount between 0 and i 


sigma_HW = @(tzero, T) sigma/a * (1 - exp(-a*(T - tzero)));
integral_B = integral(@(u) (sigma_HW(u, total_ttm).^2 - sigma_HW(u, ttm_alpha).^2), 0, ttm_alpha, 'ArrayValued', true); 

B_t_i = B_0_i.*exp(-x/sigma.*sigma_HW(0, coupons_ttm) - 0.5*integral_B);

cbp = sum(CouponValue.*B_t_i.*deltas);

end
