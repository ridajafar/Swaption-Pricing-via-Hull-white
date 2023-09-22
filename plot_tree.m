function plot_tree(year_steps, discount_dates, discounts, sigma, a, T_alpha_1, T_alpha_2, T_omega, Coupon_payment_dates_1, Coupon_payment_dates_2, swaption_price_1_jamshidian, swaption_price_2_jamshidian)
% Plots price of the swaption via tree converging to the jamshidian

swaption_price_1_tree = zeros(length(year_steps),1);
swaption_price_2_tree = zeros(length(year_steps),1);
error_1 = zeros(length(year_steps),1);
error_2 = zeros(length(year_steps),1);

for i = 1:length(year_steps)
    swaption_price_1_tree(i) = swaption_price_tree(year_steps(i), discount_dates, discounts, sigma, a, T_alpha_1 + 1, T_omega, Coupon_payment_dates_1);
    swaption_price_2_tree(i) = swaption_price_tree(year_steps(i), discount_dates, discounts, sigma, a, T_alpha_2 + 1, T_omega, Coupon_payment_dates_2);
    error_1(i) = abs(swaption_price_1_tree(i) - swaption_price_1_jamshidian);
    error_2(i) = abs(swaption_price_2_tree(i) - swaption_price_2_jamshidian);
end

void_plot(year_steps,swaption_price_1_tree, swaption_price_1_jamshidian, 1)
void_plot(year_steps,swaption_price_2_tree, swaption_price_2_jamshidian, 2)
void_plot_error(year_steps,error_1, error_2)

end