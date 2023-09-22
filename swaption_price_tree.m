function price = swaption_price_tree(year_steps, discounts_dates, discounts, sigma, a, t_alpha, t_omega, coupon_payment_dates, flag)

% Default : swaption equivalent to a put on a coupon bond with strike = 1
% flag = -1: swaption equivalent to a call on a coupon bond with strike = 1
if nargin < 9
    flag = 1;
end

% yearfrac modes
thirty360 = 6;
act365 = 3;

% Grid parameters
delta_t = 1/year_steps;
settlement_date = discounts_dates(1);
M = floor(year_steps * yearfrac(settlement_date, t_alpha, 3) + 1);
sigma_hat = sigma * sqrt((1-exp(-2*a*delta_t))/2/a);
mu = 1 - exp(-a*delta_t);
delta_x = sqrt(3)*sigma_hat;
l_max = floor((1-sqrt(2/3))/mu+1);
l_min = -l_max;
l_vec = (l_max:-1:l_min)';
x = l_vec * delta_x;

% Transition probabilities
pA = [0.5*(1/3-l_vec*mu+(l_vec*mu).^2), 2/3-(l_vec*mu).^2, 0.5*(1/3+l_vec*mu+(l_vec*mu).^2)];
pB = [1/2*(1/3+l_min*mu+(l_min*mu).^2), -1/3-2*l_min*mu-(l_min*mu).^2,  1/2*(7/3+3*l_min*mu+(l_min*mu).^2)];
pC = [1/2*(7/3-3*l_max*mu+(l_max*mu).^2), -1/3+2*l_max*mu-(l_max*mu).^2, 1/2*(1/3-l_max*mu+(l_max*mu).^2)];

p = pA;
p(1,:) = pC;
p(end,:) = pB;

% Discount factors
B_alpha = Disc_interp(discounts, discounts_dates, t_alpha);
B_zero_alpha_omega = Disc_interp(discounts, discounts_dates, t_omega)/B_alpha;
B_zero_alpha_i = Disc_interp(discounts, discounts_dates, coupon_payment_dates)/B_alpha;

% Compute BPV
delta_alpha_omega = yearfrac([t_alpha; coupon_payment_dates(1:end-1)], coupon_payment_dates, thirty360);
BPV = delta_alpha_omega' * B_zero_alpha_i; 


%% INITIALIZATION
% ATM swaption so we take strike equal to swap rate
strike = (1-B_zero_alpha_omega)/BPV;

coupons = strike * ones(length(coupon_payment_dates),1); 
coupons(end) = coupons(end)+1;

sigma_function = @(t1,t2) sigma/a * (1-exp(-a*(t2-t1)));

strike_put = 1;
cbp = zeros(length(x), 1);

for i = 1:length(x)
[cbp(i), ~] = coupon_bond_price(settlement_date, t_alpha, coupon_payment_dates, coupons, discounts, discounts_dates, sigma, a, x(i));
end

payoff = max(flag*(strike_put - cbp), 0);

Price = zeros(length(x),M);
Price(:,end) = payoff;

t_grid = round(settlement_date:delta_t*365:t_alpha)';
delta_zero_alpha = yearfrac(settlement_date, t_grid, act365);

sigma_hat_star = sigma/a*sqrt(delta_t-2*(1-exp(-a*delta_t))/a+(1-exp(-2*a*delta_t))/(2*a));


%% TREE

% TO VISUALIZE THE X TREE:

% % Initialization of the tree
% x_tree = zeros(length(x),M);
% 
% % Building the x tree
% for i = 1:floor(length(x)/2)
%     x_tree(i,l_vec(i)+1:end) = x(i);
%     x_tree(end-i+1,l_vec(i)+1:end) = x(end-i+1);
% end

% Computing the price with backward discounting
for j=M-1:-1:1
    
    % Computing forward discounts B(t0,ti,ti+1)
    B0 = Disc_interp( discounts, discounts_dates, t_grid(j+1)) / Disc_interp(discounts, discounts_dates, t_grid(j));
    
    % Computing discounts B(ti,ti+1) = B(ti,ti,ti+1)
    B = @(x) B0.*exp(-x/sigma*sigma_function(0,delta_t)+ -0.5 * integral(@(u) (sigma_function(u,delta_zero_alpha(j)+delta_t).^2 - sigma_function(u,delta_zero_alpha(j)).^2), 0, delta_zero_alpha(j), 'ArrayValued', true));  
    
    % Computing the price
    if j > l_max  % horizontal part of the tree
        for i = 1:length(x)
            if i==1
                dx = [0, -delta_x, -2*delta_x]';
                indices = i:i+2;
            elseif i==length(x)
                dx = [2*delta_x, delta_x, 0]';
                indices = i-2:i;
            else
                dx = [delta_x,0,-delta_x]';
                indices = i-1:i+1;
            end
            D = B(x(indices)) .* exp(-0.5*sigma_hat_star^2-sigma_hat_star/sigma_hat*((exp(-a*delta_t)*dx+mu*x(indices))));
            Price(i,j) = p(i,:) * (Price(indices,j+1).*D);
        end
    else % triangular part of the tree
        for i=(l_max-j)+2:length(x)-((l_max-j)+1)
            dx = [delta_x,0,-delta_x]';
            D = B(x(i-1:i+1)).*exp(-0.5*sigma_hat_star^2-sigma_hat_star/sigma_hat*((exp(-a*delta_t)*dx+mu*x(i-1:i+1))));
            Price(i,j) = p(i,:) * (Price(i-1:i+1,j+1).*D);
        end
    end
end

% Plot sparse price
spy(sparse(Price))

% Return price
price = Price(l_max+1,1);

end