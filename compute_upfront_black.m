function upfront = compute_upfront_black(party_B, party_A, d, sigma, S0, strike, dates, discounts, date_start, N)
% Computes the upfront with Black model
%
% INPUT
% alpha:        alpha of the MVM model
% party_B:      party_B struct with dates and discounts
% party_A:      party_A struct with dates and discounts
% d:            discount yeld of the underlying
% sigma:        sigma parameter of the black model
% S0:           value of the underlying at 0
% strike:       strike of the condition of the 1y coupon
% discounts:    discounts from the bootstrap
% dates:        dates of the discounts from the bootstrap
% date_start:   date when the contract starts
% N:            number of MC simulations


% Parameters for the simulation
rng(42)
act365 = 3;
ttm = party_B.ttm_B(1:end-1);
n = length(ttm);

% Setting delta times and discounts
delta_B = yearfrac([date_start; party_B.payment_dates_B(1:end-1)], party_B.payment_dates_B, act365);
discounts_B_reset = Disc_interp(discounts,dates,party_B.dates_reset_B);
forward_discounts = discounts_B_reset(1:end-1)./[1; discounts_B_reset(1:end-2)];
forward_rates = -log(forward_discounts)./party_B.delta_times_B(1:end-1);

% Initializing variables
St = zeros(N,n+1);
St(:,1) = S0;

% Simulating St
rng(42)
z = randn(N,1);
St(:,2) = St(:,1).*exp( (forward_rates - d - sigma^2/2)*party_B.delta_times_B(1) + sqrt(party_B.delta_times_B(1))*sigma*z  );

% Computing probabilities 
P = (St(:,2) < strike);

% Computing the NPVs (without considering the upfront)
NPV_B = party_B.coupon_rates(1)*delta_B(1)*party_B.discounts_B(1) * P + ...
        party_B.coupon_rates(2)*delta_B(2)*party_B.discounts_B(2) * (1-P);

NPV_A = (1-party_B.discounts_B(1) + party_A.spol_A*party_A.delta_times_A(1:4)'*party_A.discounts_A(1:4)) * P + ...
        (1-party_B.discounts_B(2) + party_A.spol_A*party_A.delta_times_A(1:8)'*party_A.discounts_A(1:8)) * (1-P);

% Computing the upfront
[upfront,~,ci] = normfit(NPV_A - NPV_B);
fprintf('Length of the confidence interval: %.4f bp \n', (ci(2)-ci(1))*1e4)

end