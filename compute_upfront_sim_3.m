function upfront = compute_upfront_sim_3(alpha, party_B, party_A, sigma, k, eta, F0, strike, date_start, N)
% Computes the upfront via MC simulation for 3 years maturity
%
% INPUT
% alpha:        alpha of the MVM model
% party_B:      party_B struct with dates and discounts
% party_A:      party_A struct with dates and discounts
% sigma:        sigma parameter of the MVM model
% k:            k parameter of the MVM model
% eta:          eta parameter of the MVM model
% F0:           value of the forward at 0
% strike:       strike of the condition of the 1y coupon
% date_start:   date when the contract starts
% N:            number of MC simulations

% Parameters for the simulation
rng(42)
act365 = 3;
n = length(party_B.delta_times_B(1:end-1));
F = zeros(N,n+1);
F(:,1) = F0;


% Setting delta times and discounts
delta_B = yearfrac([date_start; party_B.payment_dates_B(1:end-1)], party_B.payment_dates_B, act365);
delta_times = party_B.delta_times_B(1:end-1);

for i = 1:n

    % Simulating random variables
    g = randn(N,1);
    if alpha == 0
        G = gamrnd(delta_times(i)/k, k/delta_times(i), [N,1]);
    else
        G = random('Stable', alpha, 0, delta_times(i)/k, 0, N, 1);
    end

    % Simulating Ft
    laplace_exp = laplace_exponent(alpha, delta_times(i), k, sigma);
    f_t = sqrt(delta_times(i))*sigma*sqrt(G).*g - (1/2 + eta)*delta_times(i)*sigma^2*G - laplace_exp(eta);
    F(:,i+1) = F(:,i).*exp(f_t);

end

% Computing St
St = F(:,2:end);

% Computing checks 
check_1 = (St(:,1) < strike);
check_2 = (St(:,1) >= strike);
check_3 = (St(:,2) < strike);
check_4 = (St(:,2) >= strike);

% Computing the NPVs (without considering the upfront)
NPV_B = party_B.coupon_rates(1)*delta_B(1)*party_B.discounts_B(1) * check_1 + ...
        party_B.coupon_rates(2)*delta_B(2)*party_B.discounts_B(2) * check_2.*check_3 + ...
        party_B.coupon_rates(3)*delta_B(3)*party_B.discounts_B(3) * check_2.*check_4 ;

NPV_A = (1-party_B.discounts_B(1) + party_A.spol_A*party_A.delta_times_A(1:4)'*party_A.discounts_A(1:4)) * (St(:,1) < strike) + ...
        (1-party_B.discounts_B(2) + party_A.spol_A*party_A.delta_times_A(1:8)'*party_A.discounts_A(1:8)) * ((St(:,1) >= strike).*(St(:,2) < strike)) + ...
        (1-party_B.discounts_B(3) + party_A.spol_A*party_A.delta_times_A(1:end)'*party_A.discounts_A(1:end)) * ((St(:,1) >= strike).*(St(:,2) >= strike));

% Computing the upfront
[upfront,~,ci] = normfit(NPV_A - NPV_B);
fprintf('Length of the confidence interval: %.4f bp \n', (ci(2)-ci(1))*1e4)

end