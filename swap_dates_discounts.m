function [party_A, party_B] = swap_dates_discounts(today, maturity, discounts, dates, frequency_A)
% Sets swap contract dates and discounts both for party A and B
%
% INPUT
% today:        date of today
% maturity:     maturity of the swap
% discounts:    discounts from the bootstrap
% dates:        dates of the discounts from the bootstrap
% frequency_A:  frequency of payments of party_A


% yearfrac format
act360=2;
act365 = 3;
thirty360 = 6;

% Set dates and discounts for B
date_start = busdate(busdate(datenum(today)));
party_B.payment_dates_B = arrayfun(@(t) date_start + t*365, 1:maturity)';
party_B.dates_reset_B = arrayfun(@(t) busdate(busdate(t,-1),-1), party_B.payment_dates_B);
party_B.ttm_B = yearfrac(date_start,party_B.dates_reset_B, thirty360);      % 0-1, 0-2 (reset)
party_B.discounts_B = Disc_interp(discounts,dates,party_B.payment_dates_B);
dates_B = [date_start; party_B.dates_reset_B];
party_B.delta_times_B = yearfrac(dates_B(1:end-1),dates_B(2:end), act365);  % 0-1, 1-2 (reset)

% Set dates and discounts for A
party_A.frequency_A = frequency_A;
payment_dates_A = datetime(2023,02,01):calmonths(12/frequency_A):datetime(2023+maturity,02,01);
party_A.payment_dates_A = datenum(busdate(payment_dates_A))';
discounts_A = Disc_interp(discounts,dates,party_A.payment_dates_A);
party_A.discounts_A = discounts_A(2:end);
party_A.delta_times_A = yearfrac(party_A.payment_dates_A(1:end-1), party_A.payment_dates_A(2:end), act360);

end