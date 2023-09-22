function void_plot_error(year_steps,error_1, error_2)



figure
semilogy(year_steps, error_1, '-o', LineWidth = 3)
xlabel('year step', fontsize = 18) 
ylabel('error', fontsize = 18)
hold on 
semilogy( year_steps, error_2,'-or', LineWidth = 3)
title('Error on the two swaption in semilogy scale', fontsize = 15)
legend('error 1st swaption', 'error 2nd swaption', fontsize = 18)



end