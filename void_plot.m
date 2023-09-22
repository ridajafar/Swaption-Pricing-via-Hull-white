function void_plot(year_steps,swaption_price_1_tree, swaption_price_1_jamshidian, i)

figure
plot(year_steps, swaption_price_1_tree,'-ob', LineWidth=3)
xlabel('year step') 
ylabel('price')
hold on
plot(year_steps,ones(length(year_steps), 1)*swaption_price_1_jamshidian,'-or', LineWidth= 3)
xlabel('year step', fontsize = 18) 
ylabel('price', fontsize = 18)
legend('tree price', 'jamshidian price', fontsize = 18)
title(['Price behaviour with different methods swaption ', num2str(i)], fontsize = 15)



end
