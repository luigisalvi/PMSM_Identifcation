function fit_plot(d,vel_fit,torque_fit)
t=1:d;
figure(1)
stem(t,vel_fit,'b')
xticks(1:d), xlabel('Datasets'), ylabel('% of fitting'), title('Fitting su Velocità')

figure(2)
stem(t,torque_fit,'b')
xticks(1:d), xlabel('Datasets'), ylabel('% of fitting'), title('Fitting su Coppia')

end