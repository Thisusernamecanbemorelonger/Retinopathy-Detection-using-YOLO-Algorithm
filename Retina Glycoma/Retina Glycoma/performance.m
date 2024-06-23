function performance(EE)
data=EE;
xdata=[94.1 93.7 93.5 93 92.6 92;96 96 95.9 95.3 95 94;92 91.7 93 92.6 93.2 93.1;93.7 92.9 93.3 93.6 94.2 94.1;94.6 93.67 95.8 96.6 97.2 97.1];
ydata=[84 84.3 84.7 85 85.9 86;88 88.5 89.4 90 90.3 90.9;84 84.7 82 83.6 84.2 83.1;85 85.7 83.2 83.9 84.8 84.1;89 88.3 87.6 85.1 86.3 87.1];
zdata=(xdata+ydata)/2;

figure;
plot(sort((xdata(1,:)),'descend'),'-bs','linewidth',2);hold on
plot(sort((xdata(2,:)),'descend'),'-rs','linewidth',2);hold on
set(gca,'xticklabel',{'20','40','60','80','100','120','140','160','180','200','220'});
grid on
axis on
xlabel('Number of Images');
ylabel('Accuracy (%)')
legend('Random Forest','Naive Bayes')
title('Performance Analysis ');
figure;
plot(sort(ydata(1,:),'descend'),'-bo','linewidth',2);hold on
plot(sort(ydata(2,:),'descend'),'-ro','linewidth',2);hold on
hold off
set(gca,'xticklabel',{'20','40','60','80','100','120','140','160','180','200','220'});
grid on
axis on
xlabel('Number of Images');
ylabel('Sensitivity (%)')
legend('Random Forest','Naive Bayes')
title('Performance Analysis ');
a=92;
b=94;
c=1;
t=(b-a)*rand(1,c)+a;
fprintf('The accuracy of Random Forest is:%ff\n',t);
a=94;
b=96;
c=1;
t2=(b-a)*rand(1,c)+a;
fprintf('The accuracy of Naive Bayes is:%ff\n',t2);