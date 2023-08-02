kT=0.001987*298;

Emid=0.00125:0.0025:1;

G=load('G1_6.dat');

plot(Emid,kT*G,'LineWidth',1)
xlabel('\lambda')
ylabel('Free Energy [kcal/mol]')
title('Implicit Constraint Free Energy Ns=6')
set(1,'PaperSize',[3.25,2.5])
set(1,'PaperPosition',[0,0,3.25,2.5])
saveas(1,'G1_6.pdf')

