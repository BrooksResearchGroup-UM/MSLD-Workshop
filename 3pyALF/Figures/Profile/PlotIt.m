
Emid=0.00125:0.0025:1;

Gphi=Emid;
Gpsi=Emid.*(1-Emid);
Gomega=(1-Emid).*Emid./(0.017+Emid);
Gchi=(1-Emid).*(1-exp(-Emid/0.18));

f=1;
figure(f)
plot(Emid,Gphi,'LineWidth',1)
xlabel('\lambda')
set(f,'PaperSize',[3.25,2.5]/2.5)
set(f,'PaperPosition',[0,0,3.25,2.5]/2.5)
saveas(f,'Gphi.pdf')

f=2;
figure(f)
plot(Emid,Gpsi,'LineWidth',1)
xlabel('\lambda')
set(f,'PaperSize',[3.25,2.5]/2.5)
set(f,'PaperPosition',[0,0,3.25,2.5]/2.5)
saveas(f,'Gpsi.pdf')

f=3;
figure(f)
plot(Emid,Gomega,'LineWidth',1)
xlabel('\lambda')
set(f,'PaperSize',[3.25,2.5]/2.5)
set(f,'PaperPosition',[0,0,3.25,2.5]/2.5)
saveas(f,'Gomega.pdf')

f=4;
figure(f)
plot(Emid,Gchi,'LineWidth',1)
xlabel('\lambda')
set(f,'PaperSize',[3.25,2.5]/2.5)
set(f,'PaperPosition',[0,0,3.25,2.5]/2.5)
saveas(f,'Gchi.pdf')

