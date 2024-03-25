%Frh=input("Ange radialkrafter på högerlagret "); Tas bort då kraften antas
%vara mindre än axial och radial krafter på vänstra lagret. 
Fr=input("Ange radialkrafter på vänsterlagret ");
Fa=input("Ange axialkrafter på vänsterlagret ");
clearvars -except Fr Fa
clc

f0=input('Ange lagrets f0 värde(fås från SKF P.32) ');
C0=input("Ange lagrets C0 värde (fås från SKF) ");

jamf=(f0*Fa)/C0;
lager=[0.172,0.29,0.46,1.88; 
0.345,0.32,0.46,1.71;
0.689,0.36,0.46,1.52;
1.03,0.38,0.46,1.41;
1.38,0.40,0.46,1.34;
2.07,0.44,0.46,1.23;
3.45,0.49,0.46,1.10;
5.17,0.54,0.46,1.01;
6.89,0.54,0.46,1.00];
n=lager(:,1);
y1=zeros(n,1);

for i = 1:Length(n) 

    if n(i)<=jamf
        y2=n(i);
        it=i;

    elseif n(i)>=jamf
        y1(i)=n(i);
    end
end
disp y1

y1=input('Ange minsta värdet på y1 ');

