clear variables
close all

home = pwd;

% load data for S1 (main fault segment, Gofar C)
S1 = load('Gofar_3seg_C_temp_stress_press_1300_GSRM.txt');
loop_size = size(S1,1)/221;

% load data for S1, the secondary fault segment, Gofar B
S2 = load('Gofar_3seg_B_temp_stress_press_1300_GSRM.txt');
loop_size2 = size(S2,1)/221;

% load data for R1, western ridge segment
R1 = load('Gofar_3seg_R1_temp_stress_press_1300_GSRM.txt');
loop_size3 = size(R1,1)/221;

% load data for O1, the offset between Gofar C and Gofar B
O1 = load('Gofar_3seg_O1_temp_stress_press_1300_GSRM.txt');
loop_size4 = size(O1,1)/221;

z_interp = -100000:100:0;

for ii=1:loop_size
    temp = sortrows(S1((ii*221-220):ii*221,:),3);
    S1_temp_int(:,ii) = interp1(temp(:,3),temp(:,4),z_interp,'pchip');
    S1_stress_int(:,ii) = interp1(temp(:,3),temp(:,5),z_interp,'pchip');
    S1_press_int(:,ii) = interp1(temp(:,3),temp(:,6),z_interp,'pchip');
    X(ii)=temp(1,1);
end
z(:,1) = temp(:,3);

indx0 = find(X == (min(abs(X))));

% Gofar C geotherm for center of transform (depth v. temp)
Gofar_geotherm = [S1_temp_int(:,indx0),z_interp(1,:)'];


for ii=1:loop_size2
    tempS2 = sortrows(S2((ii*221-220):ii*221,:),3);
    S2_temp_int(:,ii) = interp1(tempS2(:,3),tempS2(:,4),z_interp,'pchip');
    S2_stress_int(:,ii) = interp1(tempS2(:,3),tempS2(:,5),z_interp,'pchip');
    S2_press_int(:,ii) = interp1(tempS2(:,3),tempS2(:,6),z_interp,'pchip');
    X2(ii)=tempS2(1,1);
end
z(:,1) = temp(:,3);

R1 = sortrows(R1,2);

for ii=1:loop_size3
   tempR1 = sortrows(R1((ii*221-220):ii*221,:),3);
   R1_temp_int(:,ii) = interp1(tempR1(:,3),tempR1(:,4),z_interp,'pchip');
   R1_stress_int(:,ii) = interp1(tempR1(:,3),tempR1(:,5),z_interp,'pchip');
   R1_press_int(:,ii) = interp1(tempR1(:,3),tempR1(:,6),z_interp,'pchip');
   R1Y(ii)=tempR1(1,2);
end

zR1(:,1) = tempR1(:,3);

O1 = sortrows(O1,2);

for ii=1:loop_size4
   tempO1 = sortrows(O1((ii*221-220):ii*221,:),3);
   O1_temp_int(:,ii) = interp1(tempO1(:,3),tempO1(:,4),z_interp,'pchip');
   O1_stress_int(:,ii) = interp1(tempO1(:,3),tempO1(:,5),z_interp,'pchip');
   O1_press_int(:,ii) = interp1(tempO1(:,3),tempO1(:,6),z_interp,'pchip');
   O1Y(ii)=tempO1(1,2);
end
 zO1(:,1) = tempO1(:,3);
%%
close all
clim = [0 1400]
colormap(jet(300))


subplot(2,1,2)
set(gcf,'PaperOrientation','landscape','papertype','uslegal');
set(gcf,'position',[100 250 890 495],'paperposition',[0.500    0.5000    12.5000    8.0000])

imagesc(X/1e3 + 47.5,z_interp/1e3,S1_temp_int,clim)
axis([-50 140 -30 0])
set(gca,'YDir','normal')
hold on

imagesc((R1Y/1e3),z_interp/1e3,R1_temp_int,clim)

hold on
imagesc(O1Y/1e3 + 95,z_interp/1e3,O1_temp_int,clim)

hold on
imagesc(X2/1e3 + 62.5,z_interp/1e3,S2_temp_int,clim)

set(gca,'units','points')
set(gcf,'position',[100 250 850 500])
set(gca,'position',[60 45  755  200])

v1 = vline(-0,'k')
v2 = vline(95,'k')
v3 = vline(110,'k')

set(gca,'fontsize',13)
xlabel('Distance (km)')
ylabel('Depth (km)')

c=colorbar;
set(c,'ydir','reverse','fontsize',13)
set(get(c,'title'),'String', 'Temp ^{\circ}C','fontsize',13);


subplot(2,1,1)
set(gca,'units','points')
set(gca,'position',[60 280  700  200])

% seg 1 = 95 m, seg 2 = 30, seg 3 = 46
p1 = plot([0 95],[0 0],'b','linewidth',2)
hold on
p2 = plot([95 125],[15 15],'b','linewidth',2)
p3 = plot([0 0],[0 -50],'r','linewidth',2)
p4 = plot([95 95],[0 15],'r','linewidth',2)
p5 = plot([125 125],[15 25],'--r','linewidth',2)
p6 = plot([125 140],[25 25],'color',[0 0.4470 0.7410],'linestyle','--','linewidth',2)
text(40,5,'Segment C','fontsize',13)
text(100,20,'Segment B','fontsize',13)
text(-2,-30,'EPR','fontsize',13,'rotation',90)
text(93,2,'ITSC','fontsize',13,'rotation',90)

axis([-50 140 -50 50])

set(gca,'fontsize',13)
ylabel('Across Track Distance (km)')

saveas(gcf,'Gofar_temp_structure_1300mpt','pdf')

figure(2)
plot(Gofar_geotherm(:,1),Gofar_geotherm(:,2))
set(gca,'fontsize',13)
set(gca,'XAxisLocation','top','YAxisLocation','left');
xlabel('Temperature (Deg C)')
ylabel('Depth (km)')
saveas(gcf,'Gofar_Geotherm_1300mpt','pdf')

