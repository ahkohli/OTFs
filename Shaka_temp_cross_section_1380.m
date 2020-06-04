clear variables
close all

home = pwd;

% Load data along the transform domain
S1 = load('Shaka_S1_temp_stress_press_strain_GSRM.txt');
loop_size = size(S1,1)/221;

% Load data from the western ridge segment
R1 = load('Shaka_R1_temp_stress_press_strain_GSRM.txt');
loop_size2 = size(R1,1)/221;

% Load data from eastern ridge segment
R2 = load('Shaka_R2_temp_stress_press_strain_GSRM.txt');

z_interp = -100000:100:0;

for ii=1:loop_size
   temp = sortrows(S1((ii*221-220):ii*221,:),3);
   S1_temp_int(:,ii) = interp1(temp(:,3),temp(:,4),z_interp,'pchip');
   S1_stress_int(:,ii) = interp1(temp(:,3),temp(:,5),z_interp,'pchip');
   S1_press_int(:,ii) = interp1(temp(:,3),temp(:,6),z_interp,'pchip');
   S1_strain_int(:,ii) = interp1(temp(:,3),temp(:,7),z_interp,'pchip');
   X(ii)=temp(1,1);
end
   z(:,1) = temp(:,3);
   
   indx0 = find(X == (min(abs(X))));
   
   % Shaka geotherm for center of transform (depth v. temp)
   Shaka_geotherm = [S1_temp_int(:,indx0),z_interp(1,:)'];
   
    % Shaka shear strain rate for center of transform (depth v. strain)
   Shaka_strain = [S1_strain_int(:,indx0),z_interp(1,:)'];

R1 = sortrows(R1,2);

for ii=1:loop_size2
   tempR1 = sortrows(R1((ii*221-220):ii*221,:),3);
   R1_temp_int(:,ii) = interp1(tempR1(:,3),tempR1(:,4),z_interp,'pchip');
   R1_stress_int(:,ii) = interp1(tempR1(:,3),tempR1(:,5),z_interp,'pchip');
   R1_press_int(:,ii) = interp1(tempR1(:,3),tempR1(:,6),z_interp,'pchip');
   R1_strain_int(:,ii) = interp1(tempR1(:,3),tempR1(:,7),z_interp,'pchip');
   R1Y(ii)=tempR1(1,2);
end

 zR1(:,1) = tempR1(:,3);


R2 = sortrows(R2,2);

% note that ridge segments are same size, so same number of elements
for ii=1:loop_size2
   tempR2 = sortrows(R2((ii*221-220):ii*221,:),3);
   R2_temp_int(:,ii) = interp1(tempR2(:,3),tempR2(:,4),z_interp,'pchip');
   R2_stress_int(:,ii) = interp1(tempR2(:,3),tempR2(:,5),z_interp,'pchip');
   R2_press_int(:,ii) = interp1(tempR2(:,3),tempR2(:,6),z_interp,'pchip');
   R2_strain_int(:,ii) = interp1(tempR2(:,3),tempR2(:,7),z_interp,'pchip');
   R2Y(ii)=tempR2(1,2);
end
 zR2(:,1) = tempR2(:,3);

close all
clim = [0 1200]
colormap(jet(300))


subplot(2,1,2)
set(gcf,'PaperOrientation','landscape','papertype','uslegal');
set(gcf,'position',[100 250 890 495],'paperposition',[0.500    0.5000    12.5000    8.0000])

imagesc(X/1e3 + 100,z_interp/1e3,S1_temp_int,clim)
axis([-50 250 -80 0])
set(gca,'YDir','normal')
hold on

imagesc((-1*R1Y/1e3),z_interp/1e3,R1_temp_int,clim)

hold on
imagesc(-1*R2Y/1e3 + 200,z_interp/1e3,R2_temp_int,clim)

set(gca,'units','points')
set(gcf,'position',[100 250 850 500])
set(gca,'position',[60 45  755  200])

v1 = vline(-0,'k')
v2 = vline(200,'k')

set(gca,'fontsize',13)
xlabel('Distance (km)')
ylabel('Depth (km)')

c=colorbar;
set(c,'ydir','reverse','fontsize',13)
set(get(c,'title'),'String', 'Temp ^{\circ}C','fontsize',13);


subplot(2,1,1)
set(gca,'units','points')
set(gca,'position',[60 280  700  200])

p1 = plot([0 200],[0 0],'linewidth',2)
hold on
p2 = plot([0 0],[0 50],'r','linewidth',2)
p3 = plot([200 200],[0 -50],'r','linewidth',2)
text(70,5,'Shaka Transform Fault','fontsize',13)
text(-5,20,'SWIR','fontsize',13,'rotation',90)
text(204,-20,'SWIR','fontsize',13,'rotation',270)

axis([-50 250 -50 50])

set(gca,'fontsize',13)
ylabel('Across Track Distance (km)')

saveas(gcf,'Shaka_temp_structure','pdf')

figure(2)
plot(Shaka_geotherm(:,1),Shaka_geotherm(:,2))
set(gca,'fontsize',13)
set(gca,'XAxisLocation','top','YAxisLocation','left');
xlabel('Temperature (Deg C)')
ylabel('Depth (km)')
saveas(gcf,'Shaka_Geotherm','pdf')

% figure(3)
% plot(Shaka_strain(:,1),Shaka_strain(:,2))
% set(gca,'XAxisLocation','top','YAxisLocation','left');
% set(gca,'fontsize',13)
% xlabel('Strain Rate (1/s)')
% ylabel('Depth (km)')
% saveas(gcf,'Shaka_Shear_Strain_2D','pdf')
