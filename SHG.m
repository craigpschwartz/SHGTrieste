
folder_name='experiemntal_folder';

cd (folder_name);
files=dir('*.h5');

i0=zeros(length(files));
area=zeros(length(files));
FWHM=zeros(length(files));

for i=length(files)
i0_par=mean(double(h5read(files(i).name,'/photon_diagnostics/FEL02/I0_monitor/iom_sh_a')));
i0(i)=i0_par;
area_par=mean(double(h5read(files(i).name,'addres of the variable containing the area of the peak '))); 
area(i)=area_par;
prof=double(h5read(files(i).name,'addres of the variable containing the profile of the peak ')); 
M0=max(prof)/2;
FWHM_par=sum(prof>=M0);
FWHM(i)=FWHM_par;
end


figure(1);
subplot(2,1,1);
plot(i0,area,'o'); %to plot normalized value it is sufficient to substitute area with area./i0.^2
xlabel('I0 (\muJ');
ylabel('Signal area');
title('Raw data');
subplot(2,1,2);
plot(i0,FWHM,'o');
xlabel('I0 (\muJ');
ylabel('FWHM peak (pixels)');


























