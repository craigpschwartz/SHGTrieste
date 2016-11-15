%%very basical data extraction software and plotting
folder_name='experiemntal_folder';
out_file='address/where/to/safe/namefile.txt';
out_file_bin='address/where/to/safe/namefilebin.txt';
er=1:0.01:20; %energy bins
i0_bin=nan([length(er) 1]);
area_bin=nan([length(er) 1]);
err_area_bin=nan([length(er) 1]);

FWHM_bin(i)=nan([length(er) 1]);
err_FWHM_bin(i)=nan([length(er) 1]);

cd (folder_name);
files=dir('*.h5');

i0=zeros([length(files) 1]);
area=zeros([length(files) 1]);
FWHM=zeros([length(files) 1]);

for i=1:length(files)
i0_par=mean(double(h5read(files(i).name,'/photon_diagnostics/FEL02/I0_monitor/iom_sh_a'))); %change /photon_diagnostics/FEL02/I0_monitor/iom_sh_a to /photon_diagnostics/Spectrometer/hor_area to change the 'i0' variable
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

dlmwrite(out_file,[i0,area,FWHM],'  ');

for i=1:length(er)-1
cond=i0>=er(i) & i0<er(i+1);
i0_bin(i)=mean(i0(cond));

area_bin(i)=mean(area(cond));
err_area_bin(i)=std(area(cond))./sqrt(sum(cond));

FWHM_bin(i)=mean(FWHM(cond));
err_FWHM_bin(i)=std(FWHM(cond))./sqrt(sum(cond));

end


cond2=~isnan(i0_bin);
i0_plot=i0_bin(cond2);
area_plot=area_bin(cond2);
err_area=err_area_bin(cond2);

FWHM_plot=FWHM_bin(cond2);
err_FWHM=err_FWHM_bin(cond2);


figure(2);
subplot(2,1,1);
errorbar(i0_plot,area_plot,err_area,'.');
xlabel('I0 (\muJ');
ylabel('Signal area');
subplot(2,1,2);
errorbar(i0_plot,FWHM_plot,err_FWHM,'.');
xlabel('I0 (\muJ');
ylabel('FWHM peak (pixels)');

dlmwrite(out_file_bin,[i0_plot,area_plot,err_area,FWHM_plot,err_FWHM],'  ');







