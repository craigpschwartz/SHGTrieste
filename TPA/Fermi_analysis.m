%Initialize Figure Counting
figure_count=1;

%Number of Signals (Fundamental, SHG)
num_sig = 1;

%Wavelength & Sample thickness parameterization
%find wavelength and sample thickness from subfolder directory
if strfind(R_subfolder,'4_02')>0
    wavelength=4.02;
elseif strfind(R_subfolder,'4_36')>0
    wavelength=4.36;
elseif strfind(R_subfolder,'4_76')>0
    wavelength=4.76;
end

if strfind(R_subfolder,'500')>0
    thickness=500;
elseif strfind(R_subfolder,'300')>0
    thickness=300;
elseif strfind(R_subfolder,'100')>0
    thickness=100;
elseif strfind(R_subfolder,'720')>0
    thickness=720;
end

%Define Scaling Factors (ADU/input photon) & Pixel Limits on CCD Image
if wavelength==4.02 %Wavelength=4.0278, Energy=307.86 eV
    ccd_scale_factor_fun=0.35;
    ccd_mult_factor=114;
    Al_transmission_fun=0.01912;
    Energy_eV=307.86;
    
    if thickness==300 %transmission through graphite
        C_transmission_fun=0.057562;
    elseif thickness==500
        C_transmission_fun=0.0085814;
    elseif thickness==100
        C_transmission_fun=0.38611;
    elseif thickness==720
        C_transmission_fun=0.0010576;
    end
    
    gas_cell_transmission_fun=0.035651; %N2 gas, 6 m, 1.7 mbar (1.2751 torr)
    %Scaling by gas_cell_attenuation factor is not implemented in current
    %code. Only factors that attenuate the beam after the sample are
    %included
    
    %Fundamental (a,b) and SH Signal (c,d). a & c correspond to the low
    %limits and b & d to the high limit of the ROIs.
    a=685;
    b=740;
    
elseif wavelength==4.36 %Wavelength=4.3635, Energy=284.18 eV
    ccd_scale_factor_fun=0.29;
    ccd_mult_factor=105;
    Al_transmission_fun=0.0078964;
    Energy_eV=284.18;
    
    if thickness==300 %transmission through graphite
        C_transmission_fun=0.63460;
    elseif thickness==500
        C_transmission_fun=0.46864;
    elseif thickness==100
        C_transmission_fun=0.85934;
    elseif thickness==720
        C_transmission_fun=0.33574;
    end
    
    gas_cell_transmission_fun=0.016430;
    
    a=620;
    b=670;
    
elseif wavelength==4.76 %Wavelength=4.7602, Energy=260.49 eV
    ccd_scale_factor_fun=0.23;
    ccd_mult_factor=90;
    Al_transmission_fun=0.0026371;
    Energy_eV=260.49;
    
    if thickness==300 %transmission through graphite
        C_transmission_fun=0.85287;
    elseif thickness==500
        C_transmission_fun=0.76701;
    elseif thickness==100
        C_transmission_fun=0.94833;
    elseif thickness==720
        C_transmission_fun=0.68252;
    end
    
    gas_cell_transmission_fun=0.0057603;
    
    a=540;
    b=600;
end

%Static Filter Parameters; FWHM filter is dynamically determined below
vert_sigma_filter=0.58;
hor_sigma_filter=0.006;
drain_current_threshold=0.5;

%Data Import
%import NS i0 and areas
cd (NS_Folder_Directory);
NS_files=dir('*.h5');

NS_drain_current=nan([length(NS_files) 1]);
NS_area=zeros([length(NS_files) num_sig]);
NS_area_bks_fil=nan([length(NS_files) num_sig]); %filtered data
NS_area_bks_vs_removed=nan([length(NS_files) num_sig]); %vert_sigma removed
NS_area_bks_hs_removed=nan([length(NS_files) num_sig]); %hor_sigma removed
NS_area_bks_fwhm_removed=nan([length(NS_files) num_sig]); %fwhm removed
NS_area_bks_dis_removed=nan([length(NS_files) num_sig]);

x_axis_x_proj_fun=a:b; %x axis for fwhm calculations

NS_x_proj_fwhm_fun=nan([length(NS_files) 1]); %fwhm data 

for i=1:length(NS_files)
    %filter files
    hor_sigma=mean(double(h5read(NS_files(i).name,...
        '/photon_diagnostics/Spectrometer/hor_sigma')));
    vert_sigma=mean(double(h5read(NS_files(i).name,...
        '/photon_diagnostics/Spectrometer/vert_sigma')));
    NS_drain_current(i)=double(h5read(NS_files(i).name,...
        '/eis-timex/caen_area_ch1')); %import i0 values
        %Import Chip Image
        CCD_Image=double(h5read(NS_files(i).name, '/CCD/Image'));
        
        %separate fundamental and sh signals from image
        fun_ROI=CCD_Image(a:b,:);
        
        %calculate x-projection from ROI
        x_proj_fun=sum(transpose(fun_ROI));
        
        %x-projection baseline correction
        x_proj_fun=x_proj_fun-mean(x_proj_fun(1:20));
        
        %x-projection FWHM
            %sum baseline corrected x-proj
            %index 1=fundamental; 2=signal
            NS_area_bks_fil(i,1)=sum(x_proj_fun);
            
            %scale ccd signal to account for the global response of the
            %spectrometer (ccd_scale_factor_*) and transmission through the
            %600 nm Aluminum filter
            
            NS_area_bks_fil(i,1)=NS_area_bks_fil(i,1)...
                /Al_transmission_fun/ccd_scale_factor_fun;
end

%import R i0 and areas
cd (R_Folder_Directory)
R_files=dir('*.h5');

R_drain_current=nan([length(R_files) 1]);
R_area=zeros([length(R_files) num_sig]);
R_area_bks_fil=nan([length(R_files) num_sig]);
R_area_bks_vs_removed=nan([length(R_files) num_sig]);
R_area_bks_hs_removed=nan([length(R_files) num_sig]);
R_area_bks_fwhm_removed=nan([length(R_files) num_sig]);
R_area_bks_dis_removed=nan([length(R_files) num_sig]);
R_x_proj_fwhm_fun=nan([length(R_files) 1]);

for i=1:length(R_files)
    %filter files
    hor_sigma=mean(double(h5read(R_files(i).name,...
        '/photon_diagnostics/Spectrometer/hor_sigma')));
    vert_sigma=mean(double(h5read(R_files(i).name,...
        '/photon_diagnostics/Spectrometer/vert_sigma')));
    R_drain_current(i)=double(h5read(R_files(i).name,...
        '/eis-timex/caen_area_ch1'));
        %Import Chip Image
        CCD_Image=double(h5read(R_files(i).name, '/CCD/Image'));
        
        %separate fundamental and sh signals from image
        fun_ROI=CCD_Image(a:b,:);
        
        %calculate x-projection from ROI
        x_proj_fun=sum(transpose(fun_ROI));
        
        %x-projection baseline correction
        x_proj_fun=x_proj_fun-mean(x_proj_fun(1:20));
        
        %x-projection FWHM
            %sum baseline corrected x-proj
            %index 1=fundamental; 2=signal
            R_area_bks_fil(i,1)=sum(x_proj_fun);
            
            %scale ccd signal to account for the global response of the
            %spectrometer (ccd_scale_factor_*) and transmission through the
            %600 nm Aluminum filter
            
            R_area_bks_fil(i,1)=R_area_bks_fil(i,1)...
                /Al_transmission_fun/ccd_scale_factor_fun;
        
    
end

NS_x_proj_fwhm_box_plot=boxplot(NS_x_proj_fwhm_fun,'Labels','Off-Sample');
NS_x_proj_fwhm_fun_filter=mean(get(NS_x_proj_fwhm_box_plot(3),'ydata'));
R_x_proj_fwhm_box_plot=boxplot(R_x_proj_fwhm_fun,'Labels','On-Sample');
R_x_proj_fwhm_fun_filter=mean(get(R_x_proj_fwhm_box_plot(3),'ydata'));


%Filter remaining points
%Renomalize i0 and signal axes to put everything on the same axis
NS_drain_current_renorm=NS_drain_current/max(NS_drain_current);
R_drain_current_renorm=R_drain_current/max(R_drain_current);

NS_area_bks_fil_renorm(:,1)=NS_area_bks_fil(:,1)/max(NS_area_bks_fil(:,1));

R_area_bks_fil_renorm(:,1)=R_area_bks_fil(:,1)/max(R_area_bks_fil(:,1));

%Calculate differences between sig and drain current
NS_area_bks_fil_renorm_diff(:,1)=NS_area_bks_fil_renorm(:,1)...
    -NS_drain_current_renorm;
R_area_bks_fil_renorm_diff(:,1)=R_area_bks_fil_renorm(:,1)...
    -R_drain_current_renorm;

%Plot Background Subtracted Data
figure(figure_count)
clf
figure_count=figure_count+1;
for j=1:num_sig
    subplot(2,num_sig,j)
    hold on
    plot(NS_drain_current,NS_area_bks_vs_removed(:,j),...
        '+','MarkerEdgeColor','r');
    plot(NS_drain_current,NS_area_bks_hs_removed(:,j),...
        'x','MarkerEdgeColor','m');
    plot(NS_drain_current,NS_area_bks_fwhm_removed(:,j),...
        '*','MarkerEdgeColor','g');
    plot(NS_drain_current,NS_area_bks_dis_removed(:,j),...
        '^','MarkerEdgeColor','c');
    plot(NS_drain_current,NS_area_bks_fil(:,j),'o','MarkerEdgeColor','b');
    xlabel('Drain Current (\muJ)');
    ylabel('Signal (photons)');
    if j==1
        title('NS Fundamental')
    else
        title('NS Signal')
    end
    
    subplot(2,num_sig,j+num_sig)
    hold on
    plot(R_drain_current,R_area_bks_vs_removed(:,j),...
        '+','MarkerEdgeColor','r');
    plot(R_drain_current,R_area_bks_hs_removed(:,j),...
        'x','MarkerEdgeColor','m');
    plot(R_drain_current,R_area_bks_fwhm_removed(:,j),...
        '*','MarkerEdgeColor','g');
    plot(R_drain_current,R_area_bks_dis_removed(:,j),...
        '^','MarkerEdgeColor','c');
    plot(R_drain_current,R_area_bks_fil(:,j),'o','MarkerEdgeColor','b');
    xlabel('Drain Current (\muJ)');
    ylabel('Signal (photons)');
    if j==1
        title('S Fundamental')
    else
        title('S Signal')
    end
end

%Binning Initialization

%Binning Parameters
%low end bins 'l'
%num_l_bins=7;

 bin_l_start=round((quantile(R_drain_current,0.25)...
      -1.5*iqr(R_drain_current))*2)/2;
%bin_l_start=round((quantile(R_drain_current,0.05)));
%bin_l_start=4.5;
%bin_l_end=round(quantile(R_drain_current,bin_l_value)*2)/2;
bin_l_end=0.0;
bin_l_step=.5;

%high end bins 'h'
bin_h_start=bin_l_end;
%bin_h_end=round(quantile(R_drain_current,bin_h_value)*2)/2;
bin_h_end=15;
bin_h_step=.5;

numbins=(bin_l_end-bin_l_start)/bin_l_step + ...
    (bin_h_end-bin_h_start)/bin_h_step;
numbins_l=(bin_l_end-bin_l_start)/bin_l_step;
%specify bins based on Input Binning Parameters
bins_l=(bin_l_start+bin_l_step/2):bin_l_step:bin_l_end;
num_l_bins=length(bins_l);
bins_h=(bin_h_start+bin_h_step/2):bin_h_step:bin_h_end;

%combine 'l' and 'h' bins into single array
bins=[bins_l bins_h];

NS_area_bin=nan([length(bins) 2]);
NS_err_area_bin=nan([length(bins) 2]);
NS_dc_bin=nan([length(bins) 1]);
NS_err_dc_bin=nan([length(bins) 1]);
NS_bin_counts=zeros([length(bins) 2]);
R_area_bin=nan([length(bins) 2]);
R_err_area_bin=nan([length(bins) 2]);
R_bin_counts=zeros([length(bins) 2]);
R_bin_number=zeros(length(bins));
R_dc_bin=nan([length(bins) 1]);
R_err_dc_bin=nan([length(bins) 1]);

for i=1:numbins
    count=0;
    NS_dc_accum=nan;
    NS_area_accum=nan;
    NS_area_accum2=nan;
    
    if i<=num_l_bins
        binmin=(bin_l_start+bin_l_step*(i-1));
        binmax=(bin_l_start + bin_l_step*i);
    else
        binmin=(bin_h_start+bin_h_step*(i-num_l_bins-1));
        binmax=(bin_h_start + bin_h_step*(i-num_l_bins));
    end
    
    for j=1:length(NS_drain_current)
        if NS_drain_current(j) > binmin && NS_drain_current(j) <= binmax && isnan(NS_area_bks_fil(j,1))==0
            count=count+1;
            NS_dc_accum(count,1)=NS_drain_current(j);
            NS_area_accum(count,1)=NS_area_bks_fil(j,1);
        end
    end
    
    if length(NS_area_accum(~isnan(NS_area_accum)))>bin_count_min
        NS_dc_bin(i,1)=nanmean(NS_dc_accum);
        NS_err_dc_bin(i,1)=nanstd(NS_dc_accum)/...
            sqrt(length(NS_dc_accum(~isnan(NS_dc_accum))));
        NS_err_dc_bin(i,1)=nanstd(NS_dc_accum);
        NS_area_bin(i,1)=nanmean(NS_area_accum);
        NS_area_bin(i,2)=nanmean(NS_area_accum2);
        NS_err_area_bin(i,1)=nanstd(NS_area_accum)/...
            sqrt(length(NS_area_accum(~isnan(NS_area_accum))));
        NS_err_area_bin(i,2)=nanstd(NS_area_accum2)/...
            sqrt(length(NS_area_accum2(~isnan(NS_area_accum2))));
    end
    NS_bin_counts(i,1)=length(NS_area_accum(~isnan(NS_area_accum)));
    NS_bin_counts(i,2)=length(NS_area_accum2(~isnan(NS_area_accum2)));
    
    count=0;
    R_dc_accum=nan;
    R_area_accum=nan;
    R_area_accum2=nan;
    
    for j=1:length(R_drain_current)
        if R_drain_current(j) > binmin && R_drain_current(j) <= binmax && isnan(R_area_bks_fil(j,1))==0
            count=count+1;
            R_dc_accum(count,1)=R_drain_current(j);
            R_area_accum(count,1)=R_area_bks_fil(j,1);
        end
    end
    if length(R_area_accum(~isnan(R_area_accum))) > bin_count_min
        %R_bin_number(1,i)=length(R_area_accum);%Save Var num of points.
        R_dc_bin(i,1)=nanmean(R_dc_accum);
        R_err_dc_bin(i,1)=nanstd(R_dc_accum)/...
            sqrt(length(R_dc_accum(~isnan(R_dc_accum))));
        R_err_dc_bin(i,1)=nanstd(R_dc_accum);
        R_area_bin(i,1)=nanmean(R_area_accum);
        R_area_bin(i,2)=nanmean(R_area_accum2);
        R_err_area_bin(i,1)=nanstd(R_area_accum)/...
            sqrt(length(R_area_accum(~isnan(R_area_accum))));
        R_err_area_bin(i,2)=nanstd(R_area_accum2)/...
            sqrt(length(R_area_accum2(~isnan(R_area_accum2))));
    end
    R_bin_counts(i,1)=length(R_area_accum(~isnan(R_area_accum)));
    R_bin_counts(i,2)=length(R_area_accum2(~isnan(R_area_accum2)));
    
end

%Remove nan rows for fitting
NS_fun_bin=zeros([length(bins) 2]);
NS_fun_bin(:,1)=NS_dc_bin;
NS_fun_bin(:,2)=NS_area_bin(:,1);
NS_fun_bin=NS_fun_bin(~any(isnan(NS_fun_bin),2),:);


R_fun_bin=zeros([length(bins) 2]);
R_fun_bin(:,1)=R_dc_bin;
R_fun_bin(:,2)=R_area_bin(:,1);
R_fun_bin=R_fun_bin(~any(isnan(R_fun_bin),2),:);




NS_fun_fit=fit(NS_fun_bin(:,1),NS_fun_bin(:,2),'power2');
R_fun_fit=fit(R_fun_bin(:,1),R_fun_bin(:,2),'power2');

figure(figure_count)
clf
figure_count=figure_count+1;
for j=1:num_sig
    subplot(2,num_sig,j)
    hold on
    errorbar(NS_dc_bin,NS_area_bin(:,j),NS_err_area_bin(:,j),'o');
    if j==1
        title('NS Fundamental Binned');
        plot(NS_fun_fit,'k')
    else
        title('NS Signal Binned')
        plot(NS_sig_fit,'k')
    end
    xlabel('Drain Current (\muJ)');
    ylabel('Signal (photons)');
    
    subplot(2,num_sig,j+num_sig)
    hold on
    errorbar(R_dc_bin,R_area_bin(:,j),R_err_area_bin(:,j),'o');
    if j==1
        title('S Fundamental Binned');
        plot(R_fun_fit,'k')
    else
        title('S Signal Binned')
        plot(R_sig_fit,'k')
    end
    xlabel('Drain Current (\muJ)');
    ylabel('Signal (photons)');
end

%bin division
ratio_fun=R_area_bin(:,1)./NS_area_bin(:,1);
err_ratio_fun=zeros(1,length(bins));
for k=1:length(bins)
    err_ratio_fun(k)=ratio_fun(k)*sqrt(((NS_err_area_bin(k,1)/...
        NS_area_bin(k,1)))^2+((R_err_area_bin(k,1)/R_area_bin(k,1)))^2);
end

ratio_fundamental(:,1)=bins;
ratio_fundamental(:,2)=ratio_fun;
ratio_fundamental=ratio_fundamental(~any(isnan(ratio_fundamental),2),:);
ratio_fun_fit=fit(ratio_fundamental(:,1),ratio_fundamental(:,2),'power2');

M=max(ratio_fundamental(:,1))+0.5*bin_h_step;
N=min(ratio_fundamental(:,1))-0.5*bin_l_step;
if ratio_draw==1
    ratio=zeros(10000,3);
    temp=zeros(10000,4);
    for i=1:10000
        j=((M-N)/10000)*i+N;
        temp(i,1)=feval(R_fun_fit,j);
        temp(i,2)=feval(NS_fun_fit,j);
        ratio(i,1)=j;
        ratio(i,2)=temp(i,1)/temp(i,2);
    end
end

figure(figure_count)
clf
figure_count=figure_count+1;
hold on
errorbar(bins,ratio_fun,err_ratio_fun, 'o');
%plot(ratio_fun_fit,'k');
if ratio_draw==1
    hold on
    plot(ratio(:,1),ratio(:,2),'r','LineWidth',1.5)
end
if Draw_lines==1
    hold on
    plot([N M],[C_transmission_fun C_transmission_fun])
end
title('Fundamental');
xlabel('Drain Current (\muJ)');
ylabel('Signal Ratio (I^{0}_{S}/I^{0}_{NS})');
%SHG Signal - Scale and Subtract

%Calculate Powers and power densities
Pulse_Energy(:,1)=feval(NS_fun_fit,R_dc_bin)*Energy_eV*1.60218e-13*ccd_mult_factor;
Photons_per_pulse(:,1)=feval(NS_fun_fit,R_dc_bin)*ccd_mult_factor;
Power_per_areasec(:,1)=Pulse_Energy(:,1)/(pulse_length*spot_size);
Photons_per_areasec(:,1)=Photons_per_pulse(:,1)/(pulse_length*spot_size);


%setup line for eyes
Photons_per_areasec(~any(~isnan(Photons_per_areasec), 2),:)=[];%remove nan
ratio_fun(~any(~isnan(ratio_fun), 2),:)=[];%remove nan
err_ratio_fun=err_ratio_fun.';
err_ratio_fun(~any(~isnan(err_ratio_fun), 2),:)=[];
L=length(Photons_per_areasec);

%count=0;%UGLY ASS KLUGE TO DEAL WITH SPARSE BINS
%for i=1:L
%    if Photons_per_areasec(i) >= 0
%        count = count +1;
%        Photons_per_area_sec_trim(count,1)=Photons_per_areasec(i);
%    end
%end    
%L=length(Photons_per_area_sec_trim);

L_short=round(L/2);
ratio_pperas_fit=fit(Photons_per_areasec(cutoff:L),ratio_fun(cutoff:L),'poly1');
%setup straight line to show decay
line_pperas=mean(ratio_fun(cutoff:L_short));


%setup plotting
%lower limits
%M=Photons_per_areasec(1)*0.9;For Variable upper limit
M=0;
%upperlimit
N=Photons_per_areasec(L)*1.1;
curves=zeros(10000,3);
for i=1:10000
    j=((M-N)/10000)*i+N;
    curves(i,1)=j;%index
    curves(i,2)=feval(ratio_pperas_fit,j);%not flat line
    curves(i,3)=line_pperas;%flat line
end

%lowerplot limit
M_plot=Photons_per_areasec(cutoff)*0.9;

figure(figure_count)
clf
figure_count=figure_count+1;
hold on
errorbar(Photons_per_areasec,ratio_fun(1:L),err_ratio_fun(1:L), 'o');
%plot(ratio_fun_fit,'k');
hold on
plot(curves(:,1),curves(:,2),'r','LineWidth',1.5)
hold on
plot(curves(:,1),curves(:,3),'k','LineWidth',1.5)
xlim([M_plot,N])
title('Fundamental');
xlabel('Pulse energy/area second (J/cm^2s)');
ylabel('Signal Ratio (I^{0}_{S}/I^{0}_{NS})');

output=horzcat(Photons_per_areasec,ratio_fun(1:L),err_ratio_fun);
curves=flipud(curves);





