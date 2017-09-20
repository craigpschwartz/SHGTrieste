e%Initialize Figure Counting
figure_count=1;

%Number of Signals (Fundamental, SHG)
num_sig = 2;

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
    ccd_scale_factor_sh=1.30;
    Al_transmission_fun=0.01912;
    Al_transmission_sh=0.48817;
    Energy_eV=307.86;
    
    if thickness==300 %transmission through graphite
        C_transmission_fun=0.057562;
        C_transmission_sh=0.58625;
    elseif thickness==500
        C_transmission_fun=0.0085814;
        C_transmission_sh=0.41065;
    elseif thickness==100
        C_transmission_fun=0.38611;
        C_transmission_sh=0.83694;
    elseif thickness==720
        C_transmission_fun=0.0010576;
        C_transmission_sh=0.27758;
    end
    
    gas_cell_transmission_fun=0.035651; %N2 gas, 6 m, 1.7 mbar (1.2751 torr)
    %Scaling by gas_cell_attenuation factor is not implemented in current
    %code. Only factors that attenuate the beam after the sample are
    %included
    
    %Fundamental (a,b) and SH Signal (c,d). a & c correspond to the low
    %limits and b & d to the high limit of the ROIs.
    a=685;
    b=740;
    c=1165;
    d=1225;
    
elseif wavelength==4.36 %Wavelength=4.3635, Energy=284.18 eV
    ccd_scale_factor_fun=0.29;
    ccd_scale_factor_sh=1.15;
    Al_transmission_fun=0.0078964;
    Al_transmission_sh=0.41195;
    Energy_eV=284.18;
    
    if thickness==300 %transmission through graphite
        C_transmission_fun=0.63460;
        C_transmission_sh=0.52008;
    elseif thickness==500
        C_transmission_fun=0.46864;
        C_transmission_sh=0.33635;
    elseif thickness==100
        C_transmission_fun=0.85934;
        C_transmission_sh=0.80419;
    elseif thickness==720
        C_transmission_fun=0.33574;
        C_transmission_sh=0.20825;
    end
    
    gas_cell_transmission_fun=0.016430;
    
    a=620;
    b=670;
    c=1120;
    d=1170;
    
elseif wavelength==4.76 %Wavelength=4.7602, Energy=260.49 eV
    ccd_scale_factor_fun=0.23;
    ccd_scale_factor_sh=0.99;
    Al_transmission_fun=0.0026371;
    Al_transmission_sh=0.33014;
    Energy_eV=260.49;
    
    if thickness==300 %transmission through graphite
        C_transmission_fun=0.85287;
        C_transmission_sh=0.44278;
    elseif thickness==500
        C_transmission_fun=0.76701;
        C_transmission_sh=0.25723;
    elseif thickness==100
        C_transmission_fun=0.94833;
        C_transmission_sh=0.76219;
    elseif thickness==720
        C_transmission_fun=0.68252;
        C_transmission_sh=0.14153;
    end
    
    gas_cell_transmission_fun=0.0057603;
    
    a=540;
    b=600;
    c=1060;
    d=1120;
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
x_axis_x_proj_sh=c:d;

NS_x_proj_fwhm_fun=nan([length(NS_files) 1]); %fwhm data 
NS_x_proj_fwhm_sh=nan([length(NS_files) 1]);

for i=1:length(NS_files)
    %filter files
    hor_sigma=mean(double(h5read(NS_files(i).name,...
        '/photon_diagnostics/Spectrometer/hor_sigma')));
    vert_sigma=mean(double(h5read(NS_files(i).name,...
        '/photon_diagnostics/Spectrometer/vert_sigma')));
    NS_drain_current(i)=double(h5read(NS_files(i).name,...
        '/eis-timex/caen_area_ch1')); %import i0 values
    if NS_drain_current(i)>drain_current_threshold
        %Import Chip Image
        CCD_Image=double(h5read(NS_files(i).name, '/CCD/Image'));
        
        %separate fundamental and sh signals from image
        fun_ROI=CCD_Image(a:b,:);
        sh_ROI=CCD_Image(c:d,:);
        
        %calculate x-projection from ROI
        x_proj_fun=sum(transpose(fun_ROI));
        x_proj_sh=sum(transpose(sh_ROI));
        
        %x-projection baseline correction
        x_proj_fun=x_proj_fun-mean(x_proj_fun(1:20));
        x_proj_sh=x_proj_sh-mean(x_proj_sh(1:20));
        
        %x-projection FWHM
        NS_x_proj_fwhm_fun(i)=fwhm(x_axis_x_proj_fun,x_proj_fun);
        NS_x_proj_fwhm_sh(i)=fwhm(x_axis_x_proj_sh,x_proj_sh);
        if vert_sigma>=vert_sigma_filter && hor_sigma<=hor_sigma_filter
            %sum baseline corrected x-proj
            %index 1=fundamental; 2=signal
            NS_area_bks_fil(i,1)=sum(x_proj_fun);
            NS_area_bks_fil(i,2)=sum(x_proj_sh);
            
            %scale ccd signal to account for the global response of the
            %spectrometer (ccd_scale_factor_*) and transmission through the
            %600 nm Aluminum filter
            
            NS_area_bks_fil(i,1)=NS_area_bks_fil(i,1)...
                /Al_transmission_fun/ccd_scale_factor_fun;
            NS_area_bks_fil(i,2)=NS_area_bks_fil(i,2)...
                /Al_transmission_sh/ccd_scale_factor_sh;
        elseif vert_sigma<vert_sigma_filter && hor_sigma<=hor_sigma_filter
            %sum baseline corrected x-proj
            %index 1=fundamental; 2=signal
            NS_area_bks_vs_removed(i,1)=sum(x_proj_fun);
            NS_area_bks_vs_removed(i,2)=sum(x_proj_sh);
            
            %scale ccd signal to account for the global response of the
            %spectrometer (ccd_scale_factor_*) and transmission through the
            %600 nm Aluminum filter
            
            NS_area_bks_vs_removed(i,1)=NS_area_bks_vs_removed(i,1)...
                /Al_transmission_fun/ccd_scale_factor_fun;
            NS_area_bks_vs_removed(i,2)=NS_area_bks_vs_removed(i,2)...
                /Al_transmission_sh/ccd_scale_factor_sh;
        else
                        %sum baseline corrected x-proj
            %index 1=fundamental; 2=signal
            NS_area_bks_hs_removed(i,1)=sum(x_proj_fun);
            NS_area_bks_hs_removed(i,2)=sum(x_proj_sh);
            
            %scale ccd signal to account for the global response of the
            %spectrometer (ccd_scale_factor_*) and transmission through the
            %600 nm Aluminum filter
            
            NS_area_bks_hs_removed(i,1)=NS_area_bks_hs_removed(i,1)...
                /Al_transmission_fun/ccd_scale_factor_fun;
            NS_area_bks_hs_removed(i,2)=NS_area_bks_hs_removed(i,2)...
                /Al_transmission_sh/ccd_scale_factor_sh;
        end
    end
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
    if R_drain_current(i)>drain_current_threshold
        %Import Chip Image
        CCD_Image=double(h5read(R_files(i).name, '/CCD/Image'));
        
        %separate fundamental and sh signals from image
        fun_ROI=CCD_Image(a:b,:);
        sh_ROI=CCD_Image(c:d,:);
        
        %calculate x-projection from ROI
        x_proj_fun=sum(transpose(fun_ROI));
        x_proj_sh=sum(transpose(sh_ROI));
        
        %x-projection baseline correction
        x_proj_fun=x_proj_fun-mean(x_proj_fun(1:20));
        x_proj_sh=x_proj_sh-mean(x_proj_sh(1:20));
        
        %x-projection FWHM
        R_x_proj_fwhm_fun(i)=fwhm(x_axis_x_proj_fun,x_proj_fun);
        if vert_sigma>=vert_sigma_filter && hor_sigma<=hor_sigma_filter
            %sum baseline corrected x-proj
            %index 1=fundamental; 2=signal
            R_area_bks_fil(i,1)=sum(x_proj_fun);
            R_area_bks_fil(i,2)=sum(x_proj_sh);
            
            %scale ccd signal to account for the global response of the
            %spectrometer (ccd_scale_factor_*) and transmission through the
            %600 nm Aluminum filter
            
            R_area_bks_fil(i,1)=R_area_bks_fil(i,1)...
                /Al_transmission_fun/ccd_scale_factor_fun;
            R_area_bks_fil(i,2)=R_area_bks_fil(i,2)...
                /Al_transmission_sh/ccd_scale_factor_sh;
        elseif vert_sigma<vert_sigma_filter && hor_sigma<=hor_sigma_filter
            %sum baseline corrected x-proj
            %index 1=fundamental; 2=signal
            R_area_bks_vs_removed(i,1)=sum(x_proj_fun);
            R_area_bks_vs_removed(i,2)=sum(x_proj_sh);
            
            %scale ccd signal to account for the global response of the
            %spectrometer (ccd_scale_factor_*) and transmission through the
            %600 nm Aluminum filter
            
            R_area_bks_vs_removed(i,1)=R_area_bks_vs_removed(i,1)...
                /Al_transmission_fun/ccd_scale_factor_fun;
            R_area_bks_vs_removed(i,2)=R_area_bks_vs_removed(i,2)...
                /Al_transmission_sh/ccd_scale_factor_sh;
        else
                        %sum baseline corrected x-proj
            %index 1=fundamental; 2=signal
            R_area_bks_hs_removed(i,1)=sum(x_proj_fun);
            R_area_bks_hs_removed(i,2)=sum(x_proj_sh);
            
            %scale ccd signal to account for the global response of the
            %spectrometer (ccd_scale_factor_*) and transmission through the
            %600 nm Aluminum filter
            
            R_area_bks_hs_removed(i,1)=R_area_bks_hs_removed(i,1)...
                /Al_transmission_fun/ccd_scale_factor_fun;
            R_area_bks_hs_removed(i,2)=R_area_bks_hs_removed(i,2)...
                /Al_transmission_sh/ccd_scale_factor_sh;
        end
    end
end

%Plot Histograms of FWHM distribution for the fundamental
figure(figure_count)
clf
figure_count=figure_count+1;

subplot(2,2,1)
NS_x_proj_fwhm_hist=histogram(NS_x_proj_fwhm_fun,1000);
xlabel('FWHM (Pixels)');
ylabel('Count');
title('Off-Sample Fundamental Histogram')

subplot(2,2,2)
NS_x_proj_fwhm_box_plot=boxplot(NS_x_proj_fwhm_fun,'Labels','Off-Sample');
NS_x_proj_fwhm_fun_filter=mean(get(NS_x_proj_fwhm_box_plot(3),'ydata'));

ylabel('FWHM (Pixels)');

title('Off-Sample Fundamental Box Plot')
subplot(2,2,3)
R_x_proj_fwhm_hist=histogram(R_x_proj_fwhm_fun,1000);
xlabel('FWHM (Pixels)');
ylabel('Count');
title('On-Sample Fundamental Histogram')

subplot(2,2,4)
R_x_proj_fwhm_box_plot=boxplot(R_x_proj_fwhm_fun,'Labels','On-Sample');
R_x_proj_fwhm_fun_filter=mean(get(R_x_proj_fwhm_box_plot(3),'ydata'));

ylabel('FWHM (Pixels)');
title('On-Sample Fundamental Box Plot')

%filter by fwhm
for i=1:length(NS_files)
    if NS_x_proj_fwhm_fun(i) >= NS_x_proj_fwhm_fun_filter
        NS_area_bks_fwhm_removed(i,1)=NS_area_bks_fil(i,1);
        NS_area_bks_fwhm_removed(i,2)=NS_area_bks_fil(i,2);
        NS_area_bks_fil(i,1)=nan;
        NS_area_bks_fil(i,2)=nan;
    end
end

for i=1:length(R_files)
    if R_x_proj_fwhm_fun(i) >= R_x_proj_fwhm_fun_filter
        R_area_bks_fwhm_removed(i,1)=R_area_bks_fil(i,1);
        R_area_bks_fwhm_removed(i,2)=R_area_bks_fil(i,2);
        R_area_bks_fil(i,1)=nan;
        R_area_bks_fil(i,2)=nan;
    end
end

%Filter remaining points
%Renomalize i0 and signal axes to put everything on the same axis
NS_drain_current_renorm=NS_drain_current/max(NS_drain_current);
R_drain_current_renorm=R_drain_current/max(R_drain_current);

NS_area_bks_fil_renorm(:,1)=NS_area_bks_fil(:,1)/max(NS_area_bks_fil(:,1));
NS_area_bks_fil_renorm(:,2)=NS_area_bks_fil(:,2)/max(NS_area_bks_fil(:,2));

R_area_bks_fil_renorm(:,1)=R_area_bks_fil(:,1)/max(R_area_bks_fil(:,1));
R_area_bks_fil_renorm(:,2)=R_area_bks_fil(:,2)/max(R_area_bks_fil(:,2));

%Calculate differences between sig and drain current
NS_area_bks_fil_renorm_diff(:,1)=NS_area_bks_fil_renorm(:,1)...
    -NS_drain_current_renorm;
NS_area_bks_fil_renorm_diff(:,2)=NS_area_bks_fil_renorm(:,2)...
    -NS_drain_current_renorm;

R_area_bks_fil_renorm_diff(:,1)=R_area_bks_fil_renorm(:,1)...
    -R_drain_current_renorm;
R_area_bks_fil_renorm_diff(:,2)=R_area_bks_fil_renorm(:,2)...
    -R_drain_current_renorm;

%Calculate distances between (drain_current_renorm,sig) and (0.2,0.75)
NS_area_bks_fil_renorm_dis(:,1)=sqrt((NS_drain_current_renorm-.2).^2+...
    (NS_area_bks_fil_renorm(:,1)-.75).^2);
NS_area_bks_fil_renorm_dis(:,2)=sqrt((NS_drain_current_renorm-.2).^2+...
    (NS_area_bks_fil_renorm(:,2)-.75).^2);

R_area_bks_fil_renorm_dis(:,1)=sqrt((R_drain_current_renorm-.2).^2+...
    (R_area_bks_fil_renorm(:,1)-.75).^2);
R_area_bks_fil_renorm_dis(:,2)=sqrt((R_drain_current_renorm-.2).^2+...
    (R_area_bks_fil_renorm(:,2)-.75).^2);

NS_x_proj_fwhm_dis_filter_l=quantile(NS_area_bks_fil_renorm_dis(:,2),0.25)...
    -1.5*iqr(NS_area_bks_fil_renorm_dis(:,2));
NS_x_proj_fwhm_dis_filter_h=quantile(NS_area_bks_fil_renorm_dis(:,2),0.75)...
    +1.5*iqr(NS_area_bks_fil_renorm_dis(:,2));

R_x_proj_fwhm_dis_filter_l=quantile(R_area_bks_fil_renorm_dis(:,2),0.25)...
    -1.5*iqr(R_area_bks_fil_renorm_dis(:,2));
R_x_proj_fwhm_dis_filter_h=quantile(R_area_bks_fil_renorm_dis(:,2),0.75)...
    +1.5*iqr(R_area_bks_fil_renorm_dis(:,2));

for i=1:length(NS_files)
    if NS_area_bks_fil_renorm_dis(i,2) < NS_x_proj_fwhm_dis_filter_l ||...
            NS_area_bks_fil_renorm_dis(i,2) > NS_x_proj_fwhm_dis_filter_h
        NS_area_bks_dis_removed(i,1)=NS_area_bks_fil(i,1);
        NS_area_bks_dis_removed(i,2)=NS_area_bks_fil(i,2);
        NS_area_bks_fil(i,1)=nan;
        NS_area_bks_fil(i,2)=nan;
    end
end

for i=1:length(R_files)
    if R_area_bks_fil_renorm_dis(i,2) < R_x_proj_fwhm_dis_filter_l ||...
          R_area_bks_fil_renorm_dis(i,2) > R_x_proj_fwhm_dis_filter_h  
        R_area_bks_dis_removed(i,1)=R_area_bks_fil(i,1);
        R_area_bks_dis_removed(i,2)=R_area_bks_fil(i,2);
        R_area_bks_fil(i,1)=nan;
        R_area_bks_fil(i,2)=nan;
    end
end

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
bin_l_end=round(quantile(R_drain_current,0.25)*2)/2;
%bin_l_end=8.5;
bin_l_step=.5;

%high end bins 'h'
bin_h_start=bin_l_end;
bin_h_end=round(quantile(R_drain_current,.95)*2)/2;
bin_h_end=12;
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
            NS_area_accum2(count,1)=NS_area_bks_fil(j,2);
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
            R_area_accum2(count,1)=R_area_bks_fil(j,2);
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

NS_sig_bin=zeros([length(bins) 2]);
NS_sig_bin(:,1)=NS_dc_bin;
NS_sig_bin(:,2)=NS_area_bin(:,2);
NS_sig_bin=NS_sig_bin(~any(isnan(NS_sig_bin),2),:);

R_fun_bin=zeros([length(bins) 2]);
R_fun_bin(:,1)=R_dc_bin;
R_fun_bin(:,2)=R_area_bin(:,1);
R_fun_bin=R_fun_bin(~any(isnan(R_fun_bin),2),:);

R_sig_bin=zeros([length(bins) 2]);
R_sig_bin(:,1)=R_dc_bin;
R_sig_bin(:,2)=R_area_bin(:,2);
R_sig_bin=R_sig_bin(~any(isnan(R_sig_bin),2),:);



NS_fun_fit=fit(NS_fun_bin(:,1),NS_fun_bin(:,2),'power2');
NS_sig_fit=fit(NS_sig_bin(:,1),NS_sig_bin(:,2),'power2');

R_fun_fit=fit(R_fun_bin(:,1),R_fun_bin(:,2),'power2');
R_sig_fit=fit(R_sig_bin(:,1),R_sig_bin(:,2),'power2');

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
ratio_sig=R_area_bin(:,2)./NS_area_bin(:,2);
err_ratio_fun=zeros(1,length(bins));
err_ratio_sig=zeros(1,length(bins));
for k=1:length(bins)
    err_ratio_fun(k)=ratio_fun(k)*sqrt(((NS_err_area_bin(k,1)/...
        NS_area_bin(k,1)))^2+((R_err_area_bin(k,1)/R_area_bin(k,1)))^2);
    err_ratio_sig(k)=ratio_sig(k)*sqrt(((NS_err_area_bin(k,2)/...
        NS_area_bin(k,2)))^2+((R_err_area_bin(k,2)/R_area_bin(k,2)))^2);
end

ratio_fundamental(:,1)=bins;
ratio_fundamental(:,2)=ratio_fun;
ratio_fundamental=ratio_fundamental(~any(isnan(ratio_fundamental),2),:);
ratio_fun_fit=fit(ratio_fundamental(:,1),ratio_fundamental(:,2),'power2');

ratio_signal(:,1)=bins;
ratio_signal(:,2)=ratio_sig;
ratio_signal=ratio_signal(~any(isnan(ratio_signal),2),:);
ratio_sig_fit=fit(ratio_signal(:,1),ratio_signal(:,2),'power2');

M=max(ratio_fundamental(:,1))+0.5*bin_h_step;
N=min(ratio_fundamental(:,1))-0.5*bin_l_step;
if ratio_draw==1
    ratio=zeros(10000,3);
    temp=zeros(10000,4);
    for i=1:10000
        j=((M-N)/10000)*i+N;
        temp(i,1)=feval(R_fun_fit,j);
        temp(i,2)=feval(NS_fun_fit,j);
        temp(i,3)=feval(R_sig_fit,j);
        temp(i,4)=feval(NS_sig_fit,j);
        ratio(i,1)=j;
        ratio(i,2)=temp(i,1)/temp(i,2);
        ratio(i,3)=temp(i,3)/temp(i,4);
    end
end

figure(figure_count)
clf
figure_count=figure_count+1;
subplot(2,1,1)
hold on
errorbar(bins,ratio_fun,err_ratio_fun, 'o');
plot(ratio_fun_fit,'k');
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
subplot(2,1,2)
hold on
errorbar(bins,ratio_sig,err_ratio_sig, 'o');
plot(ratio_sig_fit,'k');
if ratio_draw==1
    hold on
    plot(ratio(:,1),ratio(:,3),'r','LineWidth',1.5)
end
if Draw_lines==1
    hold on
    plot([N M],[C_transmission_sh C_transmission_sh])
end
title('Signal');
xlabel('Drain Current (\muJ)');
ylabel('Signal ratio (I^{1}_{S}/I^{1}_{NS})');

%SHG Signal - Scale and Subtract

ratio_sig_fit_coeff=coeffvalues(ratio_sig_fit);
%SHG_scale_factor=ratio_sig_fit_coeff(3);
SHG_scale_factor=nanmean(ratio_sig(1:numbins_l))
SHG_scale_factor=.7*C_transmission_sh;
%SHG_scale_factor=0.119;
R_area_bin_scaled=R_area_bin(:,2)/SHG_scale_factor;
SHG_sig=R_area_bin_scaled-feval(NS_sig_fit,R_dc_bin);
Pulse_Energy(:,1)=feval(NS_fun_fit,R_dc_bin)*Energy_eV*1.60218e-13;

%Remove nan values for SHG_sig_fit
SHG_PE_sig=zeros([length(bins) 2]);
SHG_PE_sig(:,1)=Pulse_Energy;
SHG_PE_sig(:,2)=SHG_sig;
SHG_PE_sig=SHG_PE_sig(~any(isnan(SHG_PE_sig),2),:);
SHG_sig_fit=fit(SHG_PE_sig(:,1),SHG_PE_sig(:,2),'power2');
SHG_sig_y_err=R_err_area_bin(:,2)/SHG_scale_factor;

%SHG_sig_x_err(:,1): Negative Err
SHG_sig_x_err(:,1)=Pulse_Energy-feval(NS_fun_fit,R_dc_bin-R_err_dc_bin)...
    *Energy_eV*1.60218e-13;
%SHG_sig_x_err(:,2): Positive Err
SHG_sig_x_err(:,2)=feval(NS_fun_fit,R_dc_bin+R_err_dc_bin)...
    *Energy_eV*1.60218e-13-Pulse_Energy;

%Plot SHG_signal vs Input Energy (uJ); calculated from photon count &
%energy
figure(figure_count)
clf
figure_count=figure_count+1;
title('Power Spectrum')
hold on
errorbar(Pulse_Energy,SHG_sig,SHG_sig_y_err,SHG_sig_y_err,SHG_sig_x_err(:,1),SHG_sig_x_err(:,2),'o')
plot(SHG_sig_fit,'k')
xlabel('Pulse Energy (\muJ)')
ylabel('SH Power (Photons/Pulse)')

mirror_drain_current=2:14;
photons=feval(NS_fun_fit,mirror_drain_current);

