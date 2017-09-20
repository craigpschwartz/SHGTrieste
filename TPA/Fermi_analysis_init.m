clear
close all
%Subfolder Directories [NS = No Sample; R = Raster (Sample)]
NS_subfolder = 'graph300nm/4_36nm_NS_combined5';
R_subfolder = 'graph300nm/4_36nm_R_combined5';

%Folder Directories
NS_Folder_Directory = ['/Users/craigschwartz/Documents/FERMI_SHG/' NS_subfolder];
R_Folder_Directory = ['/Users/craigschwartz/Documents/FERMI_SHG/' R_subfolder];


%eliminate bins with fewer than # points
bin_count_min=0;
%set bin minimium
bin_l_value=0.02;
bin_h_value=0.98;

%fitting cutoff
cutoff=1;

%Line Drawing Parameters (OFF=0; ON=1)
Draw_lines=0;%turn on horizontal lines - set values with fun_abs and SHG_abs
ratio_draw=0;%turn on direct ratio calc

%Laser Properties
pulse_length=35*10^-15;%Seconds
spot_size=350*10^-8;%cm^2

Fermi_analysis();