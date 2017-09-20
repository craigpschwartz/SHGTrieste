clear
close all
%Subfolder Directories [NS = No Sample; R = Raster (Sample)]
NS_subfolder = 'graph500nm/4_02nm_NS_006';
R_subfolder = 'graph500nm/4_02nm_R_006';

%Folder Directories
NS_Folder_Directory = ['/Users/rlam879/Desktop/FERMI_Analysis/' NS_subfolder];
R_Folder_Directory = ['/Users/rlam879/Desktop/FERMI_Analysis/' R_subfolder];

%eliminate bins with fewer than # points
bin_count_min=0;

%Line Drawing Parameters (OFF=0; ON=1)
Draw_lines=1;%turn on horizontal lines - set values with fun_abs and SHG_abs
ratio_draw=0;%turn on direct ratio calc

Fermi_analysis();