First open the "Etelis.Rproj" before running any R script. 
All scripts are wrapped into "wrapper.R"



# Description of directories in this folder:
	- 
	- 
	
	
# To download raw fastq files
wget -r -nH --cut-dirs=3 --header 'Authorization: Bearer opeztenR6APWV1GDwI2CnB5zlJgts5mqTnnkkY2G' https://ordering.diversityarrays.com/clients/2291/6630/Raw_Data_Index.html

	
# Variable names has been changed:
data_sites to data_stations
coord_site to data_sites
data_fishbase to data_taxo