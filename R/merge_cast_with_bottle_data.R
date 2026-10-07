# Merge cast and bottle data with data.R

# import the bottle and cast csv files
bottle_data <- read.csv("/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/194903-202105_Bottle.csv")

cast_data <- read.csv("/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/194903-202105_Cast.csv") 

names(bottle_data)
names(cast_data)

#subset cast for variables need for merge with bottle data
cast_data_subset <- subset(cast_data, select = c(Cruise_ID, Year, Date, Cst_Cnt, Sta_ID, Rpt_Line, Rpt_Sta, Lat_Dec, Lon_Dec))
cast_data_subset[1:5,]
dim(cast_data_subset)
#subset bottle data for variables to merge with cast data
bottle_data_subset <- subset(bottle_data, select = c(Cst_Cnt, Sta_ID, R_Depth, R_TEMP, R_Sal, R_Oxy_mol.Kg))
bottle_data_subset[1:5,]
dim(bottle_data_subset)

#merge cast and bottle data by Cst_Cnt, Sta_ID  
merged_data <- merge(cast_data_subset, bottle_data_subset, by = c("Cst_Cnt", "Sta_ID"))
merged_data[1:10,]
dim(merged_data)

# extract one year
merged_2020 <-merged_data[merged_data$Year == "2020",c("Cruise_ID", "Lat_Dec", "Lon_Dec", "R_Depth", "Rpt_Line", "Rpt_Sta", "R_TEMP")]
merged_2020[1:10,]

# extract one cruise in that year
cruise <-unique(merged_2020$Cruise_ID)
merged_2020_cruise <- merged_2020[merged_2020$Cruise_ID == cruise[1], c("Lon_Dec", "Lat_Dec", "R_Depth", "Rpt_Line", "Rpt_Sta", "R_TEMP")]
merged_2020_cruise[1:10,]

# extract 50m depth level
merged_2020_cruise_depth_50m <- merged_2020_cruise[merged_2020_cruise$R_Depth == 50, c("Lon_Dec", "Lat_Dec", "R_Depth", "Rpt_Line", "Rpt_Sta", "R_TEMP")]
merged_2020_cruise_depth_50m[1:10,]

# save 50m dataset 
write.table(merged_2020_cruise_depth_50, file = "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/merged_2020_cruise_depth_50m.csv", sep = ",", row.names = FALSE)

# extract 100m depth level
merged_2020_cruise_depth_100m <- merged_2020_cruise[merged_2020_cruise$R_Depth == 100, c("Lon_Dec", "Lat_Dec", "R_Depth", "Rpt_Line", "Rpt_Sta", "R_TEMP")]
merged_2020_cruise_depth_100m[1:10,]

# save 100m dataset 
write.table(merged_2020_cruise_depth_100m, file = "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/merged_2020_cruise_depth_100m.csv", sep = ",", row.names = FALSE)