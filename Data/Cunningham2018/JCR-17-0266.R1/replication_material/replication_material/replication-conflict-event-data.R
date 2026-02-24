###########################################################################################
# Replication file for "Integrating Conflict Event Data"
# Authors: Donnay, Karsten, Eric Dunford, Erin McGrath, David Backer, and David Cunningham
# Journal of Conflict Resolution, accepted: March 2018
###########################################################################################
# Script: (1) Formats & integrates data from ACLED, UCDP-GED, SCAD, and GTD for
#            Nigeria 2011 using the meltt algorithm; 
#         (2) Recreates Figure 3.
#         (3) Recreates Table 3.
#         (4) Formats & integrates data from ACLED, UCDP-GED, SCAD, and GTD for
#            Libya 2014 and South Sudan 2015 using the meltt algorithm.
#
# Note: Tables 2 and 4 were produced through a manual review of the data and cannot 
#       be replicated automatically. See meltt_validate() for a method of calculating 
#       performance.
###########################################################################################


require(meltt) # v.4.0 (to integrate data)
require(tidyverse) # v.1.2.1 (to manage data)
require(ggmap) # v.2.6.1 (to generate map)
require(sf) #v.0.5 (for spatial data)

# Data --------------------------------------------------------------------

# Load .Rdata object of Nigeria 2011 Data from all 4 datasets. Data reflect
# operative 2015 versions. Also, loads the event, actor, and precision
# taxonomies for the Nigeria 2011 data and spatial data to generate spatial
# breakdown for Table 3.
load("Data/Nigeria_2011.Rdata")

# store all three taxonomies as lists
taxonomy = list(actor_tax=actor_tax,
                event_tax=event_tax,
                prec_tax=prec_tax)

# Formatting Data for Meltt -----------------------------------------------

acled <- 
  acled %>% 
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(EVENT_DATE),
         enddate = as.Date(EVENT_DATE),
         event_tax = EVENT_TYPE,
         actor_tax = ACTOR1,
         prec_tax = GEO_PRECISION,
         
         # Break up riots and protests
         event_tax = ifelse(event_tax =="Riots/Protests" &  FATALITIES > 0,"Riots",
                     ifelse(event_tax =="Riots/Protests" &  FATALITIES==0,"Protests",event_tax))) %>% 
  
  # convert all column names to lower case
  rename_all(tolower) %>% 
  as.data.frame

ged <- 
  ged %>%  
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(date_start),
         enddate = as.Date(date_end),
         event_tax = type_of_violence,
         actor_tax = side_a,
         prec_tax = where_prec) %>% 
  as.data.frame

scad <- 
  scad %>% 
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(startdate,"%d-%b-%y"),
         enddate = as.Date(enddate,"%d-%b-%y"),
         event_tax = etype,
         actor_tax = actor1,
         prec_tax = locnum) %>% 
  as.data.frame

gtd <- 
  gtd %>%  
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(paste0(iyear,"-",imonth,"-",iday)),
         enddate = date,
         event_tax = attacktype1,
         actor_tax = gname,
         prec_tax = specificity) %>% 
  
  # Drop entries missing geo-tagged information.
  filter(!is.na(longitude)) %>% 
  as.data.frame



# Integrate Data ----------------------------------------------------------


output = meltt(acled,ged,gtd,scad,taxonomies = taxonomy,
               twindow = 1,spatwindow = 3)

# Summary of the integration reported in the paper
summary(output)

# Figure 3 ----------------------------------------------------------------

# Generating a plot similar to panel 1 (barplot)
plot(output)

# Generating a plot similar to panel 2 (histogram)
tplot(output,time_unit="weeks",free_scale = FALSE)

# Generating a plot similar to panel 3 (map)
unique_entries = meltt_data(output,columns=c("longitude","latitude"))

duplicate_entries = meltt_duplicates(output,columns = "") 

convert = function(x,feature) x %>% 
  select(contains(feature)) %>% {colnames(.) <- c("dataset","event");.} %>% 
  mutate(dataset = ifelse(dataset!=0,feature,NA)) %>% drop_na

matches <-
  bind_rows(
    duplicate_entries %>% convert("acled"),
    duplicate_entries %>% convert("ged"),
    duplicate_entries %>% convert("gtd"),
    duplicate_entries %>% convert("scad")
  ) %>% mutate(match=1)

dat = left_join(unique_entries,matches) %>% 
  mutate(type=ifelse(is.na(match),"unique","duplicate"))

# Generate map
m <- get_map("nigeria",zoom=6,maptype = "toner") 
ggmap(m) + 
  geom_point(data = dat, 
             aes(longitude,latitude,color=dataset,
                 alpha=factor(type,levels=c("duplicate","unique"))),size=2) +
  scale_color_manual(values=c("black","#8DD3C7","#80B1D3","#FDB462")) +
  scale_alpha_manual(values=c(.5,1)) +
  labs(x = "Longitude", y = "Latitude",alpha="")


# Table 2 ----------------------------------------------------------------

# The integration results were manually validated using human coders; those
# results are reported in Table 2 of the article. For details on the procedure
# we used, please refer to online appendix D.


# Table 3 ----------------------------------------------------------------

# Statistics reported in table 3

unique_entries = meltt_data(output,columns=c("date","longitude","latitude","event_tax"))
duplicate_entries = meltt_duplicates(output,columns=c("date","longitude","latitude","event_tax")) %>% 
  mutate(id=row_number())

convert2 = function(dat,feature) dat %>% 
  select(contains(feature),id) %>% 
  rename_all(function(x){gsub(str_c(feature,"_"),"",x)}) %>% 
  select(id,dataset=dataID,event=eventID,date,longitude,latitude,event_tax) %>% 
  drop_na() %>% 
  mutate(dataset = feature) 

overlapping <- rbind(
  duplicate_entries %>% convert2("acled") %>% mutate(leader=1),
  duplicate_entries %>% convert2("ged") %>% mutate(leader=2),
  duplicate_entries %>% convert2("gtd")  %>% mutate(leader=3),
  duplicate_entries %>% convert2("scad") %>% mutate(leader=4)
) %>% mutate(match=1)



# OVERLAP BY EVENT TYPE 

# Use event taxonomy
key <- 
  event_tax %>% 
  transmute(dataset=data.source,
            event_tax = base.categories,
            etype_bin = ifelse(Level_2_text == "Violent Attack (Against Civilians)", "Violence Against Civilians",
                        ifelse(Level_2_text == "Violent Displays","Riots",
                        ifelse(Level_2_text == "Nonviolent Displays","Protests",
                        ifelse(Level_2_text %in% c("Violent Attack (No State)",
                                                   "Violent Attack",
                                                   "Violent Attack (Bombing)",
                                                   "Violent Possession"),"Violent Attack",
                                                    "Nonviolent Possession")))))
overlap_multiple <- 
  left_join(overlapping,key) %>% 
  group_by(id) %>% 
  summarize(n_datasets = n_distinct(dataset),
            etype = etype_bin[leader == min(leader)]) %>% 
  group_by(n_datasets,etype) %>% 
  count() %>% 
  ungroup

no_overlap <-
  left_join(unique_entries,key) %>% 
  left_join(.,overlapping %>% 
              select(dataset,event,match)) %>% 
  filter(is.na(match)) %>%
  group_by(etype=etype_bin) %>% 
  count() %>% 
  mutate(n_datasets=1)

# Bring together and present entries
etable<-
  bind_rows(
    no_overlap,
    overlap_multiple  
  ) %>% 
  group_by(etype) %>% 
  spread(n_datasets,n,fill = 0) %>% 
  arrange(desc(`1`))


# OVERLAP BY TIME

time_overlap <-
  overlapping %>% 
  mutate(month = lubridate::month(date),
         month_grp = ifelse(month %in% 1:3,"Jan-Mar",
                     ifelse(month %in% 4:6,"Apr-June",
                     ifelse(month %in% 7:9,"July-Sept","Oct-Dec")))) %>% 
  group_by(id) %>% 
  summarize(n_datasets = n_distinct(dataset),
            month_grp = month_grp[leader == min(leader)]) %>% 
  group_by(n_datasets,month_grp) %>% 
  count

time_no_overlap <-
  left_join(unique_entries,overlapping %>% 
              select(dataset,event,match)) %>% 
  filter(is.na(match)) %>%
  mutate(month = lubridate::month(date),
         month_grp = ifelse(month %in% 1:3,"Jan-Mar",
                     ifelse(month %in% 4:6,"Apr-June",
                     ifelse(month %in% 7:9,"July-Sept","Oct-Dec")))) %>% 
  group_by(month_grp) %>% 
  count() %>% 
  mutate(n_datasets=1)

# Bring together and present entries
ttable <-
  bind_rows(
    time_no_overlap,
    time_overlap  
  ) %>% 
  group_by(month_grp) %>% 
  spread(n_datasets,n,fill = 0) 


# OVERLAP BY REGION
map_on = function(map,dat) {st_join(st_as_sf(dat, coords = c("longitude", "latitude"),crs = 4326),map)}

space_overlap <- 
  map_on(nigeria_spatial,overlapping) %>% 
  as.data.frame() %>% 
  select(id,dataset,region,leader) %>% 
  group_by(id) %>% 
  summarize(n_datasets = n_distinct(dataset),
            region = region[leader == min(leader)]) %>% 
  group_by(n_datasets,region) %>% 
  count

space_no_overlap <- 
  left_join(unique_entries,overlapping %>% 
              select(dataset,event,match)) %>% 
  filter(is.na(match)) %>%
  map_on(nigeria_spatial,.) %>% 
  as.data.frame() %>% 
  select(dataset,region) %>% 
  mutate(region=ifelse(is.na(region),"South West",region)) %>% 
  group_by(region) %>% 
  count %>% 
  mutate(n_datasets=1)

# Bring together and present entries
stable <-
  bind_rows(
    space_no_overlap,
    space_overlap  
  ) %>% 
  group_by(region) %>% 
  spread(n_datasets,n,fill = 0) 


# Prepare and integrate data for S. Sudan 2015 --------------------------------------------------------


# Data for S. Sudan 2015 for ACLED, UCDP-GED, GTD, SCAD along with the relevant taxonomies.
load("Data/SSudan_2015.Rdata")

acled <- 
  acled %>% 
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(EVENT_DATE),
         enddate = as.Date(EVENT_DATE),
         event_tax = EVENT_TYPE,
         actor_tax = ACTOR1,
         prec_tax = GEO_PRECISION,
         descr = NOTES,
         
         # Break up riots and protests
         event_tax = ifelse(event_tax =="Riots/Protests" &  FATALITIES > 0,"Riots",
                     ifelse(event_tax =="Riots/Protests" &  FATALITIES==0,"Protests",event_tax))) %>% 
  
  # convert all column names to lower case
  rename_all(tolower) %>% 
  as.data.frame

ged <- 
  ged %>%  
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(date_start),
         enddate = as.Date(date_end),
         event_tax = type_of_violence,
         actor_tax = side_a,
         prec_tax = where_prec,
         descr = source_headline) %>% 
  as.data.frame

scad <- 
  scad %>% 
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(startdate),
         enddate = as.Date(enddate),
         event_tax = etype,
         actor_tax = actor1,
         prec_tax = locnum,
         descr = issuenote) %>% 
  as.data.frame

gtd <- 
  gtd %>%  
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(paste0(iyear,"-",imonth,"-",iday)),
         enddate = date,
         event_tax = attacktype1,
         actor_tax = gname,
         prec_tax = specificity,
         descr=summary) %>% 
  
  # Drop entries missing geo-tagged information.
  filter(!is.na(longitude)) %>% 
  as.data.frame


# (1) 3km/1day
ssudan_1 = meltt(acled,ged,gtd,scad,taxonomies=taxonomies,
                 spatwindow = 3,twindow = 1)

# (2) 50km/1day
ssudan_2 = meltt(acled,ged,gtd,scad,taxonomies=taxonomies,
                 spatwindow = 50,twindow = 1)

# Prepare and integrate data for Libya 2014 --------------------------------------------------------    


# Data for Libya 2014 for ACLED, UCDP-GED, GTD, SCAD along with the relevant taxonomies.
load("Data/Libya_2014.Rdata")


acled <- 
  acled %>% 
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(EVENT_DATE),
         enddate = as.Date(EVENT_DATE),
         event_tax = EVENT_TYPE,
         actor_tax = ACTOR1,
         prec_tax = GEO_PRECISION,
         descr = NOTES,
         
         # Break up riots and protests
         event_tax = ifelse(event_tax =="Riots/Protests" &  FATALITIES > 0,"Riots",
                     ifelse(event_tax =="Riots/Protests" &  FATALITIES==0,"Protests",event_tax))) %>% 
  
  # convert all column names to lower case
  rename_all(tolower) %>% 
  as.data.frame

ged <- 
  ged %>%  
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(date_start),
         enddate = as.Date(date_end),
         event_tax = type_of_violence,
         actor_tax = side_a,
         prec_tax = where_prec,
         descr = source_headline) %>% 
  as.data.frame

scad <- 
  scad %>% 
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(startdate),
         enddate = as.Date(enddate),
         event_tax = etype,
         actor_tax = actor1,
         prec_tax = locnum,
         descr = issuenote) %>% 
  as.data.frame

gtd <- 
  gtd %>%  
  
  # Convert dates to class date. Rename variables to match taxonomy names. 
  mutate(date = as.Date(paste0(iyear,"-",imonth,"-",iday)),
         enddate = date,
         event_tax = attacktype1,
         actor_tax = gname,
         prec_tax = specificity,
         descr=summary) %>% 
  
  # Drop entries missing geo-tagged information.
  filter(!is.na(longitude)) %>% 
  as.data.frame


# (1) 3km/1day
libya_1 = meltt(acled,ged,gtd,scad,taxonomies=taxonomies,
                spatwindow = 3,twindow = 1)

# (2) 50km/1day
libya_2 = meltt(acled,ged,gtd,scad,taxonomies=taxonomies,
                spatwindow = 50,twindow = 1)


# Table 4 ----------------------------------------------------------------

# To Validate S. Sudan integration
meltt_validate(ssudan_1,description.vars = c("descr","date","event_tax"),sample_prop=1)
meltt_validate(ssudan_2,description.vars = c("descr","date","event_tax"),sample_prop=1)

# To Validate Libya integration
meltt_validate(libya_1,description.vars = c("descr","date","event_tax"),sample_prop=.5)
meltt_validate(libya_2,description.vars = c("descr","date","event_tax"),sample_prop=.25)
