###########################################################################################
# Replication file for "Integrating Conflict Event Data"
# Authors: Donnay, Karsten, Eric Dunford, Erin McGrath, David Backer, and David Cunningham
# Journal of Conflict Resolution, accepted: March 2018
###########################################################################################
# Script: The following analysis reproduces the simulation tests reported in Table 1.
###########################################################################################


# Dependencies
require(meltt) # v.0.4.0 (to integrate data)
require(tidyverse) # v.1.2.1 (to manage data)
require(sf) # v.0.5 (for imprecision sim)
require(stringr) # v.1.3.0 (for string handling)


# Generating Simulated Data -----------------------------------------------

# The following function simulates a data generating process where multiple
# event datasets are produced with precise overlap. Given that the true
# matches are known, the production of accuracy statistics (e.g. true/false
# postive/negatives) follows accordingly.

sim.melt = function(N=100,known.events=45,known.episodes=5,
                    known.episodes.events=0,n.datasets=4,
                    t.distort=1,d.distort=2,
                    episode.window = 5,
                    min_date=as.Date("1990-01-01"),
                    max_date=as.Date("1991-01-22"),
                    lat.bounds=c(-90,90),lon.bounds=c(-180,180),
                    second.crit=data.frame(k=c(rep(1,4),rep(2,3),rep(3,4)),
                                           sc = c(10,5,3,2, 13,6,3, 15,10,5,2))){
  
  # Purpose: to simulate spatio-temporal input data to test the melt algorithm
  
  # INPUTS:
  # N: simulated sample size
  
  # known.events: # of known event-event matches
  
  # known.episodes: # of known episode-episode matches
  
  # known.episodes.events: # of known episode-event matches
  
  # n.datasets: number of simulated datasets
  
  # ambiguity: boolean -- if true, a duplicate set of points generated from 
  # the same location will be produced but with different secondary criteria. 
  # This is done to simulate a unique observation that is cotemporaneous with
  # another known observation to test disambiguation. 
  
  # t.distort: time distortion around the known values
  
  # d.distort: spatial distortion around the known values (in km)
  
  # min_date: minimum date range
  
  # max_date: maximum date range
  
  # lat.bounds: latitude boundaries 
  
  # lon.bounds: longitude boundaries
  
  # second.crit: vector of max values of the secondary criteria. Entered as
  # matrix (column 1 = k, column 2 = tax levels)
  
  # Output: (as list)
  # data.frame of simulated data (ordered by date)
  # matrix of known matches (matches are recorded by row)
  # emulated data and taxonomies for input into MELTT
  
  # Essential Checks ----- 
  if (known.episodes > 0 | known.episodes.events > 0){
    if (max_date-min_date < episode.window){
      stop("Date range cannot be less than the episode.window when simulating episodes.")
    }
  }
  if (known.events == 0 & known.episodes == 0 & known.episodes.events == 0 ){
    stop("No known events specified!")
  }
  if ( ((known.events+known.episodes+known.episodes.events)*n.datasets) >= N ){
    stop("Insufficient N assigned given the number of known matches and datasets.")
  }
  
  # Boundary
  boundaries = data.frame(lon=lon.bounds,lat=lat.bounds)
  
  # Secondary criteria
  sets = second.crit[,2]
  
  # GENERATE taxonomies that are internally consisent.
  K = unique(second.crit[,1])
  taxonomies = list() # Used in simulation
  out.taxonomies = list() # Output
  for (k in K){
    sub = second.crit[second.crit[,1]==k,2]
    M = matrix(1:sub[1])
    if (length(sub) > 1){
      for (i in 2:length(sub)){
        M = cbind(M,round(runif(nrow(M),min = 1,sub[i])))
      }
    }
    taxonomies[[paste0("var_",k)]] = M
    # Create output taxonomies
    M2 = c()
    for(d in 1:n.datasets){
      data.source = paste0("D",d)
      base.categories = paste0(data.source,"_v_",M[,1])
      M2 = rbind(M2,cbind(data.source,base.categories,as.data.frame(M)))
    }
    out.taxonomies[[paste0("var_",k)]] = M2
  }
  
  # CREATE event-event ----------------------------
  if (known.events > 0){
    # Identifying the "known" points
    lon.locs = runif(known.events,boundaries[1,1],boundaries[2,1])
    lat.locs = runif(known.events,boundaries[1,2],boundaries[2,2])
    known.locs = data.frame(lon=lon.locs,lat=lat.locs)
    
    # Identifying the temporal window
    window = 1:as.integer(max_date - min_date)
    max.window = max(window)
    known.locs$date = sample(window,known.events,replace=T)
    known.locs$enddate = known.locs$date
    
    # Range of the temporal distortion on the known observations
    if(t.distort>0){
      t.distort = c((t.distort:1)*-1,0,1:t.distort)
    }
    
    # Calculating points that are dependent on the known points (dataset 1 is
    # always the known)
    output = NULL
    for (i in 1:known.events){
      dist <- rbind(known.locs[i,1:2],geosphere::destPoint(known.locs[i,1:2], 
                                                           b=sample(c(0,90,180,270),n.datasets-1,replace=TRUE),
                                                           d=runif(n.datasets-1,0,d.distort)*1000/sqrt(2)))
      date <- c(known.locs[i,3],known.locs[i,3] + sample(t.distort,n.datasets-1,replace=T))
      date = ifelse(date<=0,1,date) # Make sure there are no 0 or negative date terms (due to time distortion)
      enddate <- date
      
      # Taxonomy levels....
      crit = c()
      for (k in K){
        v = taxonomies[[k]][taxonomies[[k]][,1]==sample(taxonomies[[k]][,1],1),]
        v2 = rep(v,n.datasets)
        temp.crit = t(matrix(v2,ncol=n.datasets,nrow=ncol(taxonomies[[k]])))
        crit = cbind(crit,temp.crit)
      }
      temp = data.frame(dataset=1:n.datasets,dist,date=date,enddate=enddate,crit,match=i)
      output = rbind(output,temp)
    }
  }
  
  
  # CREATE episode-episode ----------------------------
  if (known.episodes > 0){ # Episodes-episodes known matches
    # Generate episodal known Points
    ep.lon.locs = runif(known.episodes,boundaries[1,1],boundaries[2,1])
    ep.lat.locs = runif(known.episodes,boundaries[1,2],boundaries[2,2])
    ep.known.locs = data.frame(lon=ep.lon.locs,lat=ep.lat.locs)
    ep.known.locs$date = sample(window,known.episodes,replace=T)
    # Generate a temporal window: random from 1 to specified window
    ep.known.locs$enddate = round((ep.known.locs$date + runif(known.episodes,1,episode.window)))
    # Episode cannot exceed the specified temporal window 
    ep.known.locs$enddate = ifelse(ep.known.locs$enddate > max.window,max.window,ep.known.locs$enddate)
    
    # Calculating points (episodes) that are dependent on the known points
    # (dataset 1 is always the known) 
    ep.output = NULL
    for(i in 1:known.episodes){
      dist <- rbind(ep.known.locs[i,1:2],geosphere::destPoint(ep.known.locs[i,1:2], 
                                                              b=sample(c(0,90,180,270),n.datasets-1,replace=TRUE),
                                                              d=runif(n.datasets-1,0,d.distort)*1000/sqrt(2)))
      date <- c(ep.known.locs[i,3],ep.known.locs[i,3] + sample(t.distort,n.datasets-1,replace=T))
      date = ifelse(date<=0,1,date) # Make sure there are no 0 or negative date terms (due to time distortion)
      enddate <- c(ep.known.locs[i,4],ep.known.locs[i,4] + sample(t.distort,n.datasets-1,replace=T))
      enddate = ifelse(enddate<=date,date+1,enddate) # ensure the enddate is at least one day more than the start date. 
      enddate = ifelse(enddate > max.window,max.window,enddate) # make sure the event does not exceed the maximum window
      
      # Taxonomy levels....
      crit = c()
      for (k in K){
        v = taxonomies[[k]][taxonomies[[k]][,1]==sample(taxonomies[[k]][,1],1),]
        v2 = rep(v,n.datasets)
        temp.crit = t(matrix(v2,ncol=n.datasets,nrow=ncol(taxonomies[[k]])))
        
        crit = cbind(crit,temp.crit)
      }
      temp = data.frame(dataset=1:n.datasets,dist,date=date,enddate=enddate,crit,match=i+30000)
      # For simplicity, I build the ambiguity feature out of the episodal portion. 
      ep.output = rbind(ep.output,temp)
    }
  }
  
  # CREATE episode-event ----------------------------
  if (known.episodes.events > 0){ # Episodes-Event Known Matches
    # Generate episodal known Points
    epev.lon.locs = runif(known.episodes.events,boundaries[1,1],boundaries[2,1])
    epev.lat.locs = runif(known.episodes.events,boundaries[1,2],boundaries[2,2])
    epev.known.locs = data.frame(lon=epev.lon.locs,lat=epev.lat.locs)
    epev.known.locs$date = sample(window,known.episodes.events,replace=T)
    # Generate a temporal window: random from 1 to specified window
    epev.known.locs$enddate = round((epev.known.locs$date + runif(known.episodes.events,1,episode.window)))
    # Episode cannot exceed the specified temporal window 
    epev.known.locs$enddate = ifelse(epev.known.locs$enddate > max.window,max.window,epev.known.locs$enddate)
    
    
    # Calculating points that are dependent on the known points (dataset 1 is
    # always the known) in the episode. 
    epev.output = NULL
    for (i in 1:known.episodes.events){
      dist <- rbind(epev.known.locs[i,1:2],geosphere::destPoint(epev.known.locs[i,1:2], 
                                                                b=sample(c(0,90,180,270),n.datasets-1,replace=TRUE),
                                                                d=runif(n.datasets-1,0,d.distort)*1000/sqrt(2)))
      # Range of dates that event can occur
      date.range = range(epev.known.locs[i,3],epev.known.locs[i,4])
      gen.dates = sample(date.range[1]:date.range[2],n.datasets-1,replace=T)
      date <- c(epev.known.locs[i,3],gen.dates)
      date = ifelse(date<=0,1,date) # Make sure there are no 0 or negative date terms (due to time distortion)
      enddate <- c(epev.known.locs[i,4],gen.dates)
      enddate = ifelse(enddate > max.window,max.window,enddate) # make sure the event does not exceed the maximum window
      
      # Taxonomy levels....
      crit = c()
      for (k in K){
        v = taxonomies[[k]][taxonomies[[k]][,1]==sample(taxonomies[[k]][,1],1),]
        v2 = rep(v,n.datasets)
        temp.crit = t(matrix(v2,ncol=n.datasets,nrow=ncol(taxonomies[[k]])))
        
        crit = cbind(crit,temp.crit)
      }
      
      temp = data.frame(dataset=1:n.datasets,dist,date=date,enddate=enddate,crit,match=i+40000)
      # For simplicity, I build the ambiguity feature out of the episodal portion. 
      epev.output = rbind(epev.output,temp)
    }
  }
  
  
  # random points ---------
  # Calculating random points that are dependent on a randomly selected point within the boundary
  if(known.events > 0){ N = N - (known.events*n.datasets) }  # Adjusting the number of random observations that need to be calculated
  if(known.episodes > 0 ){ N = N - (known.episodes*n.datasets) } # Adjustment for episodes-episodes
  if(known.episodes.events > 0 ){ N = N - (known.episodes.events*n.datasets) } # Adjustment for episodes-events
  # Selecting Random locations from within the boundary
  r.locs = data.frame(lon = runif(N,boundaries[1,1],boundaries[2,1]),
                      lat = runif(N,boundaries[1,2],boundaries[2,2])) # scatter points around the random point
  r.locs$date <- sample(window,N,replace = T) # scatter temporal locations given the window range
  if (known.episodes > 0 | known.episodes.events > 0){
    # Assign Enddates
    draw = as.numeric(sample(row.names(r.locs),size = round(nrow(r.locs)/2),replace = F))
    # Events 
    r.locs[draw,"enddate"] <- r.locs[draw,"date"]
    # Episodes
    r.locs[draw*-1,"enddate"] <- round((r.locs[draw*-1,"date"] + runif(length(r.locs[draw*-1,"date"]),1,episode.window)))
    r.locs[draw*-1,"enddate"] = ifelse(r.locs[draw*-1,"enddate"] > max.window,max.window,r.locs[draw*-1,"enddate"])
  }else{
    r.locs$enddate <- r.locs$date
  }
  
  r.locs$dataset <- sample(1:n.datasets,N,replace = T) # assign a "dataset"
  
  # Taxonomy levels....
  crit = c()
  for (k in K){
    v = taxonomies[[k]][sample(taxonomies[[k]][,1],N,replace = T),]
    crit = cbind(crit,v)
  }
  r.locs2 = data.frame(r.locs,crit);r.locs2$match=0
  
  # Bringing random and dependent together
  ev = c()
  ep = c()
  epev = c()
  if (known.events > 0){ev = output}
  if (known.episodes > 0){ep = ep.output}
  if (known.episodes.events > 0){epev = epev.output}
  
  full.data = rbind(ev,ep,epev,r.locs2)
  ordered.out = NULL
  for(i in 1:n.datasets){
    f = full.data[full.data$dataset==i,]
    f = f[order(f$date),]; f$event <- 1:nrow(f)
    f = cbind(f[,c("dataset","event","date","enddate","lat","lon","match")],f[,c(grep("X",colnames(f)))]) 
    ordered.out = rbind(ordered.out,f)
  }
  
  # Identify the Known matches and report them
  ordered.out$id = paste(ordered.out$dataset,ordered.out$event,sep="-")
  true.matches = list()
  detect.known = unique(ordered.out$match[ordered.out$match!=0])
  true.matches = as.data.frame(t(apply(matrix(detect.known),1,function(x)ordered.out[ordered.out$match==x,"id"])),stringsAsFactors = F)
  
  known.set = rep(NA,length(detect.known)) 
  known.set[which(detect.known>=30000 & detect.known < 40000)] = "episode-episode"
  known.set[which(detect.known>=40000)] = "episode-event"
  known.set[which(detect.known<30000)] = "event-event"
  true.matches$type = known.set
  
  # Convert Sim Data to emulate real data -----------
  user.data = ordered.out[,c("dataset","event","date","enddate","lat","lon")]
  
  # Recover Date Range 
  date.range = seq(min_date,max_date,by="day")
  user.data$date = date.range[ordered.out$date]
  user.data$enddate = date.range[ordered.out$enddate]
  # Rename Columns
  colnames(user.data)[5] = "latitude"
  colnames(user.data)[6] = "longitude"
  # assign dataset names a unique character identifier
  user.data$dataset = paste0("D",user.data$dataset)
  
  s = unique(second.crit[,1])
  s2 = paste0("var_",s) # variable place holders
  taxs = ordered.out[,grep("X",colnames(ordered.out))]
  datasets = ordered.out[,grep("data",colnames(ordered.out))]
  for(i in s){
    sub = second.crit[second.crit[,1]==i,]
    # Subset the taxonomy
    s.taxs = as.matrix(taxs[,1:nrow(sub)])
    var = paste0("v_",s.taxs[,1]) # create a base.category using one taxonomy layer
    user.data[,s2[i]] = var # map it onto the user data
    user.data[,s2[i]] = paste(user.data$dataset,user.data[,s2[i]],sep="_") # make unique to each dataset
    taxs = taxs[,1:nrow(sub)*-1] # clean the used taxonomies
  }
  
  # Return true.matches and output frame
  out.all = list()
  # Return emulated data
  out.all[["user.data"]] <- user.data
  out.all[["user.taxonomies"]] <- out.taxonomies
  out.all[["true.matches"]] <- true.matches
  return(out.all)
}




# Aux. Functions -------------------------

# The following function streamlines the handling the simulated data.

process = function(out,sim.data,show_sets=F){
  # SET UP ---------------
  known = sim.data$true.matches
  known2 = known[,grep("type",colnames(known))*-1] # Remove Type from Known
  N = nrow(sim.data$user.data) # no. of observations
  
  
  # GET the matched list from MELTT in comparable order ---------
  matched = meltt_duplicates(out,columns = "")
  if(nrow(matched)==0|length(matched)==0){
    stop("\n\nNo Matches were Located!")
}
  c = 1:ncol(matched) * 1:ncol(matched) %% 2
  c = c[c!=0]
  m1 = matched[,c]
  m2 = matched[,c*-1]
  out = c()
  for(i in 1:length(c)){
    out=cbind(out,paste(m1[,i],m2[,i],sep="-"))  
}
  matches = data.frame(out,stringsAsFactors = F)
  
  # ASSESS performance --------------
  # TRUE POSITIVES
  located = matrix(NA,nrow=nrow(known2),ncol=ncol(known2))
  for (i in 1:nrow(known2)){
    for (j in 1:ncol(known2)){
      located[i,j] = any(known2[i,j] == matches)
    }
  }
  true.positives = sum(located)/(nrow(known2)*ncol(known2))
  # FALSE POSITIVES 
  false_positives = matrix(NA,nrow=nrow(matches),ncol=ncol(matches))
  for (i in 1:nrow(matches)){
    for (j in 1:ncol(matches)){
      false_positives[i,j] = !any(matches[i,j] == known2) & matches[i,j] !="0-0"
    }
  }
  false.positives = sum(false_positives)/(N-nrow(known2))
  
  # RETURN
  full.out = list()
  full.out[["true.pos"]] = true.positives
  full.out[["false.pos"]] = false.positives
  if (show_sets){
    full.out[["located.matches"]] = matches
    full.out[["known.matches"]] = known
  }
  return(full.out)
}


expand = function(sim.data){
  # Expands out the simulated data sets and assigns them to the global
  # environment
  types = unique(sim.data$user.data$dataset)
  df =  sim.data$user.data
  for(t in types){
    assign(t,df[df$dataset==t,-2],envir = globalenv())
  }
  cat("Data in Global Environment:",types)
}


# Simulations Correspond with Entry 1 in Table 1-------------------

# Conditions: flat taxonomy and little distance/time distortion (+-2km/+-1day) in matching
# entries; examine the influence of the order of the data on the respective output.

# Generate the data

# Secondary Criteria
sc = data.frame(k=c(rep(1,2),rep(2,2),rep(3,2)),
                sc = c(9,6,5,2,5,3))

order = gtools::permutations(n = 4, r = 4, v = c("D1","D2","D3","D4"))

# Simulate Data
set.seed(123)
sim.data = sim.melt(N=1000,known.events = 45,known.episodes = 5,
                    known.episodes.events = 0,
                    n.datasets = 4,t.distort = 1,d.distort = 2,
                    episode.window = 5,
                    min_date = as.Date("2011-01-01"),max_date =as.Date("2012-01-01"),
                    lat.bounds = c(4.5,14),lon.bounds =c(3,14),
                    second.crit = sc)
expand(sim.data) # Map data onto the environment.

# Permute the order of the data and integrate
# - select fuzziness parameters that can account for full uncertainty
test = list(); twindow = 2; spatwindow = 4
test[[1]] = meltt(D1,D2,D3,D4,taxonomies = sim.data$user.taxonomies,
                  twindow = twindow,spatwindow = spatwindow)
test[[2]] = meltt(D1,D2,D4,D3,taxonomies = sim.data$user.taxonomies,
                  twindow = twindow,spatwindow = spatwindow)
test[[3]]  = meltt(D1,D3,D2,D4,taxonomies = sim.data$user.taxonomies,
                   twindow = twindow,spatwindow = spatwindow)
test[[4]]  = meltt(D1,D3,D4,D2,taxonomies = sim.data$user.taxonomies,
                   twindow = twindow,spatwindow = spatwindow)
test[[5]]  = meltt(D1,D4,D2,D3,taxonomies = sim.data$user.taxonomies,
                   twindow = twindow,spatwindow = spatwindow)
test[[6]]  = meltt(D1,D4,D3,D2,taxonomies = sim.data$user.taxonomies,
                   twindow = twindow,spatwindow = spatwindow)
test[[7]]  = meltt(D2,D1,D3,D4,taxonomies = sim.data$user.taxonomies,
                   twindow = twindow,spatwindow = spatwindow)
test[[8]]  = meltt(D2,D1,D4,D3,taxonomies = sim.data$user.taxonomies,
                   twindow = twindow,spatwindow = spatwindow)
test[[9]]  = meltt(D2,D3,D1,D4,taxonomies = sim.data$user.taxonomies,
                   twindow = twindow,spatwindow = spatwindow)
test[[10]]  = meltt(D2,D3,D4,D1,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[11]]  = meltt(D2,D4,D1,D3,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[12]]  = meltt(D2,D4,D3,D1,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[13]]  = meltt(D3,D1,D2,D4,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[14]]  = meltt(D3,D1,D4,D2,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[15]]  = meltt(D3,D2,D1,D4,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[16]]  = meltt(D3,D2,D4,D1,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[17]]  = meltt(D3,D4,D1,D2,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[18]]  = meltt(D3,D4,D2,D1,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[19]]  = meltt(D4,D1,D2,D3,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[20]]  = meltt(D4,D1,D3,D2,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[21]]  = meltt(D4,D2,D1,D3,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[22]]  = meltt(D4,D2,D3,D1,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[23]]  = meltt(D4,D3,D1,D2,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)
test[[24]]  = meltt(D4,D3,D2,D1,taxonomies = sim.data$user.taxonomies,
                    twindow = twindow,spatwindow = spatwindow)

true = sim.data$true.matches
no.order = gtools::permutations(n = 4, r = 4, v = 1:4)
entry1 = list()
true.pos = c()
false.pos = c()
for(j in 1:24){
  # Need to adjust the ordering in the known set to reflect the permutation
  true.temp = true[,c(no.order[j,],5)]
  true.temp[,1] = gsub(paste0(no.order[j,1],"-"),"1-",true.temp[,1])
  true.temp[,2] = gsub(paste0(no.order[j,2],"-"),"2-",true.temp[,2])
  true.temp[,3] = gsub(paste0(no.order[j,3],"-"),"3-",true.temp[,3])
  true.temp[,4] = gsub(paste0(no.order[j,4],"-"),"4-",true.temp[,4])
  holder = sim.data
  holder$true.matches = true.temp
  ppp = process(test[[j]],holder)
  true.pos = c(true.pos,ppp$true.pos)
  false.pos = c(false.pos,ppp$false.pos)
}
entry1[["Manipulation"]] <- "Order of the Datasets"
entry1[["manipulated.parameter"]] <- order
entry1[["true.pos"]] <- true.pos
entry1[["false.pos"]] <- false.pos


# Simulations Correspond with Entries 2 - 4 in Table 1 -------------------

# Conditions: flat taxonomy and little distance/time distortion (+-2km/+-1day); 
# MELTT: expand the (a) distance fuzziness, (b) temporal fuzziness, (c) both concurrently.

# Generate the data

# Secondary Criteria
sc = data.frame(k=c(rep(1,2),rep(2,2),rep(3,2)),
                sc = c(9,6,5,2,5,3))
# Simulate Data
set.seed(123)
sim.data = sim.melt(N=1000,known.events = 45,known.episodes = 5,
                    known.episodes.events = 0,
                    n.datasets = 4,t.distort = 1,d.distort = 1,
                    episode.window = 5,
                    min_date = as.Date("2011-01-01"),max_date =as.Date("2012-01-01"),
                    lat.bounds = c(4.5,14),lon.bounds =c(3,14),
                    second.crit = sc)
expand(sim.data) # Map data onto the environment. 

# (a) distance fuzziness 
entry2 = list()
true.pos = c()
false.pos = c()
distance = seq(2,21,by=1)
for(i in 1:length(distance)){
  cat(paste0("\nMANIPULATION ",i,"\n"))
  # Run MELT
  out = meltt(D1,D2,D3,D4,taxonomies = sim.data$user.taxonomies,
              twindow = 2,spatwindow = distance[i])
  ppp = process(out,sim.data,show_sets = FALSE)
  true.pos = c(true.pos,ppp$true.pos)
  false.pos = c(false.pos,ppp$false.pos)
}
entry2[["Manipulation"]] <- "SIM: distance distortion +-1km, temporal distortion +-1day; MELTT: expand
        manipulating the distance fuzziness by 1km, time == 2."
entry2[["manipulated.parameter"]] <- distance
entry2[["true.pos"]] <- true.pos
entry2[["false.pos"]] <- false.pos

# (b) temporal fuziness
entry3 = list()
true.pos = c()
false.pos = c()
distance = seq(2,21,by=1)
for(i in 1:length(distance)){
  cat(paste0("\nMANIPULATION ",i,"\n"))
  
  # Run MELTT
  out = meltt(D1,D2,D3,D4,taxonomies = sim.data$user.taxonomies,
              twindow = distance[i],spatwindow = 2)
  ppp = process(out,sim.data,show_sets = FALSE)
  true.pos = c(true.pos,ppp$true.pos)
  false.pos = c(false.pos,ppp$false.pos)
}
entry3[["Manipulation"]] <- "SIM: distance distortion +-1km, temporal distortion +-1km; MELTT: expand
        manipulating the temporal fuzziness (by day), space == 2. "
entry3[["manipulated.parameter"]] <- distance
entry3[["true.pos"]] <- true.pos
entry3[["false.pos"]] <- false.pos

# (c) both spatial and temporal fuzziness
entry4 = list()
true.pos = c()
false.pos = c()
distance = seq(2,21,by=1)
for(i in 1:length(distance)){
  cat(paste0("\nMANIPULATION ",i,"\n"))
  out = meltt(D1,D2,D3,D4,taxonomies = sim.data$user.taxonomies,
              twindow = distance[i],spatwindow = distance[i])
  ppp = process(out,sim.data,show_sets = FALSE)
  true.pos = c(true.pos,ppp$true.pos)
  false.pos = c(false.pos,ppp$false.pos)
}
entry4[["Manipulation"]] <- "SIM: distance distortion +-1km, temporal distortion +-1day; MELTT: expand
        manipulating both (sp/time) by units of 1."
entry4[["manipulated.parameter"]] <- distance
entry4[["true.pos"]] <- true.pos
entry4[["false.pos"]] <- false.pos



# Simulations Correspond with Entry 5 in Table 1 ------------------------------------

# Checking the dependence on the number of input data

sc = data.frame(k=c(rep(1,2),rep(2,2),rep(3,2)),
                sc = c(9,6,5,2,5,3))
N.DATA = 3:10
entry5 = list()
true.pos = c()
false.pos = c()
for(i in N.DATA){
  cat(paste0("\nMANIPULATION ",i-2,"\n"))
  # Simulate Data
  set.seed(123)
  sim.data = sim.melt(N=1000,known.events = 45,known.episodes = 5,
                      known.episodes.events = 0,
                      n.datasets = i,t.distort = 1,d.distort = 2,
                      episode.window = 5,
                      min_date = as.Date("2011-01-01"),max_date =as.Date("2012-01-01"),
                      lat.bounds = c(4.5,14),lon.bounds =c(3,14),
                      second.crit = sc)
  expand(sim.data) # Map data onto the environment. 
  
  # Run MELT
  # - select fuzziness parameters that can account for full uncertainty
  twindow = 2; spatwindow = 4
  if(i==3){out = meltt(D1,D2,D3,taxonomies = sim.data$user.taxonomies,twindow = twindow,spatwindow = spatwindow)}
  if(i==4){out = meltt(D1,D2,D3,D4,taxonomies = sim.data$user.taxonomies,twindow = twindow,spatwindow = spatwindow)}
  if(i==5){out = meltt(D1,D2,D3,D4,D5,taxonomies = sim.data$user.taxonomies,twindow = twindow,spatwindow = spatwindow)}
  if(i==6){out = meltt(D1,D2,D3,D4,D5,D6,taxonomies = sim.data$user.taxonomies,twindow = twindow,spatwindow = spatwindow)}
  if(i==7){out = meltt(D1,D2,D3,D4,D5,D6,D7,taxonomies = sim.data$user.taxonomies,twindow = twindow,spatwindow = spatwindow)}
  if(i==8){out = meltt(D1,D2,D3,D4,D5,D6,D7,D8,taxonomies = sim.data$user.taxonomies,twindow = twindow,spatwindow = spatwindow)}
  if(i==9){out = meltt(D1,D2,D3,D4,D5,D6,D7,D8,D9,taxonomies = sim.data$user.taxonomies,twindow = twindow,spatwindow = spatwindow)}
  if(i==10){out = meltt(D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,taxonomies = sim.data$user.taxonomies,twindow = twindow,spatwindow = spatwindow)}
  ppp = process(out,sim.data,show_sets = FALSE)
  true.pos = c(true.pos,ppp$true.pos)
  false.pos = c(false.pos,ppp$false.pos)
}
entry5[["Manipulation"]] <- "Impact of Increasing the Number of Input Dataframes"
entry5[["manipulated.parameter"]] <- N.DATA
entry5[["true.pos"]] <- true.pos
entry5[["false.pos"]] <- false.pos




# Simulations Correspond with Entry 6 in Table 1 ---------------------

# Condition: little distance/small time distortion (+-2km/+-1day), varying
# taxonomy (many to few secondary criteria)

# Taxonomy Manipulation 
number.of.criteria = list()
sc = data.frame(k=1,sc=c(9,6,2))
number.of.criteria[[1]] = sc
for(i in 2:10){
  sc = data.frame(k=1,sc=c(9,6,2))
  for(j  in 2:i){
    sc = rbind(sc,data.frame(k=j,sc=c(9,6,2)))
  }
  number.of.criteria[[i]] = sc
}

entry6 = list()
true.pos = c()
false.pos = c()
for(i in 1:10){
  cat(paste0("\nMANIPULATION ",i,"\n"))
  # Simulate Data
  set.seed(123)
  sim.data = sim.melt(N=1000,known.events = 45,known.episodes = 5,
                      known.episodes.events = 0,
                      n.datasets = 4,t.distort = 1,d.distort = 2,
                      episode.window = 5,
                      min_date = as.Date("2011-01-01"),max_date =as.Date("2012-01-01"),
                      lat.bounds = c(4.5,14),lon.bounds =c(3,14),
                      second.crit = number.of.criteria[[i]])
  expand(sim.data) # Map data onto the environment. 
  
  # Run MELT
  # - select fuzziness parameters that can account for full uncertainty
  out = meltt(D1,D2,D3,D4,taxonomies = sim.data$user.taxonomies,
              twindow = 2,spatwindow = 4)
  ppp = process(out,sim.data,show_sets = FALSE)
  true.pos = c(true.pos,ppp$true.pos)
  false.pos = c(false.pos,ppp$false.pos)
}
entry6[["Manipulation"]] <- "SIM: little distance and time distortion (+-2km/+-1day), varying taxonomy (few (1) to many (10)
    secondary criteria -- all the same depth: 9,6,2); MELTT: fuzzy (4km/2days)"
entry6[["manipulated.parameter"]] <- number.of.criteria
entry6[["true.pos"]] <- true.pos
entry6[["false.pos"]] <- false.pos


# Simulations Correspond with Entry 7 in Table 1 -------------------------------------

# Condition: little distance/time uncertainty in matches (+-2km/+-1day),
# varying taxonomy depth (increasing the depth in the single taxonomy from 2 to 10)

# Taxonomy manipulation
r = 3:11
criteria.depth = list()
for(i in 1:9){
  sc = data.frame(k=1,sc=c(r[i],2))
  criteria.depth[[i]] = sc
} 

entry7 = list()
true.pos = c()
false.pos = c()
for(i in 1:9){
  cat(paste0("\nMANIPULATION ",i,"\n"))
  # Simulate Data
  set.seed(123)
  sim.data = sim.melt(N=1000,known.events = 45,known.episodes = 5,
                      known.episodes.events = 0,
                      n.datasets = 4,t.distort = 1,d.distort = 2,
                      episode.window = 5,
                      min_date = as.Date("2011-01-01"),max_date =as.Date("2012-01-01"),
                      lat.bounds = c(4.5,14),lon.bounds =c(3,14),
                      second.crit = criteria.depth[[i]])
  expand(sim.data) # Map data onto the environment. 
  
  # Run MELT
  # - select fuzziness parameters that can account for full uncertainty
  out = meltt(D1,D2,D3,D4,taxonomies = sim.data$user.taxonomies,
              twindow = 2,spatwindow = 4)
  ppp = process(out,sim.data,show_sets = FALSE)
  true.pos = c(true.pos,ppp$true.pos)
  false.pos = c(false.pos,ppp$false.pos)
}
entry7[["Manipulation"]] <- "SIM: little distance and time distortion (+-2km/+-1day), varying taxonomy depth
    (increasing the depth in the single taxonomy from 2 to 10): MELTT: fuzzy (4km/2days)"
entry7[["manipulated.parameter"]] <- number.of.criteria
entry7[["true.pos"]] <- true.pos
entry7[["false.pos"]] <- false.pos


# Simulations Correspond with Entries 8 and 9 in Table 1 ------------------------------------
# Simulating Urban Concentration - many events, tight spatial concentration
# with ambiguity, flat taxonomy, and occuring close in time.

# (A) Strong Taxonomy (entry 8)
# (B) Weak Taxonomy (entry 9)

# calculate lat/lon bounds such that they correspond to given density (true
# near equator and for not too large distances)
N = 1000
density = exp(seq(-9,2,0.1))
distance = unlist(lapply(1:length(density),
                         function(x) 1/2*sqrt(N/(density[x]*111.111^2))))
sc = data.frame(k=c(rep(1,2),rep(2,2),rep(3,2)),
                sc = c(9,6,5,2,5,3))

entry8a = list() # Smartmatch
entry8b = list() # Exactmatch
for(t in 1:2){
  if(t==1){
    sm = T
  }else{sm = F}
  true.pos = c()
  false.pos = c()
  for(i in 1:length(distance)){
    cat(paste0("\nMANIPULATION ",i,"\n"))
    # Simulate Data
    set.seed(123)
    sim.data = sim.melt(N=1000,known.events = 45,known.episodes = 5,
                        known.episodes.events = 0,
                        n.datasets = 4,t.distort = 1,d.distort = 2,
                        episode.window = 5,
                        min_date = as.Date("2011-01-01"),max_date =as.Date("2012-01-01"),
                        lat.bounds = c(-distance[i],distance[i]),
                        lon.bounds= c(-distance[i],distance[i]),
                        second.crit = sc)
    expand(sim.data) # Map data onto the environment. 
    
    # Run MELT
    if(sm){
      smartmatch = T
      certainty = c(0,0,0)
    } else{
      smartmatch = F
      certainty = c(1,1,1)
    }
    out = meltt(D1,D2,D3,D4,taxonomies = sim.data$user.taxonomies,
                twindow = 2,spatwindow = 4,smartmatch = smartmatch,
                certainty = certainty)
    ppp = process(out,sim.data,show_sets = F)
    true.pos = c(true.pos,ppp$true.pos)
    false.pos = c(false.pos,ppp$false.pos)
  }
  if(sm){
    entry8a[["Manipulation"]] <- "Simulating Urban Concentration - many events, tight spatial 
        concentration (expanding lon-lat from +/-.01 to +/-2 by .01), 
        MELTT: fuzzy (4km/2day) - deep taxonomy, smartmatch"
    entry8a[["manipulated.parameter"]] <- distance
    entry8a[["true.pos"]] <- true.pos
    entry8a[["false.pos"]] <- false.pos
  } else {
    entry8b[["Manipulation"]] <- "Simulating Urban Concentration - many events, tight spatial 
        concentration (expanding lon-lat from +/-.01 to +/-2 by .01), 
        MELTT: fuzzy (4km/2day) - deep taxonomy, exact match"
    entry8b[["manipulated.parameter"]] <- distance
    entry8b[["true.pos"]] <- true.pos
    entry8b[["false.pos"]] <- false.pos
  }
}



# Shallow Tax
sc = data.frame(k=c(rep(1,2)),
                sc=c(2,1))

entry9a = list()
entry9b = list()
for(t in 1:2){
  if(t==1){
    sm = T
  }else{sm = F}
  true.pos = c()
  false.pos = c()
  for(i in 1:length(distance)){
    cat(paste0("\nMANIPULATION ",i,"\n"))
    # Simulate Data
    set.seed(123)
    sim.data = sim.melt(N=1000,known.events = 45,known.episodes = 5,
                        known.episodes.events = 0,
                        n.datasets = 4,t.distort = 1,d.distort = 2,
                        episode.window = 5,
                        min_date = as.Date("2011-01-01"),max_date =as.Date("2012-01-01"),
                        lat.bounds = c(-distance[i],distance[i]),
                        lon.bounds= c(-distance[i],distance[i]),
                        second.crit = sc)
    expand(sim.data) # Map data onto the environment. 
    
    # Run MELT
    if(sm){
      smartmatch = T
      certainty = c(0,0)
    } else{
      smartmatch = F
      certainty = c(1)
    }
    out = meltt(D1,D2,D3,D4,taxonomies = sim.data$user.taxonomies,
                twindow = 2,spatwindow = 4,smartmatch = smartmatch,
                certainty = certainty)
    ppp = process(out,sim.data,show_sets = FALSE)
    true.pos = c(true.pos,ppp$true.pos)
    false.pos = c(false.pos,ppp$false.pos)
  }
  if(sm){
    entry9a[["Manipulation"]] <- "Simulating Urban Concentration - many events, tight spatial 
        concentration (expanding lon-lat from +/-.01 to +/-2 by .01), 
        MELTT: fuzzy (4km/1day) - shallow taxonomy, smartmatch"
    entry9a[["manipulated.parameter"]] <- distance
    entry9a[["true.pos"]] <- true.pos
    entry9a[["false.pos"]] <- false.pos
  } else {
    entry9b[["Manipulation"]] <- "Simulating Urban Concentration - many events, tight spatial 
        concentration (expanding lon-lat from +/-.01 to +/-2 by .01), 
        MELTT: fuzzy (4km/2day) - shallow taxonomy, exact match"
    entry9b[["manipulated.parameter"]] <- distance
    entry9b[["true.pos"]] <- true.pos
    entry9b[["false.pos"]] <- false.pos
  }
}



# Simulations Correspond with Entry 10 in Table 1 ------------------------

# The level of geo-spatial precision as a proportion of the total data.
# Distortions both w/r/t levels and datasets are randomly assigned, so
# multiple simulations are run and the average is then reported. Finally,
# different spatial windows are explored to assess the rememdy. 

# Simulate Data 

# Secondary Criteria
sc = data.frame(k=c(rep(1,2),rep(2,2),rep(3,2)),
                sc = c(9,6,5,2,5,3))
# Simulate Data
set.seed(123)
sim.data = sim.melt(N=1000,known.events = 45,known.episodes = 5,
                    known.episodes.events = 0,
                    n.datasets = 4,t.distort = 1,d.distort = 2,
                    episode.window = 5,
                    min_date = as.Date("2011-01-01"),max_date =as.Date("2012-01-01"),
                    lat.bounds = c(4.5,14),lon.bounds =c(3,14),
                    second.crit = sc)

D = sim.data$user.data %>% group_by(dataset) %>% nest() 
tax = sim.data$user.taxonomies


# Simulate Geographies
D_geo = D %>% unnest %>% st_as_sf(coords=c("latitude","longitude")) 

level4 = st_make_grid(D_geo,cellsize = 12) %>% 
  as.tibble %>% 
  mutate(level4 = st_centroid(geometry)) 

level3 = st_make_grid(D_geo,cellsize = 8) %>% as.tibble %>% 
  mutate(level3 = st_centroid(geometry)) 

level2 = st_make_grid(D_geo,cellsize = 2) %>% as.tibble %>% 
  mutate(level2 = st_centroid(geometry)) 

# Gather
x1 = st_intersection(D_geo,st_as_sf(level4)) %>% st_set_geometry(NULL)  %>% select(dataset,event,contains("level")) 
x2 = st_intersection(D_geo,st_as_sf(level3)) %>% st_set_geometry(NULL)  %>% select(dataset,event,contains("level")) 
x3 = st_intersection(D_geo,st_as_sf(level2)) %>% st_set_geometry(NULL)  %>% select(dataset,event,contains("level")) 

# Build coordinate system options upfront
D_geo2 = left_join(D_geo,x1) %>% left_join(x2) %>% left_join(x3)
coords <-
  bind_cols(
    st_coordinates(D_geo2$geometry) %>% as.tibble %>% select(level1_latitude=X,level1_longitude=Y),
    st_coordinates(D_geo2$level2) %>% as.tibble %>% select(level2_latitude=X,level2_longitude=Y),
    st_coordinates(D_geo2$level3) %>% as.tibble %>% select(level3_latitude=X,level3_longitude=Y),
    st_coordinates(D_geo2$level4) %>% as.tibble %>% select(level4_latitude=X,level4_longitude=Y)
  )
D_geo3 <- D_geo2 %>% as.data.frame %>% select(dataset:var_3) %>% 
  bind_cols(coords)



# Locate Matches 

# Truth 
tmatch = sim.data$true.matches 
truth = c(); for(i in 1:4) truth = c(truth,tmatch[,i]);truth = truth[truth!="0-0"]

# Locate matches
loc_match = function(output){
  B = bind_rows(output$processed$event_matched,
                output$processed$episode_matched) 
  pos = B %>% select(contains("data")) %>% colnames %>% str_replace("data","") 
  bring_together= function(x) unite(B,!!x,contains(!!x),sep="-") %>% select(!!x)
  set = pos %>% map(bring_together) %>% bind_cols
  set_vec = c(); for(i in 1:ncol(set)) set_vec = c(set_vec,set[,i])
  set_vec = set_vec[!set_vec %in% c("0-0","NA-NA")]
  set_vec
}


# Recovery rate
rates = function(found){
  tibble(tpr = sum(truth %in% found)/length(truth),
         fpr = sum(!found %in% truth)/ (1000-length(truth))
  )
}
popts = c("level1","level2","level3","level4")
distort = seq(.05,1,.05)   
windows = c(5,10,25,50,100,150,200,300,400,500)
n_sims = 30

for(s in 1:n_sims){
  if(s==1){
    entry10 = c(); pb = progress_estimated(n_sims)
  }
  for(w in windows){
    if(w==5){
      tmp_out = c()
    }
    for (i in 1:length(distort)){
      if(i == 1){
        dat = c()
      }
      # Assign Distortions
      tmp <- D_geo3
      N = nrow(D_geo3)
      ind = sample(1:N,N*distort[i],replace = T)
      tmp$reassign = 1
      tmp$latitude = NA;tmp$longitude = NA
      
      # randomly assign the imprecision unit and set coordinates to that unit.
      tmp[ind,"reassign"] = round(runif(N*distort[i],2,4)) 
      for(j in 1:length(popts)){tmp[tmp$reassign==j,c("latitude","longitude")] = tmp[tmp$reassign==j,grep(popts[j],colnames(tmp),value = T)]}
      
      tmp2 = tmp %>% 
        select(dataset:var_3,latitude,longitude,reassign) %>% 
        as.data.frame %>% 
        group_by(dataset) %>% nest
      
      # Assign
      for(k in 1:nrow(tmp2)){assign(tmp2$dataset[[k]],as.data.frame(tmp2$data[[k]]))}
      
      # Integrate
      output = meltt(D1,D2,D3,D4,taxonomies=tax,spatwindow=w,twindow=2,silent=TRUE)
      
      # Locate the matching entries -- save iteration
      dat = loc_match(output) %>% rates(.) %>% 
        mutate(distortion=paste0(distort[i]*100,"%")) %>% 
        bind_rows(dat,.)
    }
    # Save info
    tmp_out = dat %>% mutate(window = w) %>% bind_rows(tmp_out,.)
  }
  # Save info
  entry10 = tmp_out %>% mutate(simID = s) %>% bind_rows(entry10,.)
  pb$tick()$print() # track progress
}


# Need values for 0 distortion:
for(i in 1:nrow(D)){assign(D$dataset[[i]],as.data.frame(D$data[[i]]))} # Unaltered data
for(w in windows){
  if(w==5){ tmp_out = c()}
  output = meltt(D1,D2,D3,D4,taxonomies=tax,spatwindow=w,twindow=2,silent=TRUE)
  tmp_out = loc_match(output) %>% rates %>% 
    mutate(distortion="0%",window=w,simID=1) %>% 
    bind_rows(tmp_out,.)
}
entry10 = bind_rows(tmp_out,entry10) 

# Return TPR/FPR summary
entry10_summary = entry10 %>%
  group_by(window,distortion) %>%
  summarize(tpr=round(mean(tpr),3),
            fpr=round(mean(fpr),3))  %>%
  gather(key,val,-window,-distortion)
entry10_summary = cbind(entry10_summary[entry10_summary$key=="tpr",],entry10_summary[entry10_summary$key=="fpr",])[,c(1:2,4,8)]
names(entry10_summary) = c('window','distortion','TPR','FPR')
entry10_summary = arrange(entry10_summary,window,as.numeric(str_remove(distortion,"%"))/100)



# Simulations Correspond with Entry 11 in Table 1 ------------------------

# The level of geo-spatial precision as a proportion of the total data.
# Distortions both w/r/t levels and datasets are randomly assigned, so multiple
# simulations are run and the average is then reported. Finally, different
# spatial windows are explored to assess the rememdy. This simulation provides a
# variant to entry 10 by examining the effect of a weak taxonomy on this process.

# Simulate Data 

# Secondary Criteria (weak)
sc = data.frame(k=c(rep(1,2)),
                sc=c(2,1))
# Simulate Data
set.seed(123)
sim.data = sim.melt(N=1000,known.events = 45,known.episodes = 5,
                    known.episodes.events = 0,
                    n.datasets = 4,t.distort = 1,d.distort = 2,
                    episode.window = 5,
                    min_date = as.Date("2011-01-01"),max_date =as.Date("2012-01-01"),
                    lat.bounds = c(4.5,14),lon.bounds =c(3,14),
                    second.crit = sc)

D = sim.data$user.data %>% group_by(dataset) %>% nest() 
tax = sim.data$user.taxonomies


# Simulate Geographies
D_geo = D %>% unnest %>% st_as_sf(coords=c("latitude","longitude")) 

level4 = st_make_grid(D_geo,cellsize = 12) %>% 
  as.tibble %>% 
  mutate(level4 = st_centroid(geometry)) 

level3 = st_make_grid(D_geo,cellsize = 8) %>% as.tibble %>% 
  mutate(level3 = st_centroid(geometry)) 

level2 = st_make_grid(D_geo,cellsize = 2) %>% as.tibble %>% 
  mutate(level2 = st_centroid(geometry)) 

# Gather
x1 = st_intersection(D_geo,st_as_sf(level4)) %>% st_set_geometry(NULL)  %>% select(dataset,event,contains("level")) 
x2 = st_intersection(D_geo,st_as_sf(level3)) %>% st_set_geometry(NULL)  %>% select(dataset,event,contains("level")) 
x3 = st_intersection(D_geo,st_as_sf(level2)) %>% st_set_geometry(NULL)  %>% select(dataset,event,contains("level")) 

# Build coordinate system options upfront
D_geo2 = left_join(D_geo,x1) %>% left_join(x2) %>% left_join(x3)
coords <-
  bind_cols(
    st_coordinates(D_geo2$geometry) %>% as.tibble %>% select(level1_latitude=X,level1_longitude=Y),
    st_coordinates(D_geo2$level2) %>% as.tibble %>% select(level2_latitude=X,level2_longitude=Y),
    st_coordinates(D_geo2$level3) %>% as.tibble %>% select(level3_latitude=X,level3_longitude=Y),
    st_coordinates(D_geo2$level4) %>% as.tibble %>% select(level4_latitude=X,level4_longitude=Y)
  )
D_geo3 <- D_geo2 %>% as.data.frame %>% select(dataset:var_1) %>% 
  bind_cols(coords)



# Locate Matches 

# Truth 
tmatch = sim.data$true.matches 
truth = c(); for(i in 1:4) truth = c(truth,tmatch[,i]);truth = truth[truth!="0-0"]

# Locate matches
loc_match = function(output){
  B = bind_rows(output$processed$event_matched,
                output$processed$episode_matched) 
  pos = B %>% select(contains("data")) %>% colnames %>% str_replace("data","") 
  bring_together= function(x) unite(B,!!x,contains(!!x),sep="-") %>% select(!!x)
  set = pos %>% map(bring_together) %>% bind_cols
  set_vec = c(); for(i in 1:ncol(set)) set_vec = c(set_vec,set[,i])
  set_vec = set_vec[!set_vec %in% c("0-0","NA-NA")]
  set_vec
}


# Recovery rate
rates = function(found){
  tibble(tpr = sum(truth %in% found)/length(truth),
         fpr = sum(!found %in% truth)/ (1000-length(truth))
  )
}
popts = c("level1","level2","level3","level4")
distort = seq(.05,1,.05)   
windows = c(5,10,25,50,100,150,200,300,400,500)
n_sims = 30

for(s in 1:n_sims){
  if(s==1){
    entry11 = c(); pb = progress_estimated(n_sims)
  }
  for(w in windows){
    if(w==5){
      tmp_out = c()
    }
    for (i in 1:length(distort)){
      if(i == 1){
        dat = c()
      }
      # Assign Distortions
      tmp <- D_geo3
      N = nrow(D_geo3)
      ind = sample(1:N,N*distort[i],replace = T)
      tmp$reassign = 1
      tmp$latitude = NA;tmp$longitude = NA
      
      # randomly assign the imprecision unit and set coordinates to that unit.
      tmp[ind,"reassign"] = round(runif(N*distort[i],2,4)) 
      for(j in 1:length(popts)){tmp[tmp$reassign==j,c("latitude","longitude")] = tmp[tmp$reassign==j,grep(popts[j],colnames(tmp),value = T)]}
      
      tmp2 = tmp %>% 
        select(dataset:var_1,latitude,longitude,reassign) %>% 
        as.data.frame %>% 
        group_by(dataset) %>% nest
      
      # Assign
      for(k in 1:nrow(tmp2)){assign(tmp2$dataset[[k]],as.data.frame(tmp2$data[[k]]))}
      
      # Integrate
      output = meltt(D1,D2,D3,D4,taxonomies=tax,spatwindow=w,twindow=2,silent=TRUE)
      
      # Locate the matching entries -- save iteration
      dat = loc_match(output) %>% rates(.) %>% 
        mutate(distortion=paste0(distort[i]*100,"%")) %>% 
        bind_rows(dat,.)
    }
    # Save info
    tmp_out = dat %>% mutate(window = w) %>% bind_rows(tmp_out,.)
  }
  # Save info
  entry11 = tmp_out %>% mutate(simID = s) %>% bind_rows(entry11,.)
  pb$tick()$print() # track progress
}

# Need values for 0 distortion:
for(i in 1:nrow(D)){assign(D$dataset[[i]],as.data.frame(D$data[[i]]))} # Unaltered data
for(w in windows){
  if(w==5){ tmp_out = c()}
  output = meltt(D1,D2,D3,D4,taxonomies=tax,spatwindow = w,twindow=2,silent = TRUE)
  tmp_out = loc_match(output) %>% rates %>% 
    mutate(distortion="0%",window=w,simID=1) %>% 
    bind_rows(tmp_out,.)
}
entry11 = bind_rows(tmp_out,entry11) 

# Return TPR/FPR summary
entry11_summary = entry11 %>%
  group_by(window,distortion) %>%
  summarize(tpr=round(mean(tpr),3),
            fpr=round(mean(fpr),3))   %>%
  gather(key,val,-window,-distortion)
entry11_summary = cbind(entry11_summary[entry11_summary$key=="tpr",],entry11_summary[entry11_summary$key=="fpr",])[,c(1:2,4,8)]
names(entry11_summary) = c('window','distortion','TPR','FPR')
entry11_summary = arrange(entry11_summary,window,as.numeric(str_remove(distortion,"%"))/100)



