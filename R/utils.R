#save(myediteddata, file="data.rda")
#LazyData: true

#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
for (lib in c('dplyr', 'readr', 'downloader', 'logging')) {
  if (!suppressPackageStartupMessages(require(lib, character.only = TRUE)) ) stop(paste("Please install the R package: ",lib))
}
set_coloures <- function(){
  site_colours <- list()
  site_colours$Reference_genome <- "white"
  site_colours$all_features <- rgb(0, 38, 84, maxColorValue = 255) #"gray60"
  site_colours$pathway <- rgb(226, 203, 146, maxColorValue = 255) #"orange"
  return(site_colours)
}
#-------------------------------------------------------------------------------------------------
set_names <- function(){
  site_names <- list()
  # Legend entries will appear in the order given here
  site_names$all_features <- "All omic fesatures in the community"
  site_names$pathway <- "Omic features contributed in the pathway"
  return(site_names)
}

load_pathway2HMDBID <- function(path, pthahway_col_number = 3, hmdb_col_number = 7){
  set_elements_mapper <- scan(gzfile(path), character(), what="", sep="\n") # header=F, fill=TRUE,
  y <- strsplit(set_elements_mapper, "\t")

  #remove the first line
  y <- y[-1]

  mapper <- list()
  for (value in y){
    #print(mapper[value[2]])
    if (is.null(mapper[[value[3]]]))
      mapper[value[3]] <- value[7]
    else {
      temp_list <- append(mapper[[value[3]]], value[7])
      mapper[value[3]] <- NULL
      mapper[value[3]] <- list(temp_list)
    }
  }
  return (mapper)
}
downloadMappingFiles <- function(db_path){
  # Download from http://smpdb.ca/downloads
  logging::loginfo("Start downloading the SMPDB metabolites database!")
  url <- 'http://smpdb.ca/downloads/smpdb_metabolites.csv.zip'
  download(url, dest=paste(db_path,"/smpdb_metabolites.csv.zip", sep=''), mode="wb")

  # extract the zip file to smpdb_metabolites
  unzip (paste(db_path,"/smpdb_metabolites.csv.zip", sep=''), exdir = paste(db_path,"/smpdb_metabolites", sep =''))
  logging::loginfo("Downloading the SMPDB metabolites database is Done!")
}
merge_cvs <- function(input_dir_path, outputfile){
  # Download from http://smpdb.ca/downloads
  # extract the zip file to smpdb_metabolites
  #input_dir_path = "~/Downloads/smpdb_metabolites"
  #outputfile = '/Users/rah/Documents/Metabolomics/example_data/HMDB/merged_ref_db.csv'
  library(dplyr)
  library(readr)
  df <- list.files(path=input_dir_path, full.names = TRUE) %>%
    lapply(read_csv) %>%
    bind_rows
  write_csv( df, outputfile)
}
setup_smpdb_metabolites_db <- function(db_path, force = FALSE, rm_intermediate_files = TRUE ){

  #Check its existence
  if (file.exists(paste(db_path,"/smpdb_metabolites.csv", sep =''))){

    #Delete file if it exist
    if (force)
      file.remove(paste(db_path,"/smpdb_metabolites.csv", sep =''))
    else {
      print(sprintf("The %s is exist! If you wish to overwrite it p lease use force =TRUE in the function parameters.",
                      paste(db_path,"/smpdb_metabolites.csv", sep ='')))
      return(paste(db_path,"/smpdb_metabolites.csv", sep =''))
    }
  }
    # create an output folder if it does not exist
  if (!file.exists(db_path)) {
    print("Creating database path folder ...")
    dir.create(db_path)
  }
  # download files and put them to the db_path
  logging::logdebug("Downloading SMPDB metabolites database started and will located under %s.",db_path)
  downloadMappingFiles(db_path)
  logging::logdebug("Downloading process is done!")
  logging::logdebug("Merging all files ...")
  merge_cvs(paste(db_path,"/smpdb_metabolites", sep =''), paste(db_path,"/smpdb_metabolites.csv", sep =''))
  logging::logdebug("Merging all files is done!")
  # clean intermdiate data
  if (rm_intermediate_files){
    file.remove(paste(db_path,"/smpdb_metabolites.csv.zip", sep=''))
    unlink(paste(db_path,"/smpdb_metabolites", sep=''), recursive = T)
  }
  logging::logdebug("The SMPDB metabolites database i created under %s.",paste(db_path,"/smpdb_metabolites.csv", sep =''))
  return(paste(db_path,"/smpdb_metabolites.csv", sep =''))
}
