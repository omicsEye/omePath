#save(myediteddata, file="data.rda")
#LazyData: true

#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
for (lib in c('dplyr', 'readr', 'downloader', 'logging')) {
  if (!suppressPackageStartupMessages(require(lib, character.only = TRUE)))
    stop(paste("Please install the R package: ", lib))
}
set_coloures <- function() {
  site_colours <- list()
  #site_colours$Reference_genome <- "white"
  site_colours$all_features <-
    grDevices::rgb(0, 38, 84, maxColorValue = 255) #"gray60"
  site_colours$pathway <-
    'darkgoldenrod' #rgb(226, 203, 146, maxColorValue = 255) #"orange"
  return(site_colours)
  
}
#-------------------------------------------------------------------------------------------------
set_names <- function() {
  site_names <- list()
  # Legend entries will appear in the order given here
  site_names$all_features <- "Members not in the pathway"
  site_names$pathway <- "Members in the pathway"
  return(site_names)
}

load_pathway2HMDBID <-
  function(path,
           pthahway_col_number = 3,
           hmdb_col_number = 7) {
    set_elements_mapper <-
      scan(gzfile(path), character(), what = "", sep = "\n") # header=F, fill=TRUE,
    y <- strsplit(set_elements_mapper, "\t")
    
    #remove the first line
    y <- y[-1]
    
    mapper <- list()
    for (value in y) {
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
downloadMappingFiles <- function(db_path) {
  # Download from http://smpdb.ca/downloads
  logging::loginfo("Start downloading the SMPDB metabolites database!")
  url <- 'http://smpdb.ca/downloads/smpdb_metabolites.csv.zip'
  downloader::download(
    url,
    dest = paste(db_path, "/smpdb_metabolites.csv.zip", sep = ''),
    mode = "wb"
  )
  
  # extract the zip file to smpdb_metabolites
  utils::unzip(
    paste(db_path, "/smpdb_metabolites.csv.zip", sep = ''),
    exdir = paste(db_path, "/smpdb_metabolites", sep = '')
  )
  logging::loginfo("Downloading the SMPDB metabolites database is Done!")
}
merge_cvs <- function(input_dir_path, outputfile) {
  # Download from http://smpdb.ca/downloads
  # extract the zip file to smpdb_metabolites
  #input_dir_path = "~/Downloads/smpdb_metabolites"
  #outputfile = '/Users/rah/Documents/Metabolomics/example_data/HMDB/merged_ref_db.csv'
  
  df <- list.files(path = input_dir_path, full.names = TRUE) %>%
    lapply(readr::read_csv) %>%
    dplyr::bind_rows()
  readr::write_csv(df, outputfile)
}
setup_smpdb_metabolites_db <-
  function(db_path,
           force = FALSE,
           rm_intermediate_files = TRUE) {
    #Check its existence
    if (file.exists(paste(db_path, "/smpdb_metabolites.csv", sep = ''))) {
      #Delete file if it exist
      if (force)
        file.remove(paste(db_path, "/smpdb_metabolites.csv", sep = ''))
      else {
        print(
          sprintf(
            "The %s is exist! If you wish to overwrite it p lease use force =TRUE in the function parameters.",
            paste(db_path, "/smpdb_metabolites.csv", sep = '')
          )
        )
        return(paste(db_path, "/smpdb_metabolites.csv", sep = ''))
      }
    }
    # create an output folder if it does not exist
    if (!file.exists(db_path)) {
      print("Creating database path folder ...")
      dir.create(db_path)
    }
    # download files and put them to the db_path
    logging::logdebug("Downloading SMPDB metabolites database started and will located under %s.",
                      db_path)
    downloadMappingFiles(db_path)
    logging::logdebug("Downloading process is done!")
    logging::logdebug("Merging all files ...")
    merge_cvs(
      paste(db_path, "/smpdb_metabolites", sep = ''),
      paste(db_path, "/smpdb_metabolites.csv", sep = '')
    )
    logging::logdebug("Merging all files is done!")
    # clean intermdiate data
    if (rm_intermediate_files) {
      file.remove(paste(db_path, "/smpdb_metabolites.csv.zip", sep = ''))
      unlink(paste(db_path, "/smpdb_metabolites", sep = ''), recursive = T)
    }
    logging::logdebug(
      "The SMPDB metabolites database i created under %s.",
      paste(db_path, "/smpdb_metabolites.csv", sep = '')
    )
    return(paste(db_path, "/smpdb_metabolites.csv", sep = ''))
  }

#' Implementation of the GPD-based p-value estimation algorithm
#' @name GPD Permutation Test
#'
#' @param x0 an input vector of values
#' @param y  an input vector of values
#' @param minM_epdf a minimumn number of iterations for EPDF
#' @param Nexc a number of permution to performe in each step default 250
#' @param Nexc_shrink default 10
#' @param Nexc_alpha default 0.05
#' @param yfun default
#' @param ystart default 200
#' @param ygrow default 100
#' @param ymax default 5000
#' @return a `pvalue`
#' @examples
#' gpd_permutation_test(x, y)
#' @references
#' Implementation of the GPD-based p-value estimation algorithm from
#' Knijnenburg, Wessels, Reinders, and Shmulevich (2009) Fewer
#' permutations, more accurate P-values. Bioinformatics 25(12): i161–i168.
#' DOI: 10.1093/bioinformatics/btp211
#' Implemented by Jason Lloyd-Price and Ali Rahnavard

#' @export
gpd_permutation_test <- function(x0,
                                 y,
                                 minM_epdf = 10,
                                 Nexc = 250,
                                 Nexc_shrink = 10,
                                 Nexc_alpha = 0.05,
                                 yfun,
                                 ystart = 200,
                                 ygrow = 100,
                                 ymax = 5000) {
  if (missing(y)) {
    # Automatic sampling of y
    stopifnot(!missing(yfun))
    y <- sapply(rep(NA, ystart), function(i)
      yfun())
    
    while (sum(y > x0) < minM_epdf && length(y) < ymax) {
      y <- c(y, sapply(rep(NA, ygrow), function(i)
        yfun()))
    }
  }
  
  # Use the ecdf when there are at least 10 (default) exceedances
  M <- sum(y > x0)
  if (M >= minM_epdf) {
    return (M / length(y))
  }
  
  # Not enough exceedances to be confident about the p-value.. use the GPD
  while (T) {
    # Fit a GPD
    t <-
      mean(y[max(ceiling(length(y) / 2), length(y) - Nexc) + c(0, 1)])
    fit <- extRemes::fevd(y, threshold = t, type = "GP")
    
    # Test goodness-of-fit
    gof <-
      stats::ks.test(extRemes::pextRemes(fit, y[y > t]), stats::punif)
    if (Nexc <= 50 || gof$p.value >= Nexc_alpha) {
      if (Nexc <= 50) {
        warning(
          "GPD fit to tail samples fits poorly even with (<= 50) null samples. Consider adding more null samples."
        )
      }
      # Didn't reject.. fit is good enough
      break
    } else {
      Nexc <- Nexc - Nexc_shrink
    }
  }
  
  # Get the p-value estimate from the fit
  p <- mean(y > t) * extRemes::pextRemes(fit, x0, lower.tail = F)
  
  return (p)
}
