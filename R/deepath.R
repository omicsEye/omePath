#!/usr/bin/env Rscript

###############################################################################
# deepath

# Copyright (c) 2019 the Broad Institute of MIT and Harvard

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
###############################################################################

# load in the required libraries, report an error if they are not installed
#for( lib in c('gsEasy','future', 'limma', 'ggplot2')) {
#  if(! suppressPackageStartupMessages(require(lib, character.only=TRUE)) ) stop(paste("Please install the R package: ",lib))
#}

###############################################################
# If running on the command line, load other deepath modules #
###############################################################

# this evaluates to true if script is being called directly as an executable
if (identical(environment(), globalenv()) &&
    !length(grep("^source\\(", sys.calls()))) {
  # source all R in deepath package, relative to this folder
  script_options <- commandArgs(trailingOnly = FALSE)
  script_path <-
    sub("--file=", "", script_options[grep("--file=", script_options)])
  script_dir <- dirname(script_path)
  script_name <- basename(script_path)
  
  for (R_file in dir(script_dir, pattern = "*.R"))
  {
    if (!(R_file == script_name))
      source(file.path(script_dir, R_file))
  }
}


option_not_valid_error <- function(message, valid_options) {
  logging::logerror(paste(message, ": %s"), toString(valid_options))
  stop("Option not valid", call. = FALSE)
}

###########################
# Set the default options #
###########################

method_choices <- c("gset", "ks", 'wilcox')

##########################################################################################
# Main deepath function with defaults set to the same as those used on the command line #
##########################################################################################
#' @export
deepath <- function(input_data,
                    output,
                    input_metadata = NA,
                    meta = NA,
                    case_label = NA,
                    control_label = NA,
                    score_col = 'logFC',
                    pval_threshold = 0.05,
                    fdr_threshold = NA,
                    Pathway.Subject = NA,
                    method = 'ks',
                    min_member = 2,
                    mapper_file = NA,
                    do_plot = TRUE,
                    pathway_col = "Pathway",
                    feature_col = "Feature")
{
  #################################################################
  # Read in the data and metadata, create output folder, init log #
  #################################################################
  # is the metadata not provides then input data should be a score file
  if (is.character(input_data)) {
    data <-
      data.frame(
        data.table::fread(
          input_data,
          header = TRUE,
          check.names = FALSE,
          sep = "\t"
        ),
        row.names = 1
      )
    
    if (nrow(data) == 1) {
      # read again to get row name
      data <-
        read.table(
          input_data,
          header = TRUE,
          check.names = FALSE,
          row.names = 1
        )
    }
  } else {
    data <- input_data
  }
  #if meatdat is not provides the data should have column with score
  if (is.na(input_metadata)) {
    if (!score_col %in% colnames(data)) {
      print ("Please provide metadata or your score file shouls have the score column!!!")
    }
  }
  
  if (is.character(input_metadata)) {
    metadata <-
      data.frame(
        data.table::fread(
          input_metadata,
          header = TRUE,
          check.names = FALSE,
          sep = "\t"
        ),
        row.names = 1
      )
    if (nrow(metadata) == 1) {
      metadata <- read.table(
        input_metadata,
        header = TRUE,
        check.names = FALSE,
        row.names = 1
      )
    }
  } else {
    metadata <- input_metadata
  }
  ###############################################################
  # Determine orientation of data in input and reorder to match #
  ###############################################################
  
  logging::loginfo("Determining format of input files")
  samples_row_row <- intersect(rownames(data), rownames(metadata))
  if (length(samples_row_row) > 0) {
    # this is the expected formatting so do not modify data frames
    logging::loginfo(paste(
      "Input format is data samples",
      "as rows and metadata samples as rows"
    ))
  } else {
    samples_column_row <- intersect(colnames(data), rownames(metadata))
    if (length(samples_column_row) > 0) {
      logging::loginfo(paste(
        "Input format is data samples",
        "as columns and metadata samples as rows"
      ))
      # transpose data frame so samples are rows
      data <- as.data.frame(t(data))
      logging::logdebug("Transformed data so samples are rows")
    } else {
      samples_column_column <-
        intersect(colnames(data), colnames(metadata))
      if (length(samples_column_column) > 0) {
        logging::loginfo(
          paste(
            "Input format is data samples",
            "as columns and metadata samples as columns"
          )
        )
        data <- as.data.frame(t(data))
        metadata <- as.data.frame(t(metadata))
        logging::logdebug("Transformed data and metadata so samples are rows")
      } else {
        samples_row_column <-
          intersect(rownames(data), colnames(metadata))
        if (length(samples_row_column) > 0) {
          logging::loginfo(
            paste(
              "Input format is data samples",
              "as rows and metadata samples as columns"
            )
          )
          metadata <- as.data.frame(t(metadata))
          logging::logdebug("Transformed metadata so samples are rows")
        } else {
          logging::logerror(
            paste(
              "Unable to find samples in data and",
              "metadata files.",
              "Rows/columns do not match."
            )
          )
          logging::logdebug("Data rows: %s",
                            paste(rownames(data), collapse = ","))
          logging::logdebug("Data columns: %s",
                            paste(colnames(data), collapse = ","))
          logging::logdebug("Metadata rows: %s",
                            paste(rownames(metadata), collapse = ","))
          logging::logdebug("Metadata columns: %s",
                            paste(colnames(data), collapse = ","))
          stop()
        }
      }
    }
  }
  
  # replace unexpected characters in feature names
  #colnames(data) <- make.names(colnames(data))
  
  # check for samples without metadata
  extra_feature_samples <-
    setdiff(rownames(data), rownames(metadata))
  if (length(extra_feature_samples) > 0)
    logging::logdebug(
      paste(
        "The following samples were found",
        "to have features but no metadata.",
        "They will be removed. %s"
      ),
      paste(extra_feature_samples, collapse = ",")
    )
  
  # check for metadata samples without features
  extra_metadata_samples <-
    setdiff(rownames(metadata), rownames(data))
  if (length(extra_metadata_samples) > 0)
    logging::logdebug(
      paste(
        "The following samples were found",
        "to have metadata but no features.",
        "They will be removed. %s"
      ),
      paste(extra_metadata_samples, collapse = ",")
    )
  
  # get a set of the samples with both metadata and features
  intersect_samples <-
    intersect(rownames(data), rownames(metadata))
  logging::logdebug(
    "A total of %s samples were found in both the data and metadata",
    length(intersect_samples)
  )
  
  # now order both data and metadata with the same sample ordering
  logging::logdebug("Reordering data/metadata to use same sample ordering")
  data <- data[intersect_samples, , drop = FALSE]
  metadata <- metadata[intersect_samples, , drop = FALSE]
  
  
  # create an output folder and figures folder if it does not exist
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }
  
  figures_folder <- file.path(output, "figures")
  if (!file.exists(figures_folder)) {
    print("Creating output figures folder")
    dir.create(figures_folder)
  }
  
  # create log file (write info to stdout and debug level to log file)
  # set level to finest so all log levels are reviewed
  log_file <- file.path(output, "deepath.log")
  # remove log file if already exists (to avoid append)
  if (file.exists(log_file)) {
    print(paste("Warning: Deleting existing log file:", log_file))
    unlink(log_file)
  }
  logging::basicConfig(level = 'FINEST')
  logging::addHandler(logging::writeToFile,
                      file = log_file, level = "DEBUG")
  logging::setLevel(20, logging::getHandler('basic.stdout'))
  
  #####################
  # Log the arguments #
  #####################
  
  logging::loginfo("Writing function arguments to log file")
  logging::logdebug("Function arguments")
  if (is.character(input_data)) {
    logging::logdebug("Input data file: %s", input_data)
  }
  if (is.character(input_metadata)) {
    logging::logdebug("Input metadata file: %s", input_metadata)
  }
  logging::logdebug("Output folder: %s", output)
  logging::logdebug("Score type: %s", score_col)
  logging::logdebug("p-value threshold: %s", pval_threshold)
  logging::logdebug("Pathway subject: %s", Pathway.Subject)
  
  ####################################
  # Check valid options are selected #
  ####################################
  
  # Check valid normalization option selected
  logging::loginfo("Verifying options selected are valid")
  
  # check valid pvalue calculation method selected
  #if (!method %in% method_choices) {
  #  option_not_valid_error("Please select a pvalue method from the list of available options", toString(method_choices))
  #}
  
  
  ####################################
  # apply the method to the data with the correction
  ####################################
  
  #site_colours <- set_coloures()
  #site_names <- set_names()
  if (!is.na(input_metadata)) {
    stats_table <-
      test2groups(
        data ,
        metadata,
        meta,
        case_label,
        control_label,
        test_type = 'wilcox.test',
        paired = F
      )
  } else{
    stats_table <- data
  }
  # write stats table result to file
  stats_file = file.path(output, "stats_table.tsv")
  # remove stats table file if already exists (since stats table append)
  if (file.exists(stats_file)) {
    logging::logwarn("Deleting existing stats table file: %s",
                     stats_file)
    unlink(stats_file)
  }
  logging::loginfo("Writing stats table to file %s", stats_file)
  write.table(
    stats_table,
    file = stats_file,
    sep = "\t",
    quote = FALSE,
     row.names =
      FALSE
   )
  
  logging::loginfo("Running selected analysis method: %s", method)
  results <- OSEA(
    stats_table = stats_table,
    score_col = score_col,
    pval_threshold = pval_threshold,
    Pathway.Subject = Pathway.Subject,
    output = output,
    do_plot = do_plot,
    mapper_file = mapper_file,
    method = method,
    min_member = min_member,
    pathway_col = pathway_col,
    feature_col = feature_col
  )
  
  #########################
  # Write out the results #
  #########################
  
  # write stats table result to file
  enrichment_stats_file = file.path(output, "enrichment_stats.tsv")
  # remove stats table file if already exists (since stats table append)
  if (file.exists(enrichment_stats_file)) {
    logging::logwarn("Deleting existing stats table file: %s",
                     enrichment_stats_file)
    unlink(enrichment_stats_file)
  }
  logging::loginfo("Writing enrichment stats table to file %s", enrichment_stats_file)
  #write.table(
  #  results$enrichment_stats,
  #  file = enrichment_stats_file,
  #  sep = "\t",
  #  quote = FALSE,
 #   row.names =
  #    FALSE
 # )
  
  #######################################################
  # Create visualizations for results passing threshold #
  #######################################################
  
  if (do_plot) {
    logging::loginfo(
      "Writing enrichment plots (one for each significant association) to output folder: %s",
      output
    )
  }
  
  return(results)
}
