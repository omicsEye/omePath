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
if(identical(environment(), globalenv()) &&
   !length( grep( "^source\\(", sys.calls()))) {

  # source all R in deepath package, relative to this folder
  script_options <-commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=","",script_options[grep("--file=", script_options)])
  script_dir <- dirname(script_path)
  script_name <- basename(script_path)

  for(R_file in dir(script_dir, pattern = "*.R"))
  {
    if(! ( R_file == script_name )) source( file.path(script_dir, R_file) )
  }
}


option_not_valid_error <- function(message, valid_options) {
  logging::logerror(paste(message,": %s"), toString(valid_options))
  stop("Option not valid", call.=FALSE)
}

###########################
# Set the default options #
###########################

method_choices <- c("gset","ks", 'wilcox')

##########################################################################################
# Main maaslin2 function with defaults set to the same as those used on the command line #
##########################################################################################
#' @export
deepath <- function(stats_table,
                   output,
                   score_col = 'logFC',
                   pval_threshold = 0.05,
                   fdr_threshold = NA,
                   Pathway.Subject = 'Metabolic',
                   method = 'gset',
                   min_member=2,
                   mapper_file=NA,
                   do_plot = TRUE,
                   pathway_col = "Pathway",
                   feature_col = "Feature")
{
  #################################################################
  # Read in the scor set data and create output folder, init log  #
  #################################################################
  # if a character string then this is a file name, else it is a data frame
  if (is.character(stats_table)) {
    stats_table <- data.frame(data.table::fread(stats_table, header = TRUE, sep = "\t"), row.names = 1)
    if (nrow(stats_table) == 1) {
      # read again to get row name
      stats_table <- read.table(stats_table, header = TRUE, row.names = 1)
    }
  } else {
    stats_table <- stats_table
  }

  # create an output folder if it does not exist
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }

  # create log file (write info to stdout and debug level to log file)
  # set level to finest so all log levels are reviewed
  log_file <- file.path(output,"deepath.log")
  # remove log file if already exists (to avoid append)
  if (file.exists(log_file)) {
    print(paste("Warning: Deleting existing log file:", log_file))
    unlink(log_file)
  }
  logging::basicConfig(level = 'FINEST')
  logging::addHandler(logging::writeToFile,file = log_file,level = "DEBUG")
  logging::setLevel(20, logging::getHandler('basic.stdout'))

  #####################
  # Log the arguments #
  #####################

  logging::loginfo("Writing function arguments to log file")
  logging::logdebug("Function arguments")
  if (is.character(score_col)) {
    logging::logdebug("Input data file: %s", score_col)
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
  if (!method %in% method_choices) {
    option_not_valid_error("Please select a pvalue method from the list of available options", toString(method_choices))
  }

  ####################################
  # apply the method to the data with the correction
  ####################################

  site_colours <- set_coloures()
  site_names <- set_names()

  logging::loginfo("Running selected analysis method: %s", method)
  results <- OSEA(stats_table = stats_table, score_col = score_col,
                                     pval_threshold = .05,
                                     Pathway.Subject = 'Metabolic',
                                     output = output,
                                     do_plot = TRUE,
                                     mapper_file = mapper_file,
                                     method = method,
                                     min_member = min_member,
                                     pathway_col = pathway_col,
                                     feature_col = feature_col)

  #########################
  # Write out the results #
  #########################

  # write stats table result to file
  enrichment_stats_file = file.path(output, "enrichment_stats.tsv")
  # remove stats table file if already exists (since stats table append)
  if (file.exists(enrichment_stats_file)) {
    logging::logwarn("Deleting existing stats table file: %s", enrichment_stats_file)
    unlink(enrichment_stats_file)
  }
  logging::loginfo("Writing stats table to file %s", enrichment_stats_file)
  write.table(results$enrichment_stats, file = enrichment_stats_file, sep = "\t", quote = FALSE, row.names =
                FALSE)

  #######################################################
  # Create visualizations for results passing threshold #
  #######################################################

  if (do_plot) {
    logging::loginfo("Writing enrichment plots (one for each significant association) to output folder: %s", output)
  }

  return(results)
}
