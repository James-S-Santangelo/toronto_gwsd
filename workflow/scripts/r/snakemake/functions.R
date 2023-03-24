# Functions sourced by other scripts in Snakemake pipeline

#' Fit Shape Constrained cubic p-spline with 10 knots
#'
#' @param Dataframe with terms to model
#'
#' @return Fitted SCAM model
fit_scam <- function(df){
  fit <- scam(cM ~ s(pos, bs='mpi', k=20, m=2),
              family=gaussian(link="identity"),
              data = df)
  return(fit)
}

#' Interpolate predicted cM from SCAM model
#' 
#' @param fit Fitted SCAM model
#' @param sites Sites to be interpolation
interpolate <- function(fit, sites){
  newdata <- data.frame(pos = sites$pos)
  preds <- predict(fit, newdata)
  return(preds)
}

#' Fit SCAM model and add predicted cM to dataframe
#' 
#' @param chrom_name Name of chromosome to be fit
#' @param markers_df Dataframe with markers for model fitting
#' @param sites_df Dataframe with sites to be interpolated
#'
#' @return Dataframe with predicted cM at each position
get_fits <- function(chrom_name, markers_df, sites_df){
  
  print(sprintf("Fitting chromosome %s", chrom_name))
  markers <- markers_df %>% filter(chrom == chrom_name)
  sites <- sites_df %>% filter(chrom == chrom_name)
  fit <- fit_scam(df = markers)
  predicted <- interpolate(fit, sites)
  
  sites$preds <- predicted
  return(sites)
}

#' Plot SCAM fit
#'
#' @param chrom_name Name of chromosome to be fit
#' @param markers_df Dataframe with markers for model fitting
#' @param scam_fits_df Dataframe with predicted values from SCAM fit
#'
#' @return ggplot object
scam_plot <- function(chrom_name, markers_df, scam_fits_df){
  
  markers <- filter(markers_df, chrom == chrom_name)
  scamFit <- filter(scam_fits_df, chrom == chrom_name)
  
  plot <- markers %>%
    ggplot(., aes(x = pos, y = cM)) +
    geom_point(size = 2.5, color = 'black', shape = 21, fill = 'black', show.legend = FALSE) +
    geom_line(data = scamFit, aes(x = pos, y = preds), size = 1, color = 'red') +
    xlab('Physical position') + ylab('Genetic position (cM)') +
    ggtitle(chrom_name) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 13))
  
  return(plot)
}
