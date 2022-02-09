# Functions sourced by other scripts in Snakemake pipeline

#' Load results from marker BLAST
#' 
#' @param path Path to marker blast results
#' 
#' @return Tibble/Dataframe with BLAST results
load_marker_blast <- function(path){
  
    # Column names
    cols <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'slen', 'evalue', 'bitscore', 'qcovs', 'qcovhsp')

    # Read in BLAST results
    markers_blast <- read.delim(path, sep = '\t', col.names = cols) %>%
        mutate(qseqid = as.character(qseqid)) %>%
        separate(qseqid, into = c('prev_chrom', 'prev_pos'), sep = ':', remove = FALSE) %>%
        
        # Re-map chromosome names from Olsen (which are based on LG) to match correct chromosomes
        mutate(prev_chrom = case_when(prev_chrom == 1 ~ 'CM019101.1',
                                      prev_chrom == 2 ~ 'CM019103.1',
                                      prev_chrom == 3 ~ 'CM019105.1',
                                      prev_chrom == 4 ~ 'CM019107.1',
                                      prev_chrom == 5 ~ 'CM019109.1',
                                      prev_chrom == 6 ~ 'CM019111.1',
                                      prev_chrom == 7 ~ 'CM019113.1',
                                      prev_chrom == 8 ~ 'CM019115.1',
                                      prev_chrom == 9 ~ 'CM019102.1',
                                      prev_chrom == 10 ~ 'CM019104.1',
                                      prev_chrom == 11 ~ 'CM019106.1',
                                      prev_chrom == 12 ~ 'CM019108.1',
                                      prev_chrom == 13 ~ 'CM019110.1',
                                      prev_chrom == 14 ~ 'CM019112.1',
                                      prev_chrom == 15 ~ 'CM019114.1',
                                      prev_chrom == 16 ~ 'CM019116.1')) %>%
        mutate(prev_id = paste(prev_chrom, prev_pos, sep = ':')) %>%
        mutate(curr_chrom = str_extract(string = sseqid, pattern = '(?<=gb\\|)(.+)(?=\\|)'),
               curr_pos = sstart + 99,
               curr_id =  paste(curr_chrom, curr_pos, sep = ':'))
    return(markers_blast)
}

#' Function to load genetic map from Olsen et al. (2021)
#' 
#' @param path Path to genetic map
#' 
#' @return Genetic map as Dataframe
load_genMap <- function(path){
    
    # Load genetic map
    genMap <- read.csv(path) %>%
        
        # Re-map linkage group names to chromosomes
        # Mean that linkage group in marker name won't match chromosome number, but that's fine
        mutate(LG = case_when(LG == 1 ~ 'CM019101.1',
                              LG == 2 ~ 'CM019103.1',
                              LG == 3 ~ 'CM019105.1',
                              LG == 4 ~ 'CM019107.1',
                              LG == 5 ~ 'CM019109.1',
                              LG == 6 ~ 'CM019111.1',
                              LG == 7 ~ 'CM019113.1',
                              LG == 8 ~ 'CM019115.1',
                              LG == 9 ~ 'CM019102.1',
                              LG == 10 ~ 'CM019104.1',
                              LG == 11 ~ 'CM019106.1',
                              LG == 12 ~ 'CM019108.1',
                              LG == 13 ~ 'CM019110.1',
                              LG == 14 ~ 'CM019112.1',
                              LG == 15 ~ 'CM019114.1',
                              LG == 16 ~ 'CM019116.1'))
    return(genMap)
}

#' Plot markers for both populations (colors) with outliers (X)
#'
#' @param chrom_name Character vector with name of chromosome to plot
#' @param df Dataframe containing X and Y variables
#'
#' @return ggplot object
plot_markers_byPop <- function(chrom_name, df){

    df_sub <- df %>% filter(chrom == chrom_name)
    cols <- c('#00A08A', '#F98400')
    plot <- df_sub %>%
        filter(exclude != 1) %>%
        ggplot(., aes(x = pos, y = cM)) +
        geom_point(size = 2, shape = 21, aes(fill = pop), show.legend = FALSE) +
        geom_point(data = df_sub %>% filter(exclude == 1),
                   shape = 4, show.legend = FALSE, size = 2) +
        xlab('Physical position') + ylab('Genetic position (cM)') +
        ggtitle(chrom_name) +
        scale_fill_manual(values = cols) +
        theme_classic() +
        theme(axis.title = element_text(size = 15),
              axis.text = element_text(size = 13))
    return(plot)
}

#' Filter dataframe to only include markers in window
#'
#' @param window_start Start position of window
#' @param window_size Size of window in bp
#' @param step Number of base pairs to shift window
#' @param df Dataframe with markers 
#'
#' @return Filtered dataframe containing only markers in windows'
get_windows <- function(window_start, window_size, step, df){

    df_filt <- df %>%
        arrange(pos) %>%
        filter(pos >= window_start & pos < window_start + window_size) %>%

        # Get min, max, and mid position of window
        mutate(min_pos = window_start,
               max_pos = window_start + (window_size - step) - 1,
               mid_pos = window_start + (max_pos - min_pos) / 2)
    return(df_filt)
}

#' Calculate mean cM of markers in window
#' 
#' @param chrom_name Character vector with name of chromosome
#' @param window_size Size of window in bp
#' @param step Number of base pairs to shift window
#' @param df Dataframe with windowed markers
#'
#' @return Dataframe with mean cM in windows
calculate_windowed_means <- function(chrom_name, window_size, step, df){

    df_sub <- df %>% filter(chrom == chrom_name)
    all_pos <- seq(from = min(df_sub$pos), to = max(df_sub$pos), by = step)
    df_withWindows <- purrr::map_dfr(all_pos, get_windows,
                                   df = df_sub,
                                   window_size = window_size,
                                   step = step,
                                   .id = 'win_id')
    df_winMeans <- df_withWindows %>%
        group_by(win_id, min_pos, max_pos, mid_pos, pop) %>%
        summarise(num_markers = n(),
                  win_mean_cM = mean(cM), .groups = 'keep') %>%
        ungroup() %>%
        group_by(win_id, min_pos, max_pos, mid_pos) %>%
        summarise(win_mean_cM = mean(win_mean_cM), .groups = 'keep') %>%
        mutate(chrom = chrom_name)
    return(df_winMeans)
}

#' Fit Shape Constrained cubic p-spline with 10 knots
#'
#' @param Dataframe with terms to model
#'
#' @return Fitted SCAM model
fit_scam <- function(df){
    fit <- scam(win_mean_cM~s(mid_pos, bs='mpi', k=10, m=2),
              family=gaussian(link="identity"),
              data = df)
    return(fit)
}

#' Interpolate predicted cM from SCAM model
#' 
#' @param fit Fitted SCAM model
#' @param sites Sites to be interpolation
interpolate <- function(fit, sites){
    newdata <- data.frame(mid_pos = sites$pos)
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
        ggplot(., aes(x = mid_pos, y = win_mean_cM)) +
        geom_point(size = 2.5, color = 'black', shape = 21, fill = 'black', show.legend = FALSE) +
        geom_line(data = scamFit, aes(x = pos, y = preds), size = 1, color = 'red') +
        xlab('Physical position') + ylab('Genetic position (cM)') +
        ggtitle(chrom_name) +
        theme_classic() +
        theme(axis.title = element_text(size = 15),
              axis.text = element_text(size = 13))

    return(plot)
}


