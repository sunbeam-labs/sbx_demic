spread_to_numeric_matrix <- function (data, row_key, col_key, value) {
  data <- dplyr::select_(data, row_key, col_key, value)
  data_wide <- tidyr::spread(data, col_key, value, fill=0)
  data_wide <- tibble::column_to_rownames(data_wide, row_key)
  as.matrix(as.data.frame(data_wide))
}
#linear mixed-effect model functions
# tidy lmer for the lmTest object
tidy_lmer <- function(lmer_test) {
  mod <- summary(lmer_test)
  data.frame(term  = rownames(mod$tTable), mod$tTable, row.names=NULL)
}

# count based lmer (for taxa abundances) with random effects
run_lmer <- function(cts_toTest, s_toTest, form1, rep_mes_label, p_cutoff) {
  rep_mes_form <- paste("~ 1 |", rep_mes_label)
  cts_toTest[,s_toTest$SampleID] %>%
    melt() %>%
    mutate(value = value+1) %>%
    setNames(c("Taxa", "SampleID", "Abundance")) %>%
    merge(s_toTest, by="SampleID") %>%
    mutate(props = Abundance / read_counts) %>%
    mutate(props100 = props * 100) %>%
    mutate(props_logit = log(props/(1-props))) %>%
    group_by(Taxa) %>%
    do(tidy_lmer(nlme::lme(as.formula(form1), random = as.formula(rep_mes_form), data=., na.action=na.omit))) %>%
    ungroup() %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value<p_cutoff)
}

#just get the lmer
get_lmer <- function(cts_toTest, s_toTest, form1, rep_mes_label, p_cutoff) {
  rep_mes_form <- paste("~ 1 |", rep_mes_label)
  cts_toTest[,s_toTest$SampleID] %>%
    melt() %>%
    mutate(value = value+1) %>%
    setNames(c("Taxa", "SampleID", "Abundance")) %>%
    merge(s_toTest, by="SampleID") %>%
    mutate(props = Abundance / read_counts) %>%
    mutate(props100 = props * 100) %>%
    mutate(props_logit = log(props/(1-props))) %>%
    group_by(Taxa) %>%
    nlme::lme(as.formula(form1), random = as.formula(rep_mes_form), data=., na.action=na.omit) %>%
    ungroup() %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value<p_cutoff)
}

#simple linear model functions
# tidy lm for the lmTest object
tidy_lm <- function(lm_test) {
  mod <- summary(lm_test)
  data.frame(term  = rownames(mod$coefficients), mod$coefficients, row.names=NULL)
}

# count based lm (for taxa abundances)
run_lm <- function(cts_toTest, s_toTest, form1, p_cutoff) {
  cts_toTest[,s_toTest$SampleID] %>%
    melt() %>%
    mutate(value = value+1) %>%
    setNames(c("Taxa", "SampleID", "Abundance")) %>%
    merge(s_toTest, by="SampleID") %>%
    mutate(props = Abundance / read_counts) %>%
    mutate(props100 = props * 100) %>%
    mutate(props_logit = log(props/(1-props))) %>%
    group_by(Taxa) %>%
    do(tidy_lm(lm(as.formula(form1), data=., na.action=na.omit))) %>%
    setNames(c("Taxa","term","Estimate","Std.Error","t.value","p.value")) %>%
    ungroup() %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value < p_cutoff)
}

# just get the lm itself
get_lm <- function(cts_toTest, s_toTest, form1) {
  cts_toTest[,s_toTest$SampleID] %>%
    melt() %>%
    mutate(value = value+1) %>%
    setNames(c("Taxa", "SampleID", "Abundance")) %>%
    merge(s_toTest, by="SampleID") %>%
    mutate(props = Abundance / read_counts) %>%
    mutate(props100 = props * 100) %>%
    mutate(props_logit = log(props/(1-props))) %>%
    group_by(Taxa) %>%
    lm(as.formula(form1), data=., na.action=na.omit)
}

# get a glm (as kyle says: It's like lm() but with a g in front of it!
get_glm <- function(cts_toTest, s_toTest, form1) {
  cts_toTest[,s_toTest$SampleID] %>%
    melt() %>%
    mutate(value = value+1) %>%
    setNames(c("Taxa", "SampleID", "Abundance")) %>%
    merge(s_toTest, by="SampleID") %>%
    mutate(props = Abundance / read_counts) %>%
    mutate(props100 = props * 100) %>%
    mutate(props_logit = log(props/(1-props))) %>%
    group_by(Taxa) %>%
    glm(as.formula(form1), data=., na.action=na.omit, family=binomial)
}

filter_low_coverage <- function(cts, perc_cutoff, min_ab=0){
  frac_nonzero <- function (x) sum(x > min_ab) / length(x)
  apply(cts, 1, frac_nonzero) >= perc_cutoff
}

se <- function(x) sd(x)/sqrt(length(x))

top_table <- function(summed_props, s_toPlot, thre=0.8, option=1, prop_cut=0.01){

  s_props <- summed_props[,s_toPlot$SampleID]

  if (option == 1) {
    rows_to_keep <- filter_low_coverage(s_props, frac_cutoff=thre)
  } else if (option == 2) {
    rows_to_keep <- apply(s_props,1,max) >= prop_cut
  }

  s_props <- as.data.frame(s_props[rows_to_keep,])

}

### clustered taxa heatmap
heatmap_cluster_taxa <- function(summed_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1, prop_cut=0.01, satu_limit=0.4){

  #color = saturated_rainbow(101)
  color = saturated_rainbow(101, saturation_limit=satu_limit)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))

  heatmap_props <- summed_props[,heatmap_s$SampleID]

  if (option == 1) {
    rows_to_keep <- filter_low_coverage(heatmap_props, frac_cutoff=thre)
  } else if (option == 2) {
    rows_to_keep <- apply(heatmap_props,1,max) >= prop_cut
  }
  heatmap_props <- heatmap_props[rows_to_keep,]

  ## group the SampleIDs
  heatmap_s %<>% arrange_(.dots=grps)
  heatmap_props <- heatmap_props[, heatmap_s$SampleID]

  ## update the annotation
  anno <- heatmap_s[,grps] %>% as.data.frame()
  rownames(anno) <- heatmap_s$SampleID
  colnames(anno) <- grps

  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = anno, color = color, breaks = breaks, filename = fname,
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = TRUE,cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = anno, color = color, breaks = breaks,
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = TRUE,cellheight = 8, cellwidth = 8)
}

change_data_format <- function(d) {
  #for dates like 20181121
  #changes to 11-21-2018 (but should be sorted on the 20181121 version)
  #we have 2019-03-29 and we want 03-29-2019
  paste(substr(d,6,7), substr(d,9,10), substr(d,1,4), sep="-")
}

#
#  make_pcoa_plot <- function(uu, s, shape_by, color_by, title)
#  uu: distance, s: mapping file, shape_by: variable used for shape, color_by: variable used for color
#

make_pcoa_plot <- function(dm, s, shape_by, color_by) {
  dm <- usedist::dist_subset(dm, s$SampleID)
  pc <- pcoa(dm)
  pc_df <- merge(s, pc$vectors[, 1:3], by.x="SampleID", by.y="row.names")
  pc_pct <- round(pc$values$Relative_eig * 100)

  pcoa_plot = ggplot(pc_df, aes(x=Axis.1, y=Axis.2)) +
    theme_bw() +
    scale_shape_discrete(name=sub("_", " ", shape_by)) +
    scale_colour_discrete(name=sub("_", " ", color_by)) +
    labs(
      x=paste0("PCoA axis 1 (", pc_pct[1], "%)"),
      y=paste0("PCoA axis 2 (", pc_pct[2], "%)")
    )

  if (is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by))))
  } else if (!is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by)), shape=factor(get(shape_by))))
  } else {
    pcoa_plot <- pcoa_plot + geom_point()
  }
  return(pcoa_plot)
}

heatmap_grouped <- function(summed_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1, prop_cut=0.01, satu_limit=0.4){

  #color = saturated_rainbow(101)
  color = saturated_rainbow(101, saturation_limit=satu_limit)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))

  heatmap_props <- summed_props[,heatmap_s$SampleID]

  if (option == 1) {
    rows_to_keep <- filter_low_coverage(heatmap_props, frac_cutoff=thre)
  } else if (option == 2) {
    rows_to_keep <- apply(heatmap_props,1,max) >= prop_cut
  }
  heatmap_props <- heatmap_props[rows_to_keep,]

  ## group the SampleIDs
  heatmap_s %<>% arrange_(.dots=grps)
  heatmap_props <- heatmap_props[, heatmap_s$SampleID]

  ## update the annotation
  annc <- heatmap_s[,grps] %>% as.data.frame()
  rownames(annc) <- heatmap_s$SampleID
  colnames(annc) <- grps

  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, filename = fname,
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks,
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
}
