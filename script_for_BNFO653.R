##### install packages #####
ip <- as.data.frame(installed.packages())
ip <- ip$Package

if (sum(ip == "rstudioapi") == 0) {
  install.packages("rstudioapi")
}

if (sum(ip == "vegan") == 0) {
  install.packages("vegan")
}

if (sum(ip == "stringr") == 0) {
  install.packages("stringr")
}


if (sum(ip == "GUniFrac") == 0) {
  install.packages("GUniFrac")
}

if (sum(ip == "tidyr") == 0) {
  install.packages("tidyr")
}


if (sum(ip == "ggplot2") == 0) {
  install.packages("ggplot2")
}

if (sum(ip == "pheatmap") == 0) {
  install.packages("pheatmap")
}



##### library #####
library(vegan) # for diversity
library(stringr) 
library(ggplot2)
library(GUniFrac) 
library(tidyr)
library(pheatmap)
##### set files path #####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

##### Read in your counts table #####
count_table = read.table('Galaxy104-[Count.seqs_on_data_97_and_data_103__count_table].mothur.count_table',sep = '\t',
                         header = T)
count_table = count_table[,-2]

align_table = read.table('Galaxy106-[Align.seqs_on_data_42_and_data_102__align.report].mothur.align.report',sep = '\t',
                         header = T)

taxonomy = read.table('Silva_reference_files_Mothur.txt',sep = '\t')
taxonomy$V2 = str_replace_all(taxonomy$V2,'\\;','.')

metadata = read.csv('metadata.csv',header = T)
##### prepare reads table and get taxa name #####
count_table = count_table[,-1]

taxa = align_table[,3]

unique_taxa = unique(taxa)

reads_table = as.data.frame(matrix(data = 0, nrow = length(unique_taxa), ncol = ncol(count_table)))
row.names(reads_table) = unique_taxa
colnames(reads_table) = colnames(count_table)

taxa_name = vector(length = length(unique_taxa))
  
for (a in 1:length(unique_taxa)) {
  n = which(taxa == unique_taxa[a])
  reads_table[a,] = colSums(count_table[n,])
  
  m = which(taxonomy$V1 == row.names(reads_table)[a])
  taxa_name[a] = taxonomy$V2[m]
}

######### alpha diversity ### alpha_diversity ### input samples in cols ###########
alpha_diversity = function(reads_table, metadata = NA, factor_name = NA, paired = F,order = NA, rarefy_to = NA) {
  
  if (is.na(metadata)[1]) {
    print('no metadata')
    return(NA)
  }
  
  if (is.na(factor_name)) {
    print('no factor name')
    return(NA)
  }
  
  # rarefy to normalize data
  reads_table = t(reads_table)
  
  if (is.na(rarefy_to)) {
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  } else {
    reads_table = Rarefy(reads_table, depth = rarefy_to)
  }
  
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  # calculate diversity
  alpha.shannon_diversity <- data.frame(diversity(reads_table))
  alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table))
  alpha.ovserved_OTU <- data.frame(colSums(t(reads_table) != 0))
  
  alpha = as.data.frame(matrix(data = NA,ncol=3,nrow = nrow(reads_table)))
  colnames(alpha) = c('alpha.shannon','alpha.evenness','alpha.ovserved_OTU')
  
  alpha$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table.
  alpha$alpha.evenness <- alpha.evenness$diversity.reads_table.
  alpha$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table.....0.
  
  metadata = cbind(metadata, alpha)
  
  colnames(metadata)[1] = 'factor'
  
  metadata = metadata[order(metadata$factor),]
  
  if (is.na(order)[1] ) {
    metadata$factor <- factor(metadata$factor , levels = unique(metadata$factor))
  } else {
    metadata$factor <- factor(metadata$factor , levels = order)
  }
  
  
  alpha.shannon = ggplot(metadata, aes(x=factor, y=alpha.shannon)) +
    geom_boxplot(aes(fill = factor),outlier.shape=NA) +
    labs(x = NULL, y = "Shannon index", fill=factor_name)+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          axis.text.x=element_blank(),
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  
  alpha.evenness = ggplot(metadata, aes(x=factor, y=alpha.evenness)) +
    geom_boxplot(aes(fill = factor),outlier.shape=NA) +
    labs(x = NULL, y = "Evenness", fill=factor_name)+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          axis.text.x=element_blank(),
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  
  alpha.ovserved_OTU = ggplot(metadata, aes(x=factor, y=alpha.ovserved_OTU)) +
    geom_boxplot(aes(fill = factor),outlier.shape=NA) +
    labs(x = NULL, y = "Observed OTU", fill=factor_name)+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          axis.text.x=element_blank(),
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  # geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5)+
  
  # calculate significance
  factor_levels = unique(metadata$factor)
  n = length(factor_levels)
  
  Shannon_sig = as.data.frame(matrix(data = NA, nrow =n, ncol = n))
  colnames(Shannon_sig) = factor_levels
  row.names(Shannon_sig) = factor_levels
  
  Evenness_sig = as.data.frame(matrix(data = NA, ncol=n, nrow = n))
  colnames(Evenness_sig) = factor_levels
  row.names(Evenness_sig) = factor_levels
  
  OTU_sig = as.data.frame(matrix(data = NA, ncol=n, nrow = n))
  colnames(OTU_sig) = factor_levels
  row.names(OTU_sig) = factor_levels
  
  for (a in 1:(n-1)) {
    for (b in (a+1) : n) {
      factor_level1 <- subset(metadata,  factor == factor_levels[a],
                              drop = TRUE)
      factor_level2 <- subset(metadata,  factor == factor_levels[b],
                              drop = TRUE)
      
      Shannon_sig[a,b] <- wilcox.test(factor_level1$alpha.shannon, 
                                      factor_level2$alpha.shannon, paired = paired)$p.value
      Evenness_sig[a,b] <- wilcox.test(factor_level1$alpha.evenness, 
                                       factor_level2$alpha.evenness, paired = paired)$p.value
      OTU_sig[a,b] <- wilcox.test(factor_level1$alpha.ovserved_OTU, 
                                  factor_level2$alpha.ovserved_OTU, paired = paired)$p.value
      
    }
  }
  output = c(list(alpha = alpha), list(shannon = alpha.shannon), 
             list(evenness =alpha.evenness) , list(ovserved_OTU =alpha.ovserved_OTU),
             list(sig_Shannon = Shannon_sig),list(sig_Evenness = Evenness_sig),
             list(sig_OTU = OTU_sig))
  
  return(output)
}




##### alpha diversity #####
alpha = alpha_diversity(reads_table, metadata = metadata$Group, factor_name = 'Group', paired = F,order = NA, rarefy_to = NA)
alpha$shannon
ggsave('shannon.pdf',,width=2, height=2.5)
alpha$evenness
ggsave('evenness.pdf',,width=2, height=2.5)
alpha$ovserved_OTU
ggsave('ovserved_OTU.pdf',,width=2, height=2.5)
write.table(alpha$sig_Shannon,'alpha.sig_Shannon.txt')
write.table(alpha$sig_Evenness,'alpha.sig_Evenness.txt')
write.table(alpha$sig_OTU,'alpha.sig_OTU.txt')
######### beta diversity ### beta_diversity ### input samples in cols; metadata and factor_name are needed; order of factors could be set; ref_group is for setting the reference of bc distance in different groups; can skip from NMDS; output bc distance, within sample distance, distance amoung groups and NMDS #####
beta_diversity = function(reads_table, metadata = NA, factor_name = NA, order = NA, NMDS_skip = T, ref_group = NA, rarefy_to = NA) {
  #  metadata = metadata$Flag
  #  factor_name = 'Flag'
  #  order = NA
  #  NMDS_skip = T
  #  ref_group = NA
  
  if (is.na(metadata)[1]) {
    print('no metadata')
    return(NA)
  }
  
  if (is.na(factor_name)[1]) {
    print('no factor name')
    return(NA)
  }
  
  # rarefy to normalize data
  reads_table = as.data.frame(t(reads_table))
  
  if (is.na(rarefy_to)) {
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  } else {
    reads_table = Rarefy(reads_table, depth = rarefy_to)
  }
  
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  metadata=as.matrix(metadata)
  
  # Bray_Curtis
  Bray_Curtis <- as.matrix(vegdist(reads_table, method = "bray", binary=FALSE))
  Bray_Curtis <- as.data.frame(Bray_Curtis)
  
  Bray_Curtis_2 = Bray_Curtis
  Bray_Curtis_2[row(Bray_Curtis_2) <= col(Bray_Curtis_2)] =NA
  
  # within sample distance
  group_dis = gather(Bray_Curtis_2)
  group_dis$key2 = rep(row.names(Bray_Curtis_2),ncol(Bray_Curtis_2))
  
  Source = matrix(data = NA, ncol = length(metadata), nrow = length(metadata))
  
  for (a in 1:length(metadata)) {
    Source[a,] = metadata
  }
  Source = gather(as.data.frame(Source))
  group_dis$Source = Source$value
  group_dis$Target = rep(metadata,length(metadata))
  
  group_dis = group_dis[!is.na(group_dis$value),]
  group_dis$Source <- as.factor(group_dis$Source)
  #  group_dis$value = as.numeric(as.character(group_dis$value))
  
  keep = group_dis$Source == group_dis$Target
  within_dis = group_dis[keep,]
  keep = within_dis$key != within_dis$key2
  within_dis = within_dis[keep,]
  #  within_dis$value = as.numeric(as.character(within_dis$value))
  
  if (!is.na(order)) {
    within_dis$Source = as.factor(within_dis$Source)
    within_dis$Source = factor(within_dis$Source, levels= order)
  }
  
  within_dis_p = ggplot(within_dis, aes(x=Source, y=value)) +
    geom_boxplot(aes(fill = Source),outlier.shape=NA) +
    labs(x = NULL, y = "Within sample distance", fill=factor_name)+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
  
  # significance among within sample distance, Wilcoxon test
  group_level = unique(within_dis$Source)
  n = length(group_level)
  
  within_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
  colnames(within_dis_sig) = group_level
  row.names(within_dis_sig) = group_level
  for (a in 1:(n-1)) {
    for (b in (a+1): n) {
      set1 <- as.matrix(subset(within_dis,  Source == group_level[a], value,
                               drop = F))
      set2 <- as.matrix(subset(within_dis,  Source == group_level[b], value,
                               drop = F))
      within_dis_sig[a,b] <- wilcox.test(set1, set2, paired = F)$p.value
      
    }
  }
  
  # distance among groups
  if (!is.na(ref_group)) {
    keep = group_dis$Source == ref_group
    group_dis_2 = group_dis[keep,]
    
    if (!is.na(order)) {
      group_dis_2$Target = as.factor(group_dis_2$Target)
      group_dis_2$Target = factor(group_dis_2$Target, levels= order)
    }
    
    group_dis_2_p = ggplot(group_dis_2, aes(x=Target, y=value)) +
      geom_boxplot(aes(fill = Target),outlier.shape=NA) +
      labs(x = NULL, y = paste0("Distance to ",ref_group))+ 
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))
    
    # significance of gourp sample distance, Wilcoxon test
    group_level = unique(group_dis_2$Target)
    n = length(group_level)
    
    group_dis_sig <- matrix(data = NA, ncol=n, nrow = 1)
    colnames(group_dis_sig) = group_level
    
    set1 <- as.matrix(subset(group_dis_2,  Target == group_level[group_level == ref_group], value,
                             drop = F))
    
    for (a in 1:n) {
      set2 <- as.matrix(subset(group_dis_2,  Target == group_level[a], value,
                               drop = F))
      group_dis_sig[a] <- wilcox.test(set1, set2, paired = F)$p.value
      
    }
    
  } else {
    group_level = unique(group_dis$Source)
    n = length(group_level)
    distance_median = matrix(data=NA, nrow = n, ncol =n)
    colnames(distance_median) = group_level
    row.names(distance_median) = group_level
    
    group_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
    colnames(group_dis_sig) = group_level
    row.names(group_dis_sig) = group_level
    
    for (a in 1: n) {
      set1 <- group_dis$value[group_dis$Source == row.names(group_dis_sig)[a] & group_dis$Target == row.names(group_dis_sig)[a]]
      for (b in a:n) {
        distance_data = group_dis$value[group_dis$Source == row.names(distance_median)[a] & group_dis$Target == colnames(distance_median)[b]]
        distance_median[a,b] <- median(distance_data)
        distance_median[b,a] <- distance_median[a,b]
        
        group_dis_sig[a,b] <- wilcox.test(set1, distance_data, paired = F)$p.value
        group_dis_sig[b,a] <- group_dis_sig[a,b]
      }
    }
    
    #   group_dis_sig_adj = adjust.p(group_dis_sig, pi0.method="bky", alpha = alpha, pz = 0.05)
    #   group_dis_sig_adj = group_dis_sig_adj$adjp
    #   group_dis_sig_adj = group_dis_sig_adj$adjusted.p
    
    group_dis_sig_2 = group_dis_sig
    group_dis_sig_2[group_dis_sig > 0.05] = ''
    group_dis_sig_2[group_dis_sig <= 0.05 & group_dis_sig > 0.001] = '*'
    group_dis_sig_2[group_dis_sig <= 0.001 & group_dis_sig > 0.0001] = '**'
    group_dis_sig_2[group_dis_sig <= 0.0001] = '***'
    
    group_dis_2_p = pheatmap(distance_median, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=T, show_colnames=T, 
                             color=colorRampPalette(c("red","white"))(100),
                             fontsize =7, display_numbers = group_dis_sig_2)
  }
  
  
  #  group_dis_sig_adj = adjust.p("group_dis_sig_adj", pi0.method="bky", alpha = alpha,pz = 0.05)
  
  
  
  
  
  # Running Nonmetric Multidimensional Scaling (NMDS) Ordination
  if (NMDS_skip == T) {
    # output
    output = c(Bray_Curtis = list(Bray_Curtis), within_dis_p = list(within_dis_p),
               within_dis_sig = list(within_dis_sig), group_dis_2_p = list(group_dis_2_p),
               group_dis_sig = list(group_dis_sig))
    return(output)
    
  } else {
    
    colnames(metadata)[1] = 'factor'
    
    NMDS <-
      metaMDS(Bray_Curtis,
              distance = "bray",
              k = 2,
              maxit = 999, 
              trymax = 20,
              wascores = TRUE)
    
    mds_data <- as.data.frame(NMDS$points)
    mds_data$factor <- metadata
    
    if (!is.na(order)[1]) {
      mds_data$factor = as.factor(mds_data$factor)
      mds_data$factor = factor(mds_data$factor, levels= order)
    }
    
    NMDS = ggplot(mds_data, aes(x = MDS1, y = MDS2, color = factor)) +
      geom_point()+
      scale_colour_discrete(factor_name)+
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))
    
    #+stat_ellipse(type = "t")
    
    
    # output
    output = c(Bray_Curtis = list(Bray_Curtis), within_dis_p = list(within_dis_p), 
               within_dis_sig = list(within_dis_sig),group_dis_2_p = list(group_dis_2_p),
               group_dis_sig = list(group_dis_sig), NMDS =list(NMDS))
    return(output)
  }
  
}







##### beta diversity #####
beta = beta_diversity(reads_table, metadata = metadata$Group, factor_name = 'Group', order = NA, NMDS_skip = F, ref_group = 'X', rarefy_to = NA)
beta$group_dis_2_p
ggsave('beta.pdf',,width=2, height=2.5)
beta$NMDS
ggsave('NMDS.pdf',,width=2, height=3)
