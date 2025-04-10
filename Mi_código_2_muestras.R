#Calculate ChIP signals genomic regions
#Genome-wide plot
rm(list=ls())

#install.packages("patchwork")
#install.packages("foreach")
#install.packages("doParallel")
library(foreach)
library(doParallel)
library(dplyr)
library(reshape)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggrepel)


#Define genomic regions
genome_regions <<- c(84170, 84170+26264, 84170+26264+17780, 154477) #Genomic regions IR1start, IR1end, IR2start, IR2end

#Function that does most of the work in the script
#plot_plastid <- function(region_start, region_end, bin_size, correction){

#Load input files
input_files <- read.csv("/home/vtrejo/Desktop/Fig1B_ChIP_patterns/input_files.csv", stringsAsFactors = FALSE)                                     #CSV file with names of input files
loci_list <- read.csv("/home/vtrejo/Desktop/Fig1B_ChIP_patterns/classes_SF.csv", header = TRUE, check.names=FALSE)         #BED file with the list of annotated features

#Create genomic bins
region_start <- 0
region_end <- 154477
bin_size <- 200
correction <- "no"
number_of_bins <- (region_end - region_start) / bin_size + 1
bed_data <- data.frame(chr = matrix("Pt", nrow = number_of_bins, ncol = 1))
bed_data$start <- seq(region_start, region_end, bin_size)
bed_data$end <- bed_data$start + bin_size
write.table(bed_data, file = "/home/vtrejo/Desktop/Fig1B_ChIP_patterns/bins.bed", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

#PARALLEL VERSION
cores=detectCores()
cl <<- makeCluster(cores[1]-0) 
registerDoParallel(cl)
system.time(intersect_result <- foreach(input_dataset = 1:nrow(input_files), .combine=cbind) %dopar% {
  #Calculate read counts for RPM normalization
  tmp_wc1 <- as.numeric(system(paste0("wc -l < /home/vtrejo/Desktop/Fig1B_ChIP_patterns/datasets/", input_files[input_dataset,2], " | tr -d ' '"), intern = TRUE))
  tmp_wc2 <- as.numeric(system(paste0("wc -l < /home/vtrejo/Desktop/Fig1B_ChIP_patterns/datasets/", input_files[input_dataset,3], " | tr -d ' '"), intern = TRUE))
  
  #perform IntersectBed using system commands
  tmp_data1 <- system(paste0("intersectBed -a /home/vtrejo/Desktop/Fig1B_ChIP_patterns/bins.bed -b /home/vtrejo/Desktop/Fig1B_ChIP_patterns/datasets/", input_files[input_dataset,2], " -c -sorted"), intern = TRUE)
  tmp_data2 <- system(paste0("intersectBed -a /home/vtrejo/Desktop/Fig1B_ChIP_patterns/bins.bed -b /home/vtrejo/Desktop/Fig1B_ChIP_patterns/datasets/", input_files[input_dataset,3], " -c -sorted"), intern = TRUE)
  
  #    #Reformat outputs to data frames
  tmp_data1 <- read.table(text = tmp_data1, sep = "\t")
  tmp_data2 <- read.table(text = tmp_data2, sep = "\t")
  
  #Calcuate ratios, add to the data frame
  new_data <- data.frame((((tmp_data1$V4+0)/tmp_wc1)) / (((tmp_data2$V4+0)/tmp_wc2)) )
  colnames(new_data) <- input_files[input_dataset,1]
  new_data
})

stopCluster(cl)

bin_list <- cbind(bed_data, intersect_result)

#Generate data for averaging analysis
bin_list_averages <- bin_list[, c("start", "end")]
bin_list_sd <- bin_list[, c("start", "end")]

#Establish a list of genotypes
data_names <- data.frame(name = input_files$dataset)
data_names_separated <- separate(data_names, name, into = c("genotype", "replicate"), remove = FALSE, sep = "_")

#Average data for all genotypes
for(current_dataset in unique(data_names_separated$genotype)){
  data_current <- bin_list[ , grepl(current_dataset, colnames(bin_list))]     #select columns that contain current genotype in their names
  if(correction != "corrected"){
    bin_list_averages[[current_dataset]] <- apply(data_current, 1, mean)        #calculate averages of all available replicates
    bin_list_sd[[current_dataset]] <- apply(data_current, 1, sd)                #calculate sd of all available replicates
  }
  if(correction == "corrected"){
    if(current_dataset == "RpoB"){
      rpoBo_signal_current <- rpoBo_signals[ , grepl(current_dataset, colnames(rpoBo_signals))]
      rpoBo_signal_current_mean <- apply(rpoBo_signal_current, 1, mean)
      bin_list_averages[[current_dataset]] <- apply(data_current/rpoBo_signal_current_mean, 1, mean)
      bin_list_sd[[current_dataset]] <- apply(data_current/rpoBo_signal_current_mean, 1, sd)
    }else{
      bin_list_averages[[current_dataset]] <- apply(data_current, 1, mean)        #calculate averages of all available replicates

          }
  }
}


#Reorganize for ggplot2

bin_list_averages$mid <- bin_list_averages$start + (bin_list_averages$end - bin_list_averages$start) / 2
bin_list_sd$mid <- bin_list_sd$start + (bin_list_sd$end - bin_list_sd$start) / 2
bin_list_averages_melted <- melt(as.data.frame(bin_list_averages), id.vars = c("start", "mid", "end"))
bin_list_sd_melted <- melt(as.data.frame(bin_list_sd), id.vars = c("start", "mid", "end"))
bin_list_averages_melted$sd <- bin_list_sd_melted$value   #hay valores distintos
max_x <- max(bin_list_averages_melted$value, na.rm=TRUE)
min_x <- min(bin_list_averages_melted$value, na.rm=TRUE)
height_x <- max_x - min_x
loci_list$strand <- ifelse(loci_list$strand == "+", 1, -1)
loci_list_label <- loci_list[loci_list$name %in% c("psbA", "trnG1_1", "rpoB", "psbC", "psaA", "trnF", "rbcL", "psbE", "psbB", "trnI1", "ftsHi", "trnL2",  "rrn16", "trnL3" ), ]

#Se añade la columna condicional "change"
loci_list_label <- loci_list_label %>%
  mutate(change = ifelse(strand != lag(strand, default = first(strand)), "YES", "NO"))

# vertical_lines <- loci_list %>%
#   filter(change == "YES") %>%
#   select(Position = start)  # Cambia 'start' si necesitas otra columna


p1 <- ggplot() +
  geom_segment(aes(x = 1, y = max_x + height_x * 0.2, xend = genome_regions[1], yend = max_x + height_x * 0.2), lineend = "round", linejoin = "mitre", size = 1.8, color = "#8DD3C7") +
  geom_segment(aes(x = genome_regions[1], y = max_x + height_x * 0.2, xend = genome_regions[2], yend = max_x + height_x * 0.2), lineend = "round", linejoin = "mitre", size = 1.8, color = "#FB8072") +
  geom_segment(aes(x = genome_regions[2], y = max_x + height_x * 0.2, xend = genome_regions[3], yend = max_x + height_x * 0.2), lineend = "round", linejoin = "mitre", size = 1.8, color = "#BEBADA") +
  geom_segment(aes(x = genome_regions[4], y = max_x + height_x * 0.2, xend = genome_regions[3], yend = max_x + height_x * 0.2), lineend = "round", linejoin = "mitre", size = 1.8, color = "#FB8072") +
  
  geom_segment(aes(x = 1, y = max_x + height_x * 0.2 - height_x * 0.04, xend = 1, yend = max_x + height_x * 0.2 + height_x * 0.04), lineend = "round", linejoin = "mitre", size = 1.8, color = "#bdbdbd") +
  geom_segment(aes(x = genome_regions[1], y = max_x + height_x * 0.2 - height_x * 0.04, xend = genome_regions[1], yend = max_x + height_x * 0.2 + height_x * 0.04), lineend = "round", linejoin = "mitre", size = 1.8, color = "#bdbdbd") +
  geom_segment(aes(x = genome_regions[2], y = max_x + height_x * 0.2 - height_x * 0.04, xend = genome_regions[2], yend = max_x + height_x * 0.2 + height_x * 0.04), lineend = "round", linejoin = "mitre", size = 1.8, color = "#bdbdbd") +
  geom_segment(aes(x = genome_regions[3], y = max_x + height_x * 0.2 - height_x * 0.04, xend = genome_regions[3], yend = max_x + height_x * 0.2 + height_x * 0.04), lineend = "round", linejoin = "mitre", size = 1.8, color = "#bdbdbd") +
  geom_segment(aes(x = genome_regions[4], y = max_x + height_x * 0.2 - height_x * 0.04, xend = genome_regions[4], yend = max_x + height_x * 0.2 + height_x * 0.04), lineend = "round", linejoin = "mitre", size = 1.8, color = "#bdbdbd") +
  
  geom_hline(yintercept = (max_x + height_x * 0.1), color = "#bdbdbd", show.legend=FALSE) +
  geom_rect(data=loci_list[loci_list$gene == "rRNA", ], aes(xmin=start, xmax=end, ymin=(max_x + height_x * 0.1), ymax=(max_x + height_x * 0.1 + height_x * 0.03 * strand)), fill = "#FC8D62", color = "black", size=0.1, show.legend=FALSE) +
  geom_rect(data=loci_list[loci_list$gene == "tRNA", ], aes(xmin=start, xmax=end, ymin=(max_x + height_x * 0.1), ymax=(max_x + height_x * 0.1 + height_x * 0.03 * strand)), fill = "#E78AC3", color = "black", size=0.1, show.legend=FALSE) +
  geom_rect(data=loci_list[loci_list$gene == "protein", ], aes(xmin=start, xmax=end, ymin=(max_x + height_x * 0.1), ymax=(max_x + height_x * 0.1 + height_x * 0.03 * strand)), fill = "#66C2A5", color = "black", size=0.1, show.legend=FALSE) +
  #geom_rect(data=operons, aes(xmin=V1, xmax=V2, ymin=(max_x + height_x * 0.1), ymax=(max_x + height_x * 0.1)), color = "#CC79A7", size=1.5, show.legend=FALSE) +
  geom_text(data=loci_list_label, aes(x=start+(end-start)/2, y=(max_x + height_x * 0.1 + height_x * 0.06 * strand), label=name), size=2, show.legend=FALSE, check_overlap = FALSE, fontface = "italic") +
  
  geom_hline(yintercept = 1, color = "#0a0a0a", show.legend=FALSE) +
  geom_hline(yintercept = 1.0941971, color = "#0a0a0a", show.legend=FALSE, size=0.3, linetype="dashed") +
  geom_hline(yintercept = 0.7622402, color = "#0a0a0a", show.legend=FALSE, size=0.3, linetype="dashed") +
  # Línea roja para las filas donde "variable" es "Day"
  geom_line(data = subset(bin_list_averages_melted), 
            aes(x = mid, y = value, group = variable, color = variable), size = 1.0) + 
  # Línea azul para las filas donde "variable" es "Night"
  #geom_line(data = subset(bin_list_averages_melted, variable == "Psoralen"), 
  #          aes(x = mid, y = value, group = variable, color = "Psoralen"), size = 0.6) + 
  # Relleno de área (ribbon) con color rojo para "Day"
  #geom_ribbon(data = subset(bin_list_averages_melted, variable == "GapR"), 
  #            aes(x = mid, ymax = value + sd, ymin = value - sd), alpha = 0.2) + 
  # Relleno de área (ribbon) con color azul para "Night"
  geom_ribbon(data = subset(bin_list_averages_melted), 
              aes(x = mid, ymax = value + sd, ymin = value - sd, fill = variable), alpha = 0.2) +
  # Configuración de colores personalizados en la leyenda
  scale_color_manual(values = c("Day" = "red", "RpoB" = "blue")) +
  scale_fill_manual(values = c("Day" = "red", "RpoB" = "blue")) +
  # Etiquetas y otros ajustes (si es necesario)
    labs(title = "Psoralen_Day_Night 200 bins sesión 2", x = "Mid", y = "Value") +
  
  coord_cartesian(ylim=c(min_x, max_x + height_x * 0.2), xlim=c(region_start, region_end)) +    
  #scale_x_continuous(expand = c(0, 0)) +  
  #scale_y_continuous(name = "ChIP enrichment", limits = c(max_RNAseq_minus/ChIP_RNAseq_ratio, 1.6 * max_x), sec.axis = sec_axis(~.*ChIP_RNAseq_ratio*1.6, breaks = c(max_RNAseq_plus, 0), name = "TSS RNA-seq", labels = label_scientific(digits=1, trim = FALSE) )) +
  #scale_color_brewer(palette = "Set1") +
  #scale_fill_brewer(palette = "Set1") +
  ylab("Unión de Psoraleno") +
  #ggtitle("WT_GapR_sig2") +
  geom_segment(data = loci_list_label[loci_list_label$change == "YES", ], 
               aes(x = start, xend = start, 
                   y = min_x, yend = max_x + height_x * 0.1), 
               color = "4B974F", size = 0.3, linetype = "dashed") +
  theme(text = element_text(size=15), 
        legend.title = element_blank(), axis.title.x=element_blank(),
        legend.position = c(.05, .80), legend.justification = c("left", "top"), legend.box.just = "left", legend.margin = margin(6, 6, 6, 6),
        legend.key = element_rect(colour = NA, fill = NA), legend.background=element_rect(colour = NA, fill = NA),
        panel.background = element_rect(fill = 'white', colour = 'white'))#, plot.margin = margin(-1, 0, 0, 0, "cm") )



print(p1)

ggsave("/home/vtrejo/Desktop/Fig1B_ChIP_patterns/MisPlots/Figura_13_con_quartiles_1st.pdf", p1, width = 20, height = 8, units = "cm")


#For BoxPlot


library(gdata)
newdata_low <- subset(cociente, cociente$cociente_day_night < 0.9298278, select=c(start,end,cociente_day_night,GapR))
newdata_med <- subset(cociente, cociente$cociente_day_night > 0.9298278 & cociente$cociente_day_night < 1.0232138, select=c(start,end,cociente_day_night,GapR))
newdata_high <- subset(cociente, cociente$cociente_day_night > 1.0232138, select=c(start,end,cociente_day_night,GapR))


newdata= cbindX(newdata_low,newdata_med,newdata_high)


names(newdata)[3] = "psoralen_low"
names(newdata)[7] = "psoralen_med"
names(newdata)[11] = "psoralen_high"

df_long <- newdata %>%
  select(psoralen_low,psoralen_med,psoralen_high)

boxplot(df_long, main="Cociente Psoraleno Día/Noche vs GapR_Noche")
#, xlab="X axis title"
#, ylab="Y axis title"
#, sub="Sub-title")

names(newdata)[4] = "GapR1"
names(newdata)[8] = "GapR2"
names(newdata)[12] = "GapR3"

df_long1 <- newdata %>%
  select(GapR1, GapR2, GapR3)

boxplot(df_long1, main="GapR")
