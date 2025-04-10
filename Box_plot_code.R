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
bin_size <- 100
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
      bin_list_sd[[current_dataset]] <- apply(data_current, 1, sd)
    }
  }
}

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


cociente <- data.frame(bin_list_averages$start,bin_list_averages$end,bin_list_averages$Day/bin_list_averages$Night,bin_list_averages$GapR)
colnames(cociente) <- c("start", "end", "cociente_day_night", "GapR")
quantile(cociente$cociente_day_night)

bin_list_averages$cociente_day_night <- bin_list_averages$Day/bin_list_averages$Night


library(gdata)
newdata_low <- subset(cociente, cociente$cociente_day_night < 0.8897250, select=c(start,end,cociente_day_night,GapR))
newdata_med <- subset(cociente, cociente$cociente_day_night > 0.8897250  & cociente$cociente_day_night < 1.0487395, select=c(start,end,cociente_day_night,GapR))
newdata_high <- subset(cociente, cociente$cociente_day_night > 1.0487395, select=c(start,end,cociente_day_night,GapR))
  #newdata_ribosome <- subset(cociente, cociente$start > 110434 & cociente$start < 128214, select=c(start,end,cociente_day_night,GapR))


newdata= cbindX(newdata_low,newdata_med,newdata_high)
  #newdata= cbindX(newdata_low,newdata_med,newdata_high,newdata_ribosome)

  
names(newdata)[3] = "Psoraleno_Baja"
names(newdata)[7] = "Psoraleno_Media"
names(newdata)[11] = "Psoraleno_Alta"
  #names(newdata)[15] = "ribosome operons"
  
  #newdata_filtered <- newdata
  #duplicados_high_ribo <- newdata$psoralen_high %in% newdata$`ribosome operons` # Identificar valores duplicados entre psoralen_high y ribosome operons
  #newdata_filtered$psoralen_high[duplicados_high_ribo] <- NA # Eliminar los valores duplicados en psoralen_high

df_long <- newdata %>%
  select(Psoraleno_Baja,Psoraleno_Media,Psoraleno_Alta)  # Seleccionar las columnas necesarias

boxplot(df_long, main="psoraleno categories (100 bins) according to quantiles")
#, xlab="X axis title"
#, ylab="Y axis title"
#, sub="Sub-title")

names(newdata)[4] = "psoralen_lowR"
names(newdata)[8] = "psoralen_medR"
names(newdata)[12] = "psoralen_highR"
  #names(newdata)[16] = "ribsosome_operons_R"

df_long1 <- newdata %>%
  select(psoralen_lowR,psoralen_medR,psoralen_highR)

boxplot(df_long1, main="GapR vs psoraleno (100 bins) according to quantiles")



#Para filtrar los datos por a) baja Desviación Estándar y b) unión de psoraleno alta/baja 

bin_list_averages$sd_Day <- bin_list_sd$Day
bin_list_averages$sd_Night <- bin_list_sd$Night

quantile(bin_list_averages$sd_Day)
quantile(bin_list_averages$sd_Night)
newdata_ds_baja <- subset(bin_list_averages, bin_list_averages$sd_Day < 0.1 & bin_list_averages$sd_Night < 0.1, select=c(start, end, cociente_day_night, GapR))
#El 0.1 se elige de acuerdo a los valores de los cuantiles; hay que descartar una DS alta

quantile(newdata_ds_baja$cociente_day_night)

#Para el gráfico de cajas y bigotes con los dos grupos
library(gdata)
newdata_lowx2 <- subset(newdata_ds_baja, newdata_ds_baja$cociente_day_night < 0.9098363, select=c(start,end,cociente_day_night,GapR))
newdata_medx2 <- subset(newdata_ds_baja, newdata_ds_baja$cociente_day_night < 0.9098363 & newdata_ds_baja$cociente_day_night > 1.25, select=c(start,end,cociente_day_night,GapR))
newdata_highx2 <- subset(newdata_ds_baja, newdata_ds_baja$cociente_day_night > 1.25, select=c(start,end,cociente_day_night,GapR))

newdata= cbindX(newdata_lowx2,newdata_medx2,newdata_highx2)
names(newdata)[3] = "psoralen_low"
names(newdata)[7] = "psoralen_med"
names(newdata)[11] = "psoralen_high"

df_long <- newdata %>%
  select(psoralen_low,psoralen_high)

boxplot(df_long, main="Psoraleno filtrado 300 bins 0.9098363 y 1.25")

names(newdata)[4] = "psoralen_lowR"
names(newdata)[8] = "psoralen_medR"
names(newdata)[12] = "psoralen_highR"

df_long1 <- newdata %>%
  select(psoralen_lowR, psoralen_highR)

boxplot(df_long1, main="GapR filtrado 300 bins 0.9098363 y 1.25")


#Gráfico de correlación 
library(ggplot2)

umbral_bajo <- 0.8897250
umbral_alto <- 1.0487395

cociente$grupo <- ifelse(cociente$cociente_day_night < umbral_bajo, "unión_baja",
                   ifelse(cociente$cociente_day_night > umbral_alto, "unión_alta", "unión_media"))

cor_value <- cor(cociente$cociente_day_night, cociente$GapR, method = "spearman", use = "complete.obs")

ggplot(cociente, aes(x = cociente_day_night, y = GapR, color = grupo)) +
  geom_point(alpha = 0.6) +  # Agregar puntos al gráfico con transparencia
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Línea de regresión
  scale_color_manual(values = c("unión_baja" = "green", "unión_alta" = "blue", "unión_media" = "gray")) +  # Colores personalizados
  labs(x = "Cociente Day/Night", y = "GapR", title = "Correlación entre Cociente Day/Night y GapR",
       color = "Categoría de Unión") +  # Etiqueta de la leyenda
  theme_minimal() +  # Estilo limpio del gráfico
  annotate("text", x = max(cociente$cociente_day_night, na.rm = TRUE), 
           y = max(cociente$GapR, na.rm = TRUE), 
           label = paste("r (Spearman) =", round(cor_value, 3)), 
           hjust = 1, vjust = 1, size = 5, color = "black") + # Agregar la correlación
           theme(legend.position = c(0.85, 0.75))  # Posición de la leyenda dentro del gráfico


#Líneas de tendencia por grupo 
cor_values <- cociente %>%
  group_by(grupo) %>%
  summarise(cor_spearman = round(cor(cociente_day_night, GapR, method = "spearman"), 3))

x_pos <- max(cociente$cociente_day_night, na.rm = TRUE) - 0.2  # Ligeramente a la izquierda del borde derecho
y_min <- min(cociente$GapR, na.rm = TRUE)  # Posición base en la parte inferior
y_step <- (max(cociente$GapR, na.rm = TRUE) - y_min) * 0.05  # Espaciado vertical entre textos

# Crear el gráfico
ggplot(cociente, aes(x = cociente_day_night, y = GapR)) +
  geom_point(aes(color = grupo), alpha = 0.6) +
  geom_smooth(aes(linetype = grupo), method = "lm", color = "red", se = FALSE) +  # líneas rojas
  scale_color_manual(
    values = c("unión_baja" = "green", "high" = "blue", "unión_media" = "gray")
  ) +
  labs(
    x = "Cociente Day/Night",
    y = "GapR",
    title = "Correlación entre Cociente Day/Night y GapR",
    color = "Categoría de Unión",
    linetype = "Grupo"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.90),
    legend.text = element_text(size = 13),       # tamaño de texto en la leyenda aumentado
    legend.title = element_text(size = 14, face = "bold")  # título de la leyenda más grande y en negritas (opcional)
  ) +
  geom_text(data = cor_values, aes(
    x = x_pos, 
    y = y_min + (as.numeric(factor(grupo)) - 1) * y_step,
    label = paste(grupo, ": r =", cor_spearman)
  ),
  hjust = 1, vjust = 0, size = 5, color = "black")

