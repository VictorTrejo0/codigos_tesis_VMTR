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
correction <- "corrected"
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


# Calcular la mediana de los datos
median_cociente <- median(cociente$cociente_day_night)

# Aplicar la prueba de Wilcoxon a TODO el conjunto de datos
wilcoxon_test <- wilcox.test(cociente$cociente_day_night, mu = median_cociente, exact = FALSE)

# Obtener el p-valor global
p_value_global <- wilcoxon_test$p.value
print(paste("P-valor de la prueba de Wilcoxon:", p_value_global))

# Definir nivel de significancia
alpha <- 0.05

# Filtrar outliers manualmente usando la distancia a la mediana
iqr_value <- IQR(cociente$cociente_day_night)  # Rango intercuartil
q1 <- quantile(cociente$cociente_day_night, 0.25)  # Primer cuartil
q3 <- quantile(cociente$cociente_day_night, 0.75)  # Tercer cuartil

# Limites para detectar outliers (método de Tukey)
lower_bound <- q1 - 1.5 * iqr_value
upper_bound <- q3 + 1.5 * iqr_value

# Identificar outliers
outliers_wilcoxon <- cociente[
  cociente$cociente_day_night < lower_bound | cociente$cociente_day_night > upper_bound, ]

# Filtrar datos sin outliers
filtered_data <- cociente[
  cociente$cociente_day_night >= lower_bound & cociente$cociente_day_night <= upper_bound, ]

# Mostrar resultados
print("Valores identificados como outliers:")
print(outliers_wilcoxon)


#Para generar data.frames con los datos filtrados (sin outliers) y solo con outliers
# Calcular el IQR (Rango Intercuartil) y los límites para outliers
iqr_value <- IQR(cociente$cociente_day_night)  
q1 <- quantile(cociente$cociente_day_night, 0.25)  
q3 <- quantile(cociente$cociente_day_night, 0.75)  

# Límites de Tukey para identificar outliers
lower_bound <- q1 - 1.5 * iqr_value
upper_bound <- q3 + 1.5 * iqr_value

# Data.frame sin outliers
bin_list_filtered <- cociente[
cociente$cociente_day_night >= lower_bound & 
cociente$cociente_day_night <= upper_bound, ]

# Data.frame solo con outliers
bin_list_outliers <- cociente[
cociente$cociente_day_night < lower_bound | 
cociente$cociente_day_night > upper_bound, ]

# Verificar tamaños
print(paste("Número total de datos:", nrow(cociente)))
print(paste("Datos sin outliers:", nrow(bin_list_filtered)))
print(paste("Datos considerados outliers:", nrow(bin_list_outliers)))



  #Crear los boxplot comparativos
  library(ggplot2)

  # Gráfico con outliers
  boxplot_original <- ggplot(bin_list_averages, aes(y = cociente_day_night)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red", outlier.shape = 16, outlier.size = 3) +
  labs(title = "Distribución con Outliers", y = "Cociente Día/Noche") +
  theme_minimal()

  # Gráfico sin outliers
  boxplot_sin_outliers <- ggplot(bin_list_filtered, aes(y = cociente_day_night)) +
  geom_boxplot(fill = "lightgreen", color = "black") +
  labs(title = "Distribución sin Outliers", y = "Cociente Día/Noche") +
  theme_minimal()

  # Mostrar ambos gráficos
  print(boxplot_original)
  print(boxplot_sin_outliers)
  
#BoxPlot del grupo sin outliers vs el grupo de solo outliers
library(ggplot2)

cociente$grupo <- ifelse(cociente$cociente_day_night %in% bin_list_outliers$cociente_day_night, "outlier", "normal")
  
ggplot(cociente, aes(x = grupo, y = cociente_day_night, fill = grupo)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 3) +
  scale_fill_manual(
    name = "Sub-grupo",
    values = c("outlier" = "red", "normal" = "lightblue"),
    labels = c("outlier" = "Valores Atípicos", "normal" = "Valores Normales")
  ) +
  labs(
    title = "Psoraleno",
    x = "Categoría",
    y = "Cociente Psoraleno_Día/Psoraleno_Noche"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # título centrado y más grande
    axis.title.y = element_text(size = 25),
    axis.title.x = element_text(size= 25),
    legend.text = element_text(size = 14),          # texto de la leyenda más grande
    legend.title = element_text(size = 16, face = "bold")  # título de la leyenda también
  )

    

#Para GapR
ggplot(cociente, aes(x = grupo, y = GapR, fill = grupo)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 3) +
  scale_fill_manual(
    name = "Sub-grupo",
    values = c("outlier" = "red", "normal" = "lightblue"),
    labels = c("outlier" = "Atípicos_psoraleno", "normal" = "Normales_psoraleno")
  ) +
  labs(
    title = "GapR",
    x = "Categoría",
    y = "Valores de GapR"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # título centrado y más grande
    axis.title.y = element_text(size = 25),
    axis.title.x = element_text(size= 25),
    legend.text = element_text(size = 14),          # texto de la leyenda más grande
    legend.title = element_text(size = 16, face = "bold")  # título de la leyenda también
)
