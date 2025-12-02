## **📌 Load Required Libraries**
```{r load_libraries, echo=TRUE, message=FALSE}
library(ggplot2)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library(clusterProfiler)
```

## **📌 Read Log2CPM Data**
```{r load_File Paths, echo=TRUE,results='hide',message=FALSE}
# Load feature count matrix
boxplot1 <- read.csv("data/filcpm_colnames_matrix.csv") %>% as.data.frame()

# Ensure column names are cleaned
#no need to run this for me
# colnames(boxplot1) <- trimws(gsub("^X", "", colnames(boxplot1)))
```

#now pull out some genes of interest from my data
p53_target_genes <- c("CDKN1A")


## **📌Generate Boxplots for Cardiac Genes**
```{r Cardiac, echo=TRUE,message=FALSE}
for (gene in cardiac_genes) {
  data_info <- process_gene_data(gene)
  p <- ggplot(data_info$long_data, aes(x = Condition, y = log2CPM, fill = Drug)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c("CX" = "#0000FF", "DOX" = "#e6d800", "VEH" = "#FF00FF")) +
    geom_point(aes(color = Indv), size = 2, alpha = 0.5, position = position_jitter(width = -1, height = 0)) +
    geom_text(data = data_info$significance_labels, aes(x = Condition, y = max_log2CPM + 0.5, label = Significance),
              inherit.aes = FALSE, size = 6, color = "black") +
    ggtitle(paste("Log2CPM Expression of", gene)) +
    labs(x = "Treatment", y = "log2CPM") +
    theme_bw() +
    theme(plot.title = element_text(size = rel(2), hjust = 0.5),
          axis.title = element_text(size = 15, color = "black"),
          axis.text.x = element_text(size = 10, color = "black", angle = 90, hjust = 1))

  print(p)
}
```

#########


# ----------------- Load Required Libraries -----------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
# library(clusterProfiler)
# ----------------- Load and Clean Log2CPM Matrix -----------------
boxplot1 <- read.csv("data/filcpm_colnames_matrix.csv") %>%
  as.data.frame()
# Clean column names (remove leading "X", trim whitespace)
# colnames(boxplot1) <- trimws(gsub("^X", "", colnames(boxplot1)))
# ----------------- Define Gene List -----------------
# Replace with your actual genes of interest
p53_target_genes <- c("CDKN1A", "MDM2", "BAX", "RARG", "TP53")  # Add more gene symbols as needed
# ----------------- Define Function to Process Log2CPM Data -----------------
process_gene_data <- function(gene) {
  gene_data <- boxplot1 %>% filter(SYMBOL == gene)
  long_data <- gene_data %>%
    pivot_longer(cols = -c(Entrez_ID, SYMBOL), names_to = "Sample", values_to = "log2CPM") %>%
    mutate(
      Drug = case_when(
        grepl("DOX", Sample) ~ "DOX",
        grepl("DMSO", Sample) ~ "DMSO",
        TRUE ~ NA_character_
      ),
      Timepoint = case_when(
        grepl("_24_", Sample) ~ "24",
        grepl("_24rec_", Sample) ~ "24R",
        grepl("_144rec_", Sample) ~ "144R",
        TRUE ~ NA_character_
      ),
      Indv = case_when(
        grepl("Ind1$", Sample) ~ "1",
        grepl("Ind2$", Sample) ~ "2",
        grepl("Ind3$", Sample) ~ "3",
        grepl("Ind4$", Sample) ~ "4",
        grepl("Ind5$", Sample) ~ "5",
        grepl("Ind6$", Sample) ~ "6",
        grepl("Ind6REP$", Sample) ~ "7",
        TRUE ~ NA_character_
      ),
      Condition = paste(Drug, Timepoint, sep = "_")
    )
  long_data$Condition <- factor(
    long_data$Condition,
    levels = c(
      "DOX_24", "DMSO_24", "DOX_24R", "DMSO_24R", "DOX_144R", "DMSO_144R"
    )
  )
  return(long_data)
}
# ----------------- Generate Boxplots -----------------
for (gene in p53_target_genes) {
  gene_data <- process_gene_data(gene)
  p <- ggplot(gene_data, aes(x = Condition, y = log2CPM, fill = Drug)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(color = Indv), size = 2, alpha = 0.5, position = position_jitter(width = -1, height = 0)) +
    scale_fill_manual(values = c("DOX" = "#499FBD", "DMSO" = "#BBBBBC")) +
    ggtitle(paste("Log2CPM Expression of", gene)) +
    labs(x = "Treatment", y = "log2CPM") +
    theme_bw() +
    theme(
      plot.title = element_text(size = rel(2), hjust = 0.5),
      axis.title = element_text(size = 15, color = "black"),
      axis.text.x = element_text(size = 10, color = "black", angle = 90, hjust = 1)
    )
  print(p)
}
