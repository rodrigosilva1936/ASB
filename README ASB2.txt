METADATA
Obtivemos o .csv com a metadata através do NCBI

SILVA (v138_1) #
mkdir ~/databases
cd ~/databases
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138_1.tgz 
tar xfv silva.seed_v138_1.tgz
rm silva.seed_v138_1.tgz
#Script para ir buscar o ficheiro silva.seed_v138_1. da base de dados SILVA amplamente utilizada para sequências de rRNA 16S em táxons conhecidos 



FASTQ (sratoolkit.3.2.1-ubuntu64)
set -e   # Se ocorrer algum erro o script irá sair imediatamente                                                              
sed '1d' SraRunTable2.csv | while read meta  
do
  accn=$(echo $meta | cut -d "," -f 1)
  name=$(echo $meta | cut -d "," -f 39)
  fasterq-dump ${accn}
  head -n 200000 ${accn}_1.fastq | gzip > ${name}_1.fastq.gz
  rm ${accn}_1.fastq
  head -n 200000 ${accn}_2.fastq | gzip > ${name}_2.fastq.gz
  rm ${accn}_2.fastq
done
#Este bloco automatiza o download de dados brutos de sequenciação do NCBI SRA, converte-os para o formato FASTQ, parte as reads para um número fixo (200000 reads por amostra) e compacta-as, eliminando os ficheiros temporários maiores.


PREPARAÇÃO DE AMBIENTE PARA PIPELINE E CONFIGURAR NEXTFLOW
micromamba create -n nextflow_env -c bioconda -c conda-forge nextflow
micromamba activate nextflow_env
cd ~/                           
git clone https://github.com/hawaiidatascience/metaflowmics
cd ~/metaflowmics/metaflowmics  
nano Pipeline-16S/nextflow.config

# Este bloco configura um ambiente Micromamba para o Nextflow, baixa o pipeline metaflowmics e prepara o ficheiro de configuração do pipeline para execução.




NEXTFLOW PIPELINE (Docker version 27.5.1, build 27.5.1-0ubuntu3~24.04.2, Java version openjdk 21.0.7 2025-04-15, nextflow version 25.04.2.5947)

nextflow run Pipeline-16S -profile docker \
  --reads "/home/rodrigo/metagenomics/*_{1,2}.fastq" \
  --referenceAln /home/rodrigo/databases/silva.seed_v138_1.align \
  --referenceTax /home/rodrigo/databases/silva.seed_v138_1.tax \
  --outdir ~/drought_results \
  --cpus 6 \
  --single_end true \
  -resume \
  -process.memory '12 GB'

#Este comando executa o pipeline de bioinformática metaflowmics usando Docker, processando os ficheiros FASTQ contra as bases de dados SILVA para análise de sequências 16S, guardando os resultados para uma pasta e otimizando o uso de recursos.




PREPARAÇÃO FICHEIROS ALPHA E BETA
cut -f 2 alpha-diversity.97.summary | cut -d "_" -f 1 |sed 's/group/type/' | paste alpha-diversity.97.summary - > alpha-div.tsv
sed 's/\t\t/\t/' beta-diversity.97.summary | sed 's/comparison/comparison1\tcomparison2/' > beta-div.tsv

#Este bloco formata e limpa os ficheiros de resumo de diversidade alfa e beta gerados pelo pipeline, ajustando os cabeçalhos e a estrutura para facilitar a importação e manipulação no ambiente R.




VISUALIZAÇÃO GRÁFICOS INICIAIS
firefox *.html   #Abre todos os ficheiros html gerados pelo pipeline


SCRIPT R (R version 4.5.0 2025-04-11, RStudio 2025.05.1 Build 513)

#!/usr/bin/Rscript

beta = read.csv("/home/rodrigo/drought_results/postprocessing/beta-div.tsv", sep="\t")
beta$comparison1 <- as.character(beta$comparison1)
beta$comparison2 <- as.character(beta$comparison2)
bc_mat <- beta[,c(2,3,4)]

# Get unique sample names
samples <- unique(c(bc_mat$comparison1, bc_mat$comparison2))

# Create an empty matrix with NAs
triangular_matrix <- matrix(NA, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))

# Fill in the matrix with the distance values
for (i in 1:nrow(bc_mat)) {
  row_name <- bc_mat$comparison1[i]
  col_name <- bc_mat$comparison2[i]
  value <- bc_mat$braycurtis[i]
  triangular_matrix[row_name, col_name] <- value
}

# Replace missing values (diagonal and lower triangle) with 0
triangular_matrix[is.na(triangular_matrix) | lower.tri(triangular_matrix)] <- 0
#print(triangular_matrix)

pcoa <- as.data.frame(cmdscale(triangular_matrix, k = 2))
names(pcoa) <- c("PCoA1", "PCoA2")

my_labels = gsub("[_].*", "", rownames(pcoa))
plot(pcoa$PCoA1, pcoa$PCoA2, col=as.factor(my_labels), xlab="PCoA1", ylab="PCoA2")
legend("top", unique(my_labels), pch=1, col=(unique(as.factor(my_labels))))
#text(pcoa$PCoA1, pcoa$PCoA2, samples, pos=3)

#Este script R lê uma tabela de dissimilaridade beta (Bray-Curtis) entre amostras, constrói uma matriz de distância completa e realiza uma Análise de Coordenadas Principais (PCoA) em 2 dimensões. Ele então gera um gráfico de dispersão simples da PCoA, onde cada ponto representa uma amostra e as cores são atribuídas com base nos prefixos dos nomes das amostras, com uma legenda correspondente.




library(vegan)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ape)

beta_data <- read.delim("/home/rodrigo/drought_results/postprocessing/beta-div.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

samples <- unique(c(beta_data$comparison1, beta_data$comparison2))
dist_matrix <- matrix(NA, nrow = length(samples), ncol = length(samples),
                      dimnames = list(samples, samples))

for (i in 1:nrow(beta_data)) {
    s1 <- beta_data$comparison1[i]
    s2 <- beta_data$comparison2[i]
    bc <- beta_data$braycurtis[i]
    dist_matrix[s1, s2] <- bc
    dist_matrix[s2, s1] <- bc
}
diag(dist_matrix) <- 0
dist_object <- as.dist(dist_matrix)

metadata_sra <- read.csv("/home/rodrigo/Downloads/SraRunTable2.csv", header = TRUE, stringsAsFactors = FALSE)
metadata_sra$Sample <- trimws(metadata_sra$Run)

metadata_lake <- read.csv("/home/rodrigo/Downloads/amostras_lago.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(metadata_lake)[colnames(metadata_lake) == "Lake"] <- "Lake_of_Origin"

combined_metadata <- merge(metadata_sra, metadata_lake, by = "Sample", all.x = TRUE)

combined_metadata$Host_diet_Clean <- trimws(as.character(combined_metadata$Host_diet))
combined_metadata$Lake_of_Origin_Clean <- trimws(as.character(combined_metadata$Lake_of_Origin))

combined_metadata$Grupo_Final <- paste0(combined_metadata$Host_diet_Clean, "_", combined_metadata$Lake_of_Origin_Clean)

samples_in_dist_matrix <- rownames(as.matrix(dist_object))

filtered_metadata_present_in_dist <- combined_metadata %>%
  filter(Sample %in% samples_in_dist_matrix)

final_merged_data <- filtered_metadata_present_in_dist %>%
  filter(Host_diet_Clean %in% c("benthic", "limnetic")) %>%
  filter(!is.na(Grupo_Final) & Grupo_Final != "")

samples_for_nmds <- final_merged_data %>% pull(Sample)

dist_subset <- as.dist(as.matrix(dist_object)[samples_for_nmds, samples_for_nmds])

desired_group_order <- c(
  "benthic_Paxton", "benthic_Priest", "benthic_Little Quarry",
  "limnetic_Paxton", "limnetic_Priest", "limnetic_Little Quarry"
)
actual_groups_present <- intersect(desired_group_order, unique(final_merged_data$Grupo_Final))
final_merged_data$Grupo_Final <- factor(final_merged_data$Grupo_Final, levels = actual_groups_present)

nmds_result <- metaMDS(dist_subset, autotransform = FALSE, wascores = FALSE, k = 2, trymax = 100)

print(paste("NMDS Stress:", nmds_result$stress))
plot(nmds_result, type = "t")

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$Sample <- rownames(nmds_scores)

final_nmds_data_for_plot <- merge(nmds_scores, final_merged_data[, c("Sample", "Grupo_Final", "Host_diet_Clean")], by = "Sample", all.x = TRUE)

final_nmds_data_for_plot <- final_nmds_data_for_plot %>%
    filter(!is.na(Grupo_Final))

hulls <- final_nmds_data_for_plot %>%
    group_by(Grupo_Final) %>%
    slice(chull(NMDS1, NMDS2))

num_unique_groups <- length(unique(final_nmds_data_for_plot$Grupo_Final))

colors_for_groups <- brewer.pal(min(num_unique_groups, 8), "Dark2")
if (num_unique_groups > 8) {
    colors_for_groups <- c(brewer.pal(8, "Dark2"), brewer.pal(min(num_unique_groups - 8, 8), "Set2"))
}

shapes_for_groups <- c(16, 17, 15, 18, 8, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12)[1:num_unique_groups]

ggplot(final_nmds_data_for_plot, aes(x = NMDS1, y = NMDS2, color = Grupo_Final, fill = Grupo_Final)) +
    geom_point(size = 3, aes(shape = Grupo_Final)) +
    geom_polygon(data = hulls, alpha = 0.3, show.legend = FALSE) +
    scale_color_manual(values = colors_for_groups) +
    scale_shape_manual(values = shapes_for_groups) +
    labs(title = "NMDS da Composição do Microbioma por Dieta e Lago (Benthic e Limnetic)",
         x = "NMDS1",
         y = "NMDS2",
         color = "Grupo (Dieta_Lago)",
         fill = "Grupo (Dieta_Lago)") +
    theme_minimal() +
    theme(legend.position = "right")

 
O script lê dados de dissimilaridade beta e metadados de amostras para comunidades microbianas. Ele constrói uma matriz de distância, filtra as amostras para incluir apenas aquelas com dados completos de dieta ("benthic"/"limnetic") e lago, e então executa uma análise NMDS. Os resultados da NMDS são combinados com os metadados para criar um gráfico no ggplot2 que visualiza o agrupamento das amostras por uma combinação de dieta e lago, usando pontos coloridos e formas, e desenhando contornos ao redor dos grupos.

GITHUB
git init
git remote add origin git@github.com:rodrigosilva1936/ASB.git
git add .
git commit -m "Adiciona scripts iniciais e dados"
git branch -M main
git push -u origin main
hash commit : f6fbb498c7df0893a24171e7510820e3526d6227

SOFTWARES USADOS:
SRA Toolkit (especificamente fasterq-dump versão: sratoolkit.3.2.1-ubuntu64)
Nextflow: (nextflow version 25.04.2.5947)
Docker: (Docker version 27.5.1, build 27.5.1-0ubuntu3~24.04.2)
Java: (Java version openjdk 21.0.7 2025-04-15)
R: (R version 4.5.0 2025-04-11) 
RStudio: (RStudio 2025.05.1 Build 513)
Firefox 138.0.4(64-bit)
Micromamba 2.0.8

