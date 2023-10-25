FROM rocker/rstudio:4.3.1

RUN R -e 'install.packages("devtools")'

RUN R -e 'install.packages("BiocManager")'
RUN R -e 'install.packages("remotes")'

RUN apt-get update
RUN apt install -y zlib1g-dev
RUN R -e 'BiocManager::install(c("SingleR","limma","cmapR","celldex"))'
RUN R -e 'install.packages("Seurat")'

WORKDIR /home/rstudio

COPY . .

RUN R -e 'install.packages(".", repos = NULL, type = "source")'

# WORKDIR /home/rstudio/build
# RUN mkdir -p /home/rstudio/build/DrugReference
# RUN R -e 'library("Asgard"); PrepareReference(cell.info="GSE70138_Broad_LINCS_cell_info_2017-04-28.txt", gene.info="GSE70138_Broad_LINCS_gene_info_2017-03-06.txt", GSE70138.sig.info = "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt", GSE92742.sig.info = "GSE92742_Broad_LINCS_sig_info.txt", GSE70138.gctx = "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx", GSE92742.gctx = "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", Output.Dir = "DrugReference/")'

# RUN mv DrugReference /home/rstudio/.

# WORKDIR /home/rstudio

# RUN rm -rf /home/rstudio/build
