####################################################
# 
### BUILD IMAGE
# docker build -t icrsc/locus-explorer:latest .
#
### RUN IMAGE
# docker run -p 3838:3838 --rm --name locus-explorer -v  /mnt/c/Users/ralcraft/Documents/shiny-proxy-data/locus-explorer:/app/mnt icrsc/locus-explorer:latest
# docker run -p 3838:3838 --rm --name locus-explorer icrsc/locus-explorer:latest
#
# Should then be visible on http://localhost:3838/
#
# PUSH IMAGE
# docker push icrsc/locus-explorer:latest
# 
# docker image prune
####################################################

#https://github.com/lescailab/r-ggbio-reporting/blob/main/Dockerfile

#FROM r-base
FROM rocker/r-ver:4.3.2

# IGRAPH - https://r.igraph.org/articles/installation-troubleshooting
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \ 
    gfortran \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libmkl-rt \
    libglpk-dev \
    libxml2-dev \
    apt-file \    
    && rm -rf /var/lib/apt/lists/*


RUN apt-get update
RUN apt-get install -y \
procps \
libcurl4-openssl-dev \
libfontconfig1-dev \
libxml2 \
libxml2-dev \
libz-dev \
libbz2-dev \
libclang-dev \
liblzma-dev \
python3 \
python3-pip \
wget \
make \
gcc \
libssl-dev \
pandoc

WORKDIR /app
# Copy relevant app files
RUN mkdir www
COPY www ./www
RUN mkdir Data
COPY Data ./Data
RUN mkdir Markdown
COPY Markdown ./Markdown
RUN mkdir mnt
COPY global.R server.R ui.R README.md ./

RUN Rscript -e "install.packages('BiocManager')"

RUN Rscript -e "BiocManager::install(c(\
'tidyverse', \
'stringr', \
'biovizBase', \
'GenomicRanges', \
'ggbio', \
'knitr' \
))"

RUN Rscript -e 'install.packages(c("shiny","dplyr","tidyr","lazyeval","data.table","ggplot2","ggrepel","knitr","markdown","DT","lattice","acepack","cluster","DBI","colourpicker","igraph","visNetwork", "devtools"))'
RUN R -e 'BiocManager::install(c("TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db","rtracklayer"))'
RUN installGithub.r oncogenetics/oncofunco
              
EXPOSE 3838
CMD Rscript -e "shiny::runApp('/app', port = 3838, host = '0.0.0.0')"