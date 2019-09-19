FROM continuumio/miniconda3:4.6.14

RUN apt-get update && apt-get install -y procps tabix bcftools
RUN conda create -n env -c bioconda freebayes==1.3.1 vcflib==1.0.0_rc3 gsort==0.0.6
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
