FROM continuumio/miniconda3

ENV HOME /root
RUN apt-get install -y curl build-essential
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda create -n env python=3.7 htseq numpy pandas
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
WORKDIR ${HOME}
RUN git clone https://github.com/emrobe/peruse.git
RUN git clone https://github.com/marbl/Krona.git
RUN mkdir ${HOME}/taxonomy
RUN ${HOME}/Krona/KronaTools/updateTaxonomy.sh ${HOME}/taxonomy
