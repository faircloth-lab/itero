FROM centos:7
MAINTAINER "Brant Faircloth" <brant _at_ faircloth-lab _dot_ org>
ENV container docker

# update yum
RUN yum -y update; yum clean all

# install some packages
RUN yum -y install wget tar bzip2 git

# add test user
RUN useradd -ms /bin/bash test

# switch to test user
USER test
ENV HOME /home/test
WORKDIR /home/test/

# download conda
RUN ["/bin/bash", "-c", "wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O $HOME/miniconda.sh"]
RUN chmod 0755 $HOME/miniconda.sh
RUN ["/bin/bash", "-c", "$HOME/miniconda.sh -b -p $HOME/conda"]
ENV PATH="$HOME/conda/bin:$PATH"
RUN rm $HOME/miniconda.sh

# update conda
RUN conda update conda

# linkup bioconda
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# install ipython
RUN conda install ipython numpy biopython argcomplete schwimmbad bedtools bwa gawk samtools spades

# clone phyluce source to $HOME/git/phyluce
RUN mkdir git
RUN cd $HOME/git && git clone https://github.com/faircloth-lab/itero.git
RUN cd $HOME/git/itero && git fetch --all && git checkout -b working origin/working

ENV PYTHONPATH="/home/test/git/itero:$PYTHONPATH"