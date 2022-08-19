################# BASE IMAGE #####################
#FROM ubuntu:bionic-20200921
FROM ubuntu:kinetic-20220801
#FROM ubuntu:20.04
##site to test docker configuration files
# https://labs.play-with-docker.com/
################## METADATA #######################
LABEL base_image="ubuntu:bionic-20200921"
LABEL version="4.7.12"
LABEL software="BrumiR-Docker"
LABEL software.version="2.0"
LABEL about.summary="Container image containing all requirements for BrumiR toolkit"
LABEL about.home="https://github.com/camoragaq/"
LABEL about.documentation="https://github.com/camoragaq/BrumiR-Docker/README.md"
LABEL about.license_file="https://github.com/camoragaq/BrumiR-Docker/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER Carol Moraga <camoragaq@gmail.com>
################## INSTALLATION ######################
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
   wget \
   ca-certificates \
   perl-doc build-essential libomp-dev r-base
RUN apt-get install -y git
RUN cd /opt && git clone https://github.com/camoragaq/BrumiR.git
RUN cd /opt/BrumiR && mkdir bin && make all 
ENV PATH /opt/BrumiR/bin:$PATH
ENV BRUMIR /opt/BrumiR
RUN apt-get install -y r-cran-randomforest r-cran-argparse
RUN apt-get install -y vienna-rna zlib1g-dev cmake
RUN cd /opt && git clone --recursive https://github.com/GATB/bcalm
RUN cd /opt/bcalm && mkdir build &&  cd build &&  cmake .. &&  make -j 8 && cp /opt/bcalm/build/bcalm /opt/BrumiR/bin/
RUN cp /usr/bin/RNAfold /opt/BrumiR/bin
RUN cd /opt && git clone https://github.com/camoragaq/miRsim.git
RUN cd /opt/miRsim && make && mkdir bin && mv mirsim ./bin
ENV PATH /opt/miRsim/bin:$PATH
ENV PERL5LIB $BRUMIR/brumir-rf/:$PERL5LIB

# r-base-core r-cran-class r-cran-rngtools r-cran-caret \
#   r-cran-truncdist r-cran-randomfields r-cran-rcppgsl r-cran-markdown r-cran-clustergeneration \
#   r-cran-rlumshiny r-base-dev r-cran-ddalpha
