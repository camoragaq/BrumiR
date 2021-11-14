################# BASE IMAGE #####################
FROM ubuntu:bionic-20200921
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
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
   wget \
   ca-certificates \
   perl-doc build-essential libomp-8-dev 
RUN apt-get install -y git
RUN cd /opt && git clone https://github.com/camoragaq/BrumiR.git
RUN cd /opt/BrumiR && mkdir bin && make all 
ENV PATH /opt/BrumiR/bin:$PATH
ENV BRUMIR /opt/BrumiR
