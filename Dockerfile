# Use an official Ubuntu base image
FROM ubuntu:latest

# Update and install necessary packages
RUN apt-get update && apt-get install -y \
    perl \
    libdbi-perl \
    bioperl \
    libwww-perl \
    python3 \
    git \
    ftp \
    wget \
    vim \
    bash \
    libswitch-perl

# Install GCC 10.3.0 for sift4g
RUN apt-get install -y \
	gcc-10 \
	g++-10 \
	make

# create an alias for python3
RUN ln -s /usr/bin/python3 /usr/bin/python

# install java for SIFT4G Annotator 
RUN apt-get update && apt-get install -y openjdk-11-jdk

# Create symlinks to make GCC 10.3.0 the default
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100 \
    --slave /usr/bin/g++ g++ /usr/bin/g++-10

# Clone the sift4g repository from GitHub
WORKDIR /
RUN git clone --recursive https://github.com/rvaser/sift4g.git

# compile sift4g for CPU
WORKDIR /sift4g
RUN make

# copy code from this repo for creating the database

RUN mkdir -p /SIFT4G_Create_Genomic_DB
COPY . /SIFT4G_Create_Genomic_DB
RUN chmod 775 /SIFT4G_Create_Genomic_DB/common-utils.pl

# Set environment variables 
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
ENV PATH=$PATH:$JAVA_HOME/bin

RUN wget https://github.com/pauline-ng/SIFT4G_Annotator/raw/master/SIFT4G_Annotator.jar
RUN chmod 775 SIFT4G_Annotator.jar

# You can add any additional commands here to build or run your application
WORKDIR /SIFT4G_Create_Genomic_DB
