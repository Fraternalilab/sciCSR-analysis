FROM ubuntu:xenial
MAINTAINER sridhar <sridhar@wustl.edu>

LABEL description="Image for STAR aligner version 2.5.1.b, HISAT2 v2.2.1, Salmon v1.9.0 and subreads v2.0.3"

RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    cmake \
    default-jdk \
    git \
    libnss-sss \
    libtbb2 \
    libtbb-dev \
    ncurses-dev \
    python-dev \
    python-pip \
    tzdata \
    unzip \
    wget \
    zlib1g \
    zlib1g-dev \
    python3 \
    samtools \
    vim

ENV star_version 2.5.1b
WORKDIR /docker_main
RUN wget -O /docker_main/hista2-2.2.1-Linux_x86_64.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
RUN wget -O /docker_main/subreads-2.0.3-Linux-x86_64.tar.gz https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz/download
# Install STAR
ADD https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz /usr/bin/
RUN tar -xzf /usr/bin/${star_version}.tar.gz -C /usr/bin/
RUN cp /usr/bin/STAR-${star_version}/bin/Linux_x86_64/* /usr/local/bin

# Install HISAT2
RUN unzip /docker_main/hista2-2.2.1-Linux_x86_64.zip

# INSTALL subread
RUN cd /docker_main
RUN tar -xzvf subreads-2.0.3-Linux-x86_64.tar.gz

ENV PACKAGES git gcc make g++ libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev libcurl4-openssl-dev curl unzip autoconf apt-transport-https ca-certificates gnupg software-properties-common wget
ENV SALMON_VERSION 1.9.0

# salmon binary will be installed in /home/salmon/bin/salmon

### don't modify things below here for version updates etc.

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | apt-key add -

RUN apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'

RUN apt-get update

RUN apt-key --keyring /etc/apt/trusted.gpg del C1F34CDD40CD72DA

RUN apt-get install kitware-archive-keyring

RUN apt-get install -y libstdc++6 libssl1.1 libc6 cmake

RUN curl -k -L https://github.com/COMBINE-lab/salmon/archive/v${SALMON_VERSION}.tar.gz -o salmon-v${SALMON_VERSION}.tar.gz && \
    tar xzf salmon-v${SALMON_VERSION}.tar.gz && \
    cd salmon-${SALMON_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/salmon && make && make install

RUN apt-get update \
    && apt-get install -y --no-install-recommends libhwloc5 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=base /usr/local/salmon/ /usr/local/
ENV PATH /home/salmon-${SALMON_VERSION}/bin:${PATH}
ENV LD_LIBRARY_PATH "/usr/local/lib:${LD_LIBRARY_PATH}"

RUN echo "export PATH=$PATH" > /etc/environment
RUN echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" > /etc/environment

##lsf time stamp bug
## borrow from cmiller

RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

RUN cd /docker_main / && \
   rm -rf 2.7.0f.tar.gz && \
   apt-get autoremove -y && \
   apt-get autoclean -y  && \
   apt-get clean

# Set default working path
WORKDIR /docker_main

