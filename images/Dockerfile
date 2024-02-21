
# Use Ubuntu as the base image
FROM ubuntu:focal
# Fix timezone issue
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
# Update the system and install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libboost-all-dev \
    libnlopt-dev \
    libcgal-dev \
    libomp-dev \
    libarmadillo-dev \
    libboost-math-dev \
    libboost-program-options-dev \
    libboost-test-dev \
    libboost-serialization-dev \
    git \
    wget \
    cmake \
    && rm -rf /var/lib/apt/lists/*


RUN /bin/bash -c 'git clone https://github.com/stevengj/nlopt.git && cd nlopt && mkdir build && cd build && cmake .. && make && make install'
RUN git clone https://github.com/rmcantin/bayesopt.git \
    && cd bayesopt \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install

RUN apt-get update && apt-get install -y libensmallen-dev

RUN wget https://github.com/USCiLab/cereal/archive/refs/tags/v1.1.2.tar.gz \
    && tar -xzf v1.1.2.tar.gz \
    && mkdir -p /usr/local/include/cereal \
    && cp -r cereal-1.1.2/include/cereal/* /usr/local/include/cereal/
    
#RUN  git clone  https://github.com/mlpack/mlpack.git \
#   && cd mlpack \
#   && git checkout tags/3.2.2  \
#   && mkdir build \
#   && cd build \
#   && cmake .. \
#   &&  make install

# Set the working directory in the container
WORKDIR /usr/src/bo_robotic_palpation

RUN apt-get update && apt-get install -y \
    libglu1-mesa-dev \
    freeglut3-dev \
    mesa-common-dev \
    pkg-config \
    libstb-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y \
    libmlpack-dev \
    && apt-cache policy libmlpack-dev
# Copy your project files into the container
COPY . .


# Compile your project
# Replace this with the actual build commands for your project
RUN cd build \
    && cmake .. \
    && make


# Setzen Sie den Befehl, der beim Starten des Containers ausgeführt wird
CMD ["bash"]