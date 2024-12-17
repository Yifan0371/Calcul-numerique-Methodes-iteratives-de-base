# Use the official Ubuntu 20.04 as the base image
FROM ubuntu:20.04

# Set non-interactive mode for APT to prevent prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary dependencies: GCC compiler, BLAS, LAPACK, LAPACKE libraries, Make, Vim, and Gnuplot
# Clean up APT cache to reduce the image size
RUN apt-get update -y && \
    apt-get install -y gcc libblas-dev liblapack-dev liblapacke-dev make vim gnuplot && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory to /app where the application code will reside
WORKDIR /app

# Copy all the local project files into the container's /app directory
COPY . .

# Set HOSTNAME environment variable to "docker" for consistent execution behavior
ENV HOSTNAME=docker

# Create the bin directory and build all targets using the Makefile
RUN mkdir -p bin && \
    make all

# Define the default command to run when the container starts
# This will execute the "make run_tpPoisson1D_direct" target from the Makefile
CMD ["bash", "-c", "make run_testenv && make run_tpPoisson1D_direct"]

