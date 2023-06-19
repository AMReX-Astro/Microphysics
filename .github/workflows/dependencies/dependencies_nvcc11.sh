#!/usr/bin/env bash
#
# Copyright 2020 Axel Huebl
#
# License: BSD-3-Clause-LBNL

set -eu -o pipefail

sudo apt-get -qqq update
sudo apt-get install -y \
    build-essential     \
    ca-certificates     \
    cmake               \
    g++                 \
    gnupg               \
    libopenmpi-dev      \
    openmpi-bin         \
    pkg-config          \
    wget

sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/3bf863cc.pub
echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64 /" \
    | sudo tee /etc/apt/sources.list.d/cuda.list
sudo apt-get update
sudo apt-get install -y \
    cuda-command-line-tools-11-7 \
    cuda-compiler-11-7           \
    cuda-cupti-dev-11-7          \
    cuda-minimal-build-11-7      \
    cuda-nvml-dev-11-7           \
    cuda-nvtx-11-7               \
    libcurand-dev-11-7
sudo ln -s cuda-11.7 /usr/local/cuda

