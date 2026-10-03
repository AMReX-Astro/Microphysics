#!/usr/bin/env bash

set -eu -o pipefail

CVER=${1:-4.8}

wget https://github.com/ccache/ccache/releases/download/v${CVER}/ccache-${CVER}-linux-x86_64.tar.xz
tar xvf ccache-${CVER}-linux-x86_64.tar.xz
sudo cp -f ccache-${CVER}-linux-x86_64/ccache /usr/local/bin/
