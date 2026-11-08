#!/bin/bash
set -e

LIGHT_PURPLE='\033[1;35m'
ORIGINAL='\033[0m'

P=${LIGHT_PURPLE}
O=${ORIGINAL}


# install dependencies
echo -e "\n${P}Installing dependencies...${O}\n"
sudo apt-get update
sudo apt-get install -y samtools zlib1g-dev make gcc
echo -e "\n${P}done!${O}\n"
