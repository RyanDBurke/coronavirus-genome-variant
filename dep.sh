#!/bin/bash

LIGHT_PURPLE='\033[1;35m'
ORIGINAL='\033[0m'

P=${LIGHT_PURPLE}
O=${ORIGINAL}


# install dependencies
echo -e "\n${P}Installing dependencies...${O}\n"
sudo apt-get update 
sudo apt-get install -y samtools 
sudo apt-get install libz-dev 
sudo apt install make 
sudo apt install gcc 
echo -e "\n${P}done!${O}\n"