#!/bin/bash

#Tool installer.
# -tools a,f,s,j default all


# Default variable values
a=false
s=false
f=false
j=false

# Check if any flags are provided
if [[ $# -eq 0 ]]; then
  # No flags provided, set all variables to true
  a=true
  s=true
  f=true
  j=true
else
  # Flags are provided
  while [[ $# -gt 0 ]]; do
    flag="$1"
    case "$flag" in
      -tools)
        shift
        # Loop through each character in the tools flag argument
        for char in $(echo "$1" | fold -w 1); do
          case "$char" in
            a)
              a=true
              ;;
            s)
              s=true
              ;;
            f)
              f=true
              ;;
            j)
              j=true
              ;;
          esac
        done
        ;;
    esac
    shift
  done
fi

#INSTALL DEPENDENCIES BEFORE MAKE
Rscript ../bin/Package_install.R

if $a ; then
    rm -r ../Arriba_INS/*
    rm -r -y ../Arriba_INS/.git*
    git clone https://github.com/suhrig/arriba ../Arriba_INS
    cd ../Arriba_INS
    make
fi

if $s; then
    rm -r ../Starfusion_INS/*
    rm -r ../Starfusion_INS/.git* 
    sudo apt install cpanminus
    sudo cpanm DB_File
    sudo cpanm URI::Escape
    sudo cpanm Set::IntervalTree
    sudo cpanm Carp::Assert
    sudo cpanm JSON::XS
    sudo cpanm PerlIO::gzip
    git clone --recursive https://github.com/STAR-Fusion/STAR-Fusion ../Starfusion_INS
    cd ../Starfusion_INS
    make
fi

if $f; then
    echo -e "Prior to installation you will need the following libraries instaled:\n wget\n gawk\n gcc\n g++\n make\n cmake\n automake\n curl\n unzip\n zip\n bzip2\n tar\n gzip\n pigz\n parallel\n build-essential\n libncurses5-dev\n libc6-dev\n zlib1g\n zlib1g-dev\n libtbb-dev\n libtbb2\n python\n python-dev\n python-numpy\n python-biopython\n python-xlrd\n python-openpyxl\n default-jdk\n"
    echo -e "Type Password to install:\n"
    sudo apt-get install wget gawk gcc g++ make cmake automake curl unzip zip bzip2 tar gzip pigz parallel build-essential libncurses5-dev libc6-dev zlib1g zlib1g-dev libtbb-dev libtbb2 python python-dev python-numpy python-biopython python-xlrd python-openpyxl default-jdk
    rm -r ../FusionCatcher_INS/*
    rm -r ../FusionCatcher_INS/.git* 
    git clone https://github.com/ndaniel/fusioncatcher ../FusionCatcher_INS
    cd ../FusionCatcher_INS/tools/
    ./install_tools.sh
    cd ..
fi

if $j; then
    rm -r ../JAFFA_INS/*
    rm -r ../JAFFA_INS/.git* 
    git clone https://github.com/Oshlack/JAFFA ../JAFFA_INS
    cd ../JAFFA_INS
    ./install_linux64.sh
fi


os=$(uname -s)

# Install the JDK based on the operating system
if [[ "$os" == "Linux" ]]; then
  # Linux OS
  echo "Detected Linux OS."
  echo "Installing OpenJDK..."

  # Download and extract the JDK tarball
  wget -O ../jdk.tar.gz https://download.oracle.com/java/20/latest/jdk-20_linux-x64_bin.tar.gz
  tar -xzf ../jdk.tar.gz 

#   # Set the environment variables
#   export JAVA_HOME=/usr/lib/jvm/jdk*/   # Update the path according to your JDK version
#   export PATH=$JAVA_HOME/bin:$PATH

elif [[ "$os" == "Darwin" ]]; then
  # macOS
  echo "Detected macOS."
  echo "Installing AdoptOpenJDK..."

  # Install AdoptOpenJDK using Homebrew
  brew update
  brew tap AdoptOpenJDK/openjdk
  brew install --cask adoptopenjdk

else
  # Unsupported operating system
  echo "Unsupported operating system: $os. JDK installation cannot proceed."
  exit 1

fi

# Verify JDK installation
java_version=$(java -version 2>&1 | awk -F '"' '/version/ {print $2}')
echo "JDK installation completed. Java version: $java_version"
