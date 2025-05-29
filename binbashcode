# Download spaceranger
wget -O spaceranger.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-2.1.1.tar.gz"
tar -xzvf spaceranger.tar.gz

# Get absolute path to spaceranger folder
SPACERANGER_PATH="$PWD/spaceranger-2.1.1"

# Add spaceranger to PATH for this terminal session only
export PATH=$SPACERANGER_PATH:$PATH

# Add spaceranger to PATH permanently 
SHELL_CONFIG=""

if [ -n "$ZSH_VERSION" ]; then
  SHELL_CONFIG="$HOME/.zshrc"
elif [ -n "$BASH_VERSION" ]; then
  SHELL_CONFIG="$HOME/.bashrc"
else
  echo "Unknown shell. Please add this line manually to your shell config:"
  echo "export PATH=$SPACERANGER_PATH:\$PATH"
  exit 1
fi

grep -qxF "export PATH=$SPACERANGER_PATH:\$PATH" $SHELL_CONFIG || echo "export PATH=$SPACERANGER_PATH:\$PATH" >> $SHELL_CONFIG

echo "Added spaceranger to PATH in $SHELL_CONFIG"
echo "Please restart your terminal or run 'source $SHELL_CONFIG' to apply changes."

# Download 10x mouse reference genome
mkdir -p reference
cd reference
wget -O refdata-gex-mm10-2020-A.tar.gz "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-mm10-2020-A.tar.gz"
tar -xvzf refdata-gex-mm10-2020-A.tar.gz
cd ..

nano run_spaceranger_all.sh

#!/bin/bash - COPY INTO NANO, check paths

# config
LOCAL_CORES=14
LOCAL_MEM=30
TRANSCRIPTOME="/path/to/refdata-gex-mm10-2020-A"

# TMA 1 
echo "Running TMA1..."
spaceranger count \
--id=TMA1_output \
--transcriptome=$TRANSCRIPTOME \
--libraries=/path/to/TMA1_libraries.csv \
--sample=TMA1_sample \
--slide=TMA1_slide_id \
--area=TMA1_area_id \
--image=/path/to/TMA1_HE_image.tif \
--loupe-alignment=/path/to/TMA1_alignment.json \
--localcores=$LOCAL_CORES \
--localmem=$LOCAL_MEM \
> TMA1_output.log 2>&1

# TMA 2 
echo "Running TMA2..."
spaceranger count \
--id=TMA2_output \
--transcriptome=$TRANSCRIPTOME \
--libraries=/path/to/TMA2_libraries.csv \
--sample=TMA2_sample \
--slide=TMA2_slide_id \
--area=TMA2_area_id \
--image=/path/to/TMA2_HE_image.tif \
--loupe-alignment=/path/to/TMA2_alignment.json \
--localcores=$LOCAL_CORES \
--localmem=$LOCAL_MEM \
> TMA2_output.log 2>&1

# TMA 3 
echo "Running TMA3..."
spaceranger count \
--id=TMA3_output \
--transcriptome=$TRANSCRIPTOME \
--libraries=/path/to/TMA3_libraries.csv \
--sample=TMA3_sample \
--slide=TMA3_slide_id \
--area=TMA3_area_id \
--image=/path/to/TMA3_HE_image.tif \
--loupe-alignment=/path/to/TMA3_alignment.json \
--localcores=$LOCAL_CORES \
--localmem=$LOCAL_MEM \
> TMA3_output.log 2>&1

# TMA 4 
echo "Running TMA4..."
spaceranger count \
--id=TMA4_output \
--transcriptome=$TRANSCRIPTOME \
--libraries=/path/to/TMA4_libraries.csv \
--sample=TMA4_sample \
--slide=TMA4_slide_id \
--area=TMA4_area_id \
--image=/path/to/TMA4_HE_image.tif \
--loupe-alignment=/path/to/TMA4_alignment.json \
--localcores=$LOCAL_CORES \
--localmem=$LOCAL_MEM \
> TMA4_output.log 2>&1

echo "All runs finished"

# script executable
chmod +x run_spaceranger_all.sh

# run script 
./run_spaceranger_all.sh



