wget https://codebook.ccbr.utoronto.ca/Data/Codebook/GHT_SELEX_Metadata_web_2024_09_21.xlsx -O ghtselex_metadata.xlsx
wget https://codebook.ccbr.utoronto.ca/Data/Codebook/ChIP_Metadata_web_2024_09_25.xlsx -O chipseq_metadata.xlsl

wget https://codebook.ccbr.utoronto.ca/Data/Codebook/ChIPSeq/Merged_Peaks_narrowPeaks.tar.gz -O chipseq_raw.tar.gz
tar -xzf chipseq_raw.tar.gz
mv Merged_Peaks_narrowPeaks -O chipseq_raw

wget https://codebook.ccbr.utoronto.ca/Temp/Peaks_MAGIX_McGill.tar.gz -O ghtselex_raw.tar.gz
tar -xzf ghtselex_raw.tar.gz
mv Peaks_MAGIX_McGill -O ghtselex_raw
