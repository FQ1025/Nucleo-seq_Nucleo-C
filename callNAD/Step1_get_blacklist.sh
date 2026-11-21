if [ ! -d "blacklist" ]; then
    mkdir blacklist
fi
if [ ! -d "output" ]; then
    mkdir output
fi
if [ ! -d "bedgraph" ]; then
    mkdir bedgraph
fi

cat /lustre/user/liclab/fangq/ref/UCSC/hg19/hg19_Stalk.bed /lustre/user/liclab/fangq/ref/UCSC/hg19/hg19.centromere.bed |cut -f 1,2,3 | bedtools sort | bedtools merge -d 2000000 -i - > blacklist/bl.bed
