input=$1
ctrl=$2
output=${input}_rmBL.bdg
#output=$3
bamCompare -b1 bam/${input}_rDNA.rmdup.rmM.MAPQ30.sorted.bam -b2 bam/${ctrl}_rDNA.rmdup.rmM.MAPQ30.sorted.bam \
        --scaleFactorsMethod None \
        --skipZeroOverZero \
        --normalizeUsing RPKM \
        --effectiveGenomeSize 2864785220 \
        -bl blacklist/bl.bed \
	--binSize 10000 --smoothLength 30000 \
        -p 20 \
        -of bedgraph -o bedgraph/${output}
