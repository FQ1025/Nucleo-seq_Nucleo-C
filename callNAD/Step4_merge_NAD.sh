samples=`realpath output/*.bed`
if [ ! -d "processed_output" ]; then
    mkdir processed_output
fi

for i in $samples
do
	name=`basename ${i} .bed`
	./filter_nad.sh -i $i -o processed_output/${name}_merge.bed -d 20000 -l 50000
done
