removeMapOutliers2.py -
Removes the outlying points from a mapping between phsyical location (bp) and genetic location (cM)
Input - .map file, where each line is of the following format:
[Chromosome number]	[Position name]	[cM distance from beginning of chromosome]	[bp distance from beginning of chromosome]
For example:
1	s15529.1	99.0289	67019569

findGeneticDistanceNoOutliers.py -
Calculates the genetic distance between two points in the genome
Input - .map file mapping physical distance to genetic location (with outliers removed), format described above
Output - a function interpolating, for a given physical location, it's approximate genetic location

newFilesPairGenomeV4withGeneticDistance4B2.py -
Calculates similarity scores between pairs of samples
Input -
1)	TSV files holding variations in the genome, of the following format:
	[NCBI accession]	[Location]	[Reference]	[Alternative allele]	[Alternative allele count]	[Total coverage]
	For example:
	NC_019458.2     94604   C       C       0       3
	
	and the following header:
	CHR     POS     REF     ALT     AF      COV
2)	Corresponding VCF files holding the coverage at each point in the genome, generated from each BAM file
Output - CSV files with the similarity scores between each pair of samples, for different similarity metrics