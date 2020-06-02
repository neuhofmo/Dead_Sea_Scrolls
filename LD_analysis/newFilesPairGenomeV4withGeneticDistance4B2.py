##########################################################################################
# Calculate similarity while taking into consideration Likage Disequilibrium

# This script accompanies the paper "Illuminating Genetic Mysteries of the Dead Sea Scrolls"
# Author: Or Sagy
##########################################################################################

# Calculate similarity scores between pairs of samples.
# The scores are based on genomic locations mapped in both samples -
# the ratio of cases where the same nucleotide value was mapped in both
# samples at the location out of all of the cases where a value different
# from the reference value was mapped in at least one of the samples at the location.
# These are also "weighted" in an attempt to eliminate bias caused by linkage disequilibrium.

import gzip
import os
import itertools
import glob
from findGeneticDistanceNoOutliers import findGeneticDistance

GENETIC_DISTANCE_CUTOFF = 5
MIN_COVERAGE = 3
nucleotides = "ACGT"

def chromNameToNum(chrom):
    """Convert an OARv4 chromosome's NCBI reference sequence to a number between 1 and 27

    NCBI reference sequences for OARv4 are:
    * NC_019458.2 for chromosome I to NC_019483.2 for chromosome XXVI (1-26 returned)
    * NC_019484.2 for chromosome X (27 returned)
	"""
    if(chrom.startswith("NC_0194") and chrom.endswith(".2")):
        return int(chrom[6:9])-457
    else:
        return 0

# TSV files holding variations in the genome
filesTSV = glob.glob("~/scrolls/cov_w_alt_col/dss*.tsv.gz")

# Corresponding VCF files holding the coverage at each point in the genome, from each BAM file
filesVCF = [file.replace("cov_w_alt_col", "vcf_collapsed").replace("cov3.af_alt_cov.tsv.gz", "_SNP.vcf") for file in filesTSV]

# A sorted list of the names of scroll samples taken from their filenames
# (plus "SRR514" for the outgroup, sequencing of a Tibetan sheep)
sampleNames = [os.path.basename(filename)[:6] for filename in filesVCF] + ["SRR514"]
sampleNames = sorted(set(sampleNames))

# snpValsPerFragment will hold, for each fragment, a mapping between each sequenced location and
# how many times each value mapped to it (a 4-tuple holding the count of A, C, G, T)
# (values taken from the TSVs and VCFs)
snpValsPerFragment = {"SRR514": dict()}

# refValsPerFragment will hold, for each fragment, a mapping between each sequenced location and
# the nucleotide at that location in the reference genome
# (values taken from the TSVs)
refValsPerFragment = {"SRR514": dict()}

vals = dict()

# Going over each fragment that was sequenced
i = 0
for (vcfFilename, tsvFilename) in zip(filesVCF, filesTSV):
    sampleName = os.path.basename(vcfFilename)[:6]
    if(sampleName != os.path.basename(tsvFilename)[:6]):
        print("Weird... vcf and perbase aren't in sync -", vcfFilename, "/", tsvFilename)
        input()
        continue

    if sampleName in snpValsPerFragment:
        snpVals = snpValsPerFragment[sampleName]
    else:
        snpVals = {}

    if sampleName in refValsPerFragment:
        refVals = refValsPerFragment[sampleName]
    else:
        refVals = {}

    print(i, "/", len(filesVCF), "-", vcfFilename)
    i += 1

    sampleID = sampleNames.index(sampleName)

	# A dictionary used to verify that each location appears only once for a certainsample
    foundInVCF = dict()

    with open(vcfFilename) as vcfFile:
        line = vcfFile.readline()
		
		# Skip the header lines
        while(line.startswith("##")):
            line = vcfFile.readline()
        if(not line.startswith("#CHROM")):
            print("Unexpected line:", line)
            input()

        for line in vcfFile:
            splitLine = line.rstrip().split()
            chrom = splitLine[0]

            chromNum = chromNameToNum(chrom)
            if(chromNum == 0):
                continue

            loc = int(splitLine[1])	# Physical location (bp) on the chromosome

            if(not splitLine[7].split(";")[1].startswith("CS=")):
                print("Weird -", splitLine[7])

			# Extract the number of reads where each nucleotide value was mapped to the location
            nucleotideCnts = [int(val) for val in splitLine[7].split(";")[1][3:].split(",")]

            if((chromNum, loc) in foundInVCF):
                print("Weird -", (chromNum, loc), "already in foundInVCF, =", foundInVCF[(chromNum, loc)])
                input()
            foundInVCF[(chromNum, loc)] = sum(nucleotideCnts)

			# Increment the count of each nucleotide at the location for the relevant fragment
			# (there can be multiple samples of the same fragment)
            if (chromNum, loc) not in snpVals:
                snpVals[(chromNum, loc)] = [0 for nuc in nucleotides]
            snpVals[(chromNum, loc)] = [sum(val) for val in zip(snpVals[(chromNum, loc)], nucleotideCnts)]

    with gzip.open(tsvFilename) as tsvFile:
        tsvFile.readline()

        for line in tsvFile:
            splitLine = line.decode().rstrip().split()
            chrom = splitLine[0]
            loc = int(splitLine[1])
            (ref, alt) = splitLine[2:4]
            af = float(splitLine[4])
            cov = int(splitLine[5])

            chromNum = chromNameToNum(chrom)
            if(chromNum == 0):
                continue

			# The cases where only the reference nucleotide was mapped at the location are
			# found only in the TSV
            if ref==alt:
                if (chromNum, loc) not in snpVals:
                    snpVals[(chromNum, loc)] = [0 for nuc in nucleotides]

                snpVals[(chromNum, loc)][nucleotides.find(ref)] += cov
                if (chromNum, loc) in refVals:
                    if refVals[(chromNum, loc)] != ref:
                        refVals[(chromNum, loc)] = "X"
                else:
                    refVals[(chromNum, loc)] = ref

	# Update the general dictionaries with the new values
    snpValsPerFragment[sampleName] = snpVals
    refValsPerFragment[sampleName] = refVals

# A set of all the locations "strictly" mapped (as explained hereinafter) in at least one fragment
allSNPlocs = set()

# For each fragment, we only keep the locations that had exactly one nucleotide mapped to them
# and it was covered more than MIN_COVERAGE times, and map the location to that nucleotide value
strictVals = dict()
for sampleName in sampleNames:
    strictVals[sampleName] = dict()
    for loc in snpValsPerFragment[sampleName]:
        nucleotideCnts = snpValsPerFragment[sampleName][loc]
        if(sum(cnt!=0 for cnt in nucleotideCnts) == 1 and max(nucleotideCnts) > MIN_COVERAGE):
            strictVals[sampleName][loc] = nucleotides[nucleotideCnts.index(max(nucleotideCnts))]
            allSNPlocs.add(loc)

# Follow the same process for the outgroup (Tibetan sheep)
# It is done separately, and after all of the fragment samples,
# because we are only interested in taking locations mapped in the outgroup
# to locations that were already taken in at least one of the fragments
# (and the outgroup files are massive compared to the fragment files,
# most of which is irrelevant given this)
tsvFilename = "~/scrolls/SRR5149616-1_1_1_1_sub5_merged.Ovis_aries.local.sam.collapsedcov3.af_alt_cov.tsv.gz"
vcfFilename = "~/scrolls/SRR5149616-1_1_1_1_sub5_merged.Ovis_aries.local.sam.collapsed_SNP.vcf"

sampleName = os.path.basename(vcfFilename)[:6]
if(sampleName != os.path.basename(tsvFilename)[:6]):
    print("Weird... vcf and perbase aren't in sync -", vcfFilename, "/", tsvFilename)
    input()

if sampleName in snpValsPerFragment:
    snpVals = snpValsPerFragment[sampleName]
else:
    snpVals = {}

if sampleName in refValsPerFragment:
    refVals = refValsPerFragment[sampleName]
else:
    refVals = {}

print(i, "/", len(filesVCF), "-", vcfFilename)
i += 1

sampleID = sampleNames.index(sampleName)

foundInVCF = dict()

with open(vcfFilename) as vcfFile:
    line = vcfFile.readline()
    while(line.startswith("##")):
        line = vcfFile.readline()

    if(not line.startswith("#CHROM")):
        print("Unexpected line:", line)
        input()

    lineNum = 0
    for line in vcfFile:
        lineNum += 1
        if(lineNum % 1000000 == 0):
            print("VCF -", lineNum)

        splitLine = line.rstrip().split()
        chrom = splitLine[0]

        chromNum = chromNameToNum(chrom)
        if(chromNum == 0):
            continue

        loc = int(splitLine[1])

		# We skip locations that weren't mapped in any of the fragments
        if((chromNum, loc) not in allSNPlocs):
            continue

        if(not splitLine[7].split(";")[1].startswith("CS=")):
            print("Weird -", splitLine[7])

        nucleotideCnts = [int(val) for val in splitLine[7].split(";")[1][3:].split(",")]

        if((chromNum, loc) in foundInVCF):
            print("Weird -", (chromNum, loc), "already in foundInVCF, =", foundInVCF[(chromNum, loc)])
            input()

        foundInVCF[(chromNum, loc)] = sum(nucleotideCnts)

        if (chromNum, loc) not in snpVals:
            snpVals[(chromNum, loc)] = [0 for nuc in nucleotides]

        snpVals[(chromNum, loc)] = [sum(val) for val in zip(snpVals[(chromNum, loc)], nucleotideCnts)]

with gzip.open(tsvFilename) as tsvFile:
    tsvFile.readline()

    lineNum = 0
    for line in tsvFile:
        lineNum += 1
        if(lineNum > 1000000):
            print("TSV -", lineNum)

        splitLine = line.decode().rstrip().split()
        chrom = splitLine[0]
        loc = int(splitLine[1])
        (ref, alt) = splitLine[2:4]
        af = float(splitLine[4])
        cov = int(splitLine[5])

        chromNum = chromNameToNum(chrom)
        if(chromNum == 0):
            continue

		# We skip locations that weren't mapped in any of the fragments
        if((chromNum, loc) not in allSNPlocs):
            continue

        if ref==alt:
            if (chromNum, loc) not in snpVals:
                snpVals[(chromNum, loc)] = [0 for nuc in nucleotides]

            snpVals[(chromNum, loc)][nucleotides.find(ref)] += cov
            if (chromNum, loc) in refVals:
                if refVals[(chromNum, loc)] != ref:
                    refVals[(chromNum, loc)] = "X"
            else:
                refVals[(chromNum, loc)] = ref

snpValsPerFragment[sampleName] = snpVals
refValsPerFragment[sampleName] = refVals

strictVals[sampleName] = dict()
for loc in snpValsPerFragment[sampleName]:
    nucleotideCnts = snpValsPerFragment[sampleName][loc]
    if(sum(cnt!=0 for cnt in nucleotideCnts) == 1 and max(nucleotideCnts) > MIN_COVERAGE):
        strictVals[sampleName][loc] = nucleotides[nucleotideCnts.index(max(nucleotideCnts))]
        allSNPlocs.add(loc)

sampleNames = list(strictVals.keys())

# Dictionaries mapping pairs of samples (including the outgroup) to the...
samplesSame = {}			# ... weighted sum of locations mapped in both, to the same nucleotide value
samplesDif = {}				# ... weighted sum of locations mapped in both, but to different nucleotides
samplesSameCnt = {}			# ... number of locations mapped in both, to the same nucleotide value
samplesDifCnt = {}			# ... number of locations mapped in both, but to different nucleotides
samplesSameRatio = {}		# ... weighted sum of same nucleotides out of the total weighted sum of (same + different)
samplesSameCntRatio = {}	# ... number of locations with the same mapped nucleotide out of the number of locations mapped in both
# Note: locations where both samples had the same nucleotide mapped,
# and it is the one in the reference genome, aren't taken into account

# The weighting used here is meant to offset the bias that may be caused
# by close locations when counting absolute numbers of location similarities/differences
# i.e. two close locations that have the same/different mapped value in a pair of samples
# give us significantly less than twice as much "information" as a single location
# (because of linkage disequilibrium).
# The goal is offsetting this, by weighing in proportion to the surrounding
# "genomic area" that the location gives us novel "information" about.
# This is calculated in proportion to the distance (in cM) between a location
# and it's relevant "neighbors" - the closest mapped locations that
# also give us information.

samplePairs = list(itertools.combinations(sampleNames, 2)) + [(sample, sample) for sample in sampleNames]
for i in range(len(samplePairs)):
    print(i, "/", len(samplePairs))
    (sample1, sample2) = samplePairs[i]

	# Mapping of genomic locations to "strictly" mapped nucleotide values
    sampleVals1 = strictVals[sample1]
    sampleVals2 = strictVals[sample2]

	# Genomic locations
    sampleKeys1 = set(sampleVals1)
    sampleKeys2 = set(sampleVals2)

	# Locations where both samples had "strictly" mapped nucleotides
    intersectingKeysIncludingRefRef = sorted(sampleKeys1.intersection(sampleKeys2))

	# Take out "uninteresting" locations - cases where both samples had the genome reference mapped
    intersectingKeys = [key for key in intersectingKeysIncludingRefRef if not
						(sampleVals1[key] == sampleVals2[key] and
						key in refValsPerFragment[sample1] and
						refValsPerFragment[sample1][key] == sampleVals1[key] and
						key in refValsPerFragment[sample2] and
						refValsPerFragment[sample2][key] == sampleVals2[key])]

    same = 0
    dif = 0
    sameCnt = 0
    difCnt = 0
    for i in range(1, len(intersectingKeys)-1):
        prevKey = intersectingKeys[i-1]	# The closest relevant location to the "left"
        key = intersectingKeys[i]		# The location we are observing
        nextKey = intersectingKeys[i+1]	# The closest relevant location to the "right"

		# This calculation is relevant only when "neighboring" locations are on the same chromosome
        if prevKey[0] != key[0] or key[0] != nextKey[0]:
            continue

		# The genetic distance (in cM) between the location and its "neighbors" is approximated.
		# The "weight" for the location is the distance between its "left and right neighbors",
		# and it is cut-off at GENETIC_DISTANCE_CUTOFF (after which we assume linkage
		# disequilibrium isn't relevant enough to take into account).
        prevLoc = findGeneticDistance(*prevKey)
        nextLoc = findGeneticDistance(*nextKey)
        dist = min(nextLoc - prevLoc, GENETIC_DISTANCE_CUTOFF)

        if(dist < 0):
            continue

		# Update the relevant variables according to the observation
        if sampleVals1[key] == sampleVals2[key]:
            same += dist
            sameCnt += 1
        else:
            dif += dist
            difCnt += 1

    samplesSame[(sample1, sample2)] = same
    samplesDif[(sample1, sample2)] = dif
    samplesSameCnt[(sample1, sample2)] = sameCnt
    samplesDifCnt[(sample1, sample2)] = difCnt
    samplesSameRatio[(sample1, sample2)] = same/max(same+dif, 0.000001)	# The miniscule number is added to avoid dividing by 0
    samplesSameCntRatio[(sample1, sample2)] = sameCnt/max(sameCnt+difCnt, 0.000001)

# Output the results into a file, as a table showing the results for each pair of samples
outFile = open("vcfTsv-withoutRefRef-results-v4distsG-cutoff" + str(GENETIC_DISTANCE_CUTOFF) + ".csv", "w")

for table in [samplesSame, samplesDif, samplesSameCnt, samplesDifCnt, samplesSameRatio, samplesSameCntRatio]:
    origKeys = list(table.keys())
    for (sample1, sample2) in origKeys:
        table[(sample2, sample1)] = table[(sample1, sample2)]

    print("," + ",".join(sampleNames), file=outFile)
    for sample1 in sampleNames:
        print(sample1 + "," + ",".join(str(table[(sample1, sample2)]) for sample2 in sampleNames), file=outFile)
    print("", file=outFile)
    print("", file=outFile)

outFile.close()