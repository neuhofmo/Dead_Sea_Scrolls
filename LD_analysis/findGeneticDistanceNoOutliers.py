##########################################################################################
# For calculating the genetic distance between two points in the genome

# This script accompanies the paper "Illuminating Genetic Mysteries of the Dead Sea Scrolls"
# Author: Or Sagy
##########################################################################################

from bisect import bisect_left

chromToBpDistances = {i:[] for i in range(1, 28)}	# A dictionary from each chromosome to the ordered list of physical distances from the mapping
chromToCmDistances = {i:[] for i in range(1, 28)}	# A dictionary from each chromosome to the ordered list of genetic distances from the mapping

# Read the known mapping of physical distance to genetic distance
# The mapping is based on previous research (http://researchrepository.murdoch.edu.au/id/eprint/20751/2/02Whole.pdf)
# Its results had been remapped to OARv4 using NCBI Remap
# Outlying locations (i.e. mapped to a location beyond both locations before
# and after) had been removed - see removeMapOutliers.py
with open("remapped_Oarv4.50kSNP.removedOutliers.map") as f:
    for line in f:
        splitLine = line.split()
        chrom = int(splitLine[0])
        bpLoc = int(splitLine[3])
        cmLoc = float(splitLine[2])

		# Cases of locations not mapped to chromosomes 1-26 or chromosome X were ignored
        if chrom not in [0, 28]:
            if len(chromToBpDistances[chrom]) > 0 and cmLoc == 0:
                continue

            chromToBpDistances[chrom].append(bpLoc)
            chromToCmDistances[chrom].append(cmLoc)


def findGeneticDistance(chrom, bpDistance):
    """Approximate genetic distance from the beginning of a chromosome to a physical location
    Keyword arguments:
    chrom - the chromsome number (1-27, 27 being chromosome X)
    bpDistance - the physical location (bp from the 3' end of the chromosome)
    """
    bpDistances = chromToBpDistances[chrom]
    cmDistances = chromToCmDistances[chrom]

	# Find the index, in the physical distance list, that the point would be added after
    arrayLoc = bisect_left(bpDistances, bpDistance)

    if(arrayLoc <= 1):
        return 0

	# The mapped locations closest to the point that are to the "left" (3') of it
    leftBpDistance = bpDistances[arrayLoc-1]
    leftCmDistance = cmDistances[arrayLoc-1]

	# If the point is to the "right" (5') of all mapped locations -
	# extrapolate the location linearly based on the previous two mapped locations
    if(arrayLoc == len(bpDistances)):
        slope = (cmDistances[arrayLoc-1] - cmDistances[arrayLoc-2]) / (bpDistances[arrayLoc-1] - bpDistances[arrayLoc-2])
        return leftCmDistance + slope * (bpDistance - leftBpDistance)

	# Interpolate the genetic location of the point, between the mapped point to the
	# "left" (3') and to the "right" (5') of it
    rightBpDistance = bpDistances[arrayLoc]
    relativeBpDistanceDif = (bpDistance - leftBpDistance) / (rightBpDistance - leftBpDistance)
    rightCmDistance = cmDistances[arrayLoc]
    cmDistance = leftCmDistance + relativeBpDistanceDif * (rightCmDistance - leftCmDistance)

    return cmDistance
