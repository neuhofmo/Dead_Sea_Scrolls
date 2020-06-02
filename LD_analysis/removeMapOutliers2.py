##########################################################################################
# Calculate similarity while taking into consideration Likage Disequilibrium

# This script accompanies the paper "Illuminating Genetic Mysteries of the Dead Sea Scrolls"
# Author: Or Sagy
##########################################################################################

# Remove the outlying points from a mapping between phsyical location (bp) and genetic location (cM)
# Outlying points are defined as those who are mapped to a preceding genetic location
# compared to their physical neighbor
# Note: the input and output files are sorted, within each chromosome, in ascending location order
# The .map format is [chromosome #, ID, physical location, genetic location] separated by tabs

inF = open("remapped_Oarv4.50kSNP.map")
lines = inF.readlines()
newLines = lines.copy()

first = True
while first or len(newLines) != len(lines):
    first = False
    lines = newLines.copy()
    split = [0, 0, 0, 0]

    for line in lines:
        prevSplit = split
        split = line.rstrip().split()
		
		# If the mapped genetic location precedes the previous genetic location, it is removed
        if split[0] == prevSplit[0] and float(split[2]) < float(prevSplit[2]):
            split = prevSplit
            newLines.remove(line)


outF = open("remapped_Oarv4.50kSNP.removedOutliers.map", "w")
for line in newLines:
    outF.write(line)
outF.close()
