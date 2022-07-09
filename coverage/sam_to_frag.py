import os
import os.path
import sys
# import fileinput
import csv

path_barcode = sys.argv[1]
path_sam = sys.argv[2]
path_output = sys.argv[3]
print(path_barcode)
print(path_sam)
print(path_output)

with open(path_barcode) as f:
    reader = csv.reader(f)
    # next(reader) # skip header
    data = []
    for r in reader:
        # get the barcode list from the 1st column
        bclist = r[0]
        data.append(bclist)
with open(path_sam) as f:
    lines = f.read().splitlines()        
# total number of reads
tcount = 0
bcount = 0

# if os.path.exists(path_output):
#     sys.exit("The output file already exists.")
f = open(path_output, "w")
# for line in fileinput.input():
for line in lines:
    tcount += 1
    row = line.split("\t")
    rname = row[2]
    pos = int(row[3])
    mapq = row[4]
    seq = row[9]
    start = pos - 1
    end = pos + len(seq) - 1
    # line = line.strip()
    tags = row[11:] # tags = line.split()[11:]
    tags_sort = sorted(tags)
    for tag in tags_sort:
        tag_split = tag.split(':')
        if "CB" in tag_split:
            # get the last element
            barcode = tag_split[-1]
            if barcode in data:
                bcount += 1
                # print(str(tcount) + " " + str(bcount))
                output = rname + "\t" + str(start) + "\t" + str(end) + "\t" + barcode + "\t" + mapq
                print(output, file = f) # print >> f,line
f.close()            
pct = float(bcount)/float(tcount)
print(pct)
