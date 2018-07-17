import sys
import Levenshtein

def complement(base):
  if base == 'A':
    return 'T'
  elif base == 'C':
    return 'G'
  elif base == 'G':
    return 'C'
  elif base == 'T':
    return 'A'
  else:
    return 'N'

def reverse_complement(string):
  return ''.join(map(complement, string[::-1]))

def process_sam(filename):
  f = open(filename)
  reads = []
  for line in f:
    if line[0] == '@':
      continue
    fields = line.split()
    if int(fields[1]) & 4 == 4:
        # unmapped read
        continue
    reads.append(fields)
  f.close()
  return reads

def process_uncompressed(filename):
  f = open(filename)
  reads = []
  for line in f:
      # for now, exclude reads
      reads.append(line.split()[:-1])
  return reads

def compare(sam, uncompressed):
  sam_reads = process_sam(sam)
  reads = process_uncompressed(uncompressed)
  #assert(len(sam_reads) == len(reads))
  for i, (x, y) in enumerate(zip(sam_reads, reads)):
    if i % 100000 == 0:
      print "Checked " + str(i) + " examples"
    for k in xrange(min(len(x), len(y))):
      if x[k] != y[k]:
          print "Example: " + str(i)
          print "Field: " + str(k)
          print "Ref: " + x[k]
          print "Att: " + y[k]
          return

    """
    if (len(x) != len(y)):
        print "Wrong lengths at example " + str(i)
        assert(len(x) == len(y))

    value = Levenshtein.editops(x, y)
    dist = Levenshtein.distance(x, y)
    if dist != 0:
      print "Example: " + str(i)
      print dist
      print value
      print "Ref: " + x
      print "Att: " + y"""
  print "Checked " + str(i) + " examples in total"


compare(sys.argv[1], sys.argv[2])
