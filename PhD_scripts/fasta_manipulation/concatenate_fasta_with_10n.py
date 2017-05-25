import sys

inpath = sys.argv[1]
infile = open(inpath,'r')
#name = sys.argv[3]
print inpath

incolumns = inpath.split('/')
last = len(incolumns)-1
remove_no = incolumns[last].index('.')
name = incolumns[last][:remove_no]
#print name

blank_folder = '/'.join(incolumns[0:last])
#print blank_folder
outpath = blank_folder + '/' + name + '_concat.fasta'
#print outpath

outfile = open(outpath,'w')

outfile.write('>concat_' + name + '\n')

for line in infile:
    if '>' in line:
        outfile.write('nnnnnnnnnn')
    else:
        line = line.strip('\n')
        outfile.write(line)
outfile.write('\n')
infile.close()
outfile.close()
