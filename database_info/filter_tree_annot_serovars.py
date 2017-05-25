
from time import sleep as sl

infile = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/genometrakr_data/9-5-17_gtrackr_serovar_annot.tsv","r").readlines()

strain_list = ["aberdeen","birkenhead","chester","infantis","muenchen","saintpaul","stanley","virchow","waycross","typhimurium","enteritidis"]


for s in strain_list:
    outfile = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/genometrakr_data/9-5-17_gtrackr_"+ s +"_annot.txt","w")
    outfile.write(infile[0])
    for i in infile[1:]:
        col = i.strip('\n').split('\t')
        if s in str.lower(col[1]):
            outfile.write(col[0] + '\t' + str.lower(col[1]) + '\n')
        else:
            outfile.write(col[0] + '\t\n')
    outfile.close()



