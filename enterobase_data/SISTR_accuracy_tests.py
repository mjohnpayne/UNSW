
from time import sleep as sl

data = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/SISTR/all_SISTR_data.txt","r").readlines()

serovar_list = ["Aberdeen","Chester","Infantis","Muenchen","Saintpaul","Stanley","Virchow","Waycross","Typhimurium","Enteritidis"]

stats = {}

for i in data:
    col = i.split('\t')
    serovar = col[19]
    sistr = col[33]
    if serovar != 0 and "(Predicted)" not in serovar and serovar != "" and col[2][:7] != "traces":
        if sistr in serovar or serovar in sistr:
            if serovar not in stats:
                stats[serovar] = {"count":0,"truepos":0,"falseneg":0,"falsepos":0}
                stats[serovar]["count"] = 1
                stats[serovar]["truepos"] = 1
            else:
                stats[serovar]["count"] +=1
                stats[serovar]["truepos"] +=1
        elif sistr not in serovar or serovar not in sistr:
            if serovar not in stats:
                stats[serovar] = {"count": 0, "truepos": 0, "falseneg": 0, "falsepos": 0}
                stats[serovar]["count"] = 1
                stats[serovar]["falseneg"] = 1
            else:
                stats[serovar]["count"] += 1
                stats[serovar]["falseneg"] += 1
            if sistr not in stats:
                stats[sistr] = {"count": 0, "truepos": 0, "falseneg": 0, "falsepos": 0}
                stats[sistr]["count"] = 1
                stats[sistr]["falsepos"] = 1
            else:
                stats[sistr]["count"] += 1
                stats[sistr]["falsepos"] += 1

nums = [0,0,0,0]

for i in stats:
    nums[0] += stats[i]["count"]
    nums[1] += stats[i]["truepos"]
    nums[2] += stats[i]["falseneg"]
    nums[3] += stats[i]["falsepos"]

print nums
