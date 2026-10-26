from packages import *
import packages

def file_ext(tf):
    temp_list = tf.split('.')
    ext = temp_list[len(temp_list)-1]
    return ext
    
def csv_opener(file):    #returns counts
    csvnumbers = []
    data = []
    counts = []
    with open(file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:                #reformats each row to be a list of string instead of nested list
            for i in row:
                x = ' '.join(row)
            csvnumbers.append(x)
        for i in csvnumbers:                 #appends only numeric values, ignoring spaces
            j = i.replace(' ','')
            if j.isnumeric():
                data.append(i.strip())
        for i in data:
            binsvscounts = [int(x) for x in i.split() if x.isdigit()]
            counts.append(binsvscounts[len(binsvscounts)-1])
        return counts, len(counts)

def spe_opener(file,bin_location=11):
    data = []

    with open(file, "r") as f:
        file = f.read()
        chopped_file = file.splitlines()
        bins = int(chopped_file[bin_location].split()[1]) + 1
        data = [int(chopped_file[i]) for i in range(bin_location+1,bin_location+1+bins)]

    return data, bins
    
def fileopener(file:str):
    ext = file_ext(file)
    #dedicated function chooser
    if ext.casefold() == "spe":               #I'm using $DATA as a marker
        return spe_opener(file)
    if ext.casefold() == "csv":
        return csv_opener(file)
