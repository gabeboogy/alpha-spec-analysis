import pandas as pd 
import urllib.request
import csv
#naming convention:
#nuclide-radtype.csv

path = "/home/user/Documents/VSCode/Alpha Spectroscopy 042026/isotopedata"
livechart = "https://nds.iaea.org/relnsd/v1/data?"

def lc_pd_dataframe(url):
    req = urllib.request.Request(url)
    req.add_header('User-Agent', 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')
    return pd.read_csv(urllib.request.urlopen(req))

def element_name(A,Z):
    with open("xxelements.csv", "r") as f:
        file = list(csv.reader(f))
        return f"{A}{file[Z-1][0]}"

def check(nuclide, rad_type):
    try:
        file = open(f"{path}/{nuclide}-{rad_type}.csv")
        file.close()
        return True
    except FileNotFoundError:
        return False

#returns a pd dataframe of nuclide-radtype
def df_open(nuclide, rad_type):
    file = f"{path}/{nuclide}-{rad_type}.csv"
    #returns existing dataframe if true, creates if not
    x = check(nuclide, rad_type)
    if x:
        print("Accessing Database...")
        return pd.read_csv(file)
    else:
        print("Accessing Livechart...")
        access = lc_pd_dataframe(f"{livechart}fields=decay_rads&nuclides={nuclide}&rad_types={rad_type}")
        access.to_csv(file)
        return pd.read_csv(file)
