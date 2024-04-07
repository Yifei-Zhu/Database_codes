import pandas as pd 
import numpy as np  
import json 
import requests 
import time 
 
cid = 2244 
url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=Chemical+and+Physical+Properties' 
 
req = requests.get(url) 
req 
 
pubc_json = json.loads(req.text) 
 
pubc_json['Record']['Section'][0]['Section'][0]['Section'][2] 
 
Section  = pubc_json['Record']['Section'][0]['Section'][1]['Section'] 

for i in range(len(Section)): 
    if Section[i]['TOCHeading'] == 'Boiling Point':
        for j in range(len(Section[i]['Information'])):
            strs = Section[i]['Information'][j]['Value']
            try:
                bp = strs['StringWithMarkup'][0]['String']
            except:
                bp = strs['Number'][0]

            print(bp)


