import json
from datetime import datetime

def load_json(fname): #function for loading in json files
    with open(fname) as json_file: #assign opened file to object json_file
        data=json.load(json_file) #load data
    return data

def rewrite_globals(rewrite_dict):
    current_dict=load_json('global_variables.json')
    for key in list(rewrite_dict.keys()):
        current_dict[key]=rewrite_dict[key]
    new_dict_json=json.dumps(current_dict)
    new_dict_file=open('global_variables.json', 'w')
    new_dict_file.write(new_dict_json)
    new_dict_file.close()
    
#Get todays date
date=datetime.now()
year, month, day = str(date.year), str(date.month), str(date.day)
if len(day)==1:
    day='0'+day
if len(month)==1:
    month='0'+month
global_vars=load_json('global_variables.json')
opt_path = global_vars['optimization_data_path'][:-8]
new_opt_data_path=opt_path+str(year)+str(month)+str(day)
rewrite_dict={
    "optimization_data_path":new_opt_data_path
}
rewrite_globals(rewrite_dict)
    