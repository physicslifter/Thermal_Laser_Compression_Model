import json

geometry_dict = {
    's88773': {
        'face1':1.07*10**(-6),
        'face2':2.22*10**(-6),
        'face3':3.07*10**(-6)
    },
    's88776': {
        'face1':1.07*10**(-6),
        'face2':2.22*10**(-6),
        'face3':3.07*10**(-6)
    },
    's88780': {
        'face1':0.45*10**(-6),
        'face2': 1.7*10**(-6),
        'face3':2.8*10**(-6)
    },
    's86483': {
        'face1': 0.76*10**(-6),
        'face2': 2.76*10**(-6),
        'face3': 4.76*10**(-6)
    }
}

default_values = {
    'tf' : 30*10**-9,
    'num_steps' : 60,     
    'init_temp' : 2000,
    'peak_temp' : 30000,
    'a' : 0.03,
    'b' : 30,
    'k_1' : 100,
    'rho' : 12000,
    'c' : 450,
    'MgO_length':2*10**(-6) #window length in microns
    }

geometry_json=json.dumps(geometry_dict)
default_values_json=json.dumps(default_values)
geo_f=open('geometry_dict.json', 'w')
val_f=open('default_values_dict.json', 'w')
geo_f.write(geometry_json)
val_f.write(default_values_json)
geo_f.close()
val_f.close()