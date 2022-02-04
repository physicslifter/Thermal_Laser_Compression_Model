#%%
from matplotlib import pyplot as plt
#Import function from run_model.py
import run_model

#Define runs to get
runs=['s88773',
      's88776',
      's88780',
      's88783'
      ]

#Define variables for the FEM model
num_steps=60
a=0.03
b=30
peak_temp=30000

#Run the model
output=run_model.main(runs, num_steps, a, b, peak_temp)

#Plot output
fig=plt.figure()

plt.subplot(2,2,1)
plt.plot(output['s88773']['face1'][0],output['s88773']['face1'][1], label='face 1')
plt.plot(output['s88773']['face2'][0],output['s88773']['face2'][1], label='face 2')
plt.plot(output['s88773']['face3'][0],output['s88773']['face3'][1], label='face 3')
plt.title('s88773')