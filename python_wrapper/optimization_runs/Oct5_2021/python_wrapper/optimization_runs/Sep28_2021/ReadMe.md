# Summary of Optimizations September 28, 2021

## Run 1
- initial parameters:
    -  params0=[36000.0,0.017,30,2.0*10**-8]
- changed tlist so that it begins back at 20ns
    - fits looked much better than they did yesterday. Have to look more into the Matlab code
    - fits are better than yesterday, but we're running into the age old problem

## Run 2
- initial parameters:
    -  params0=[36000.0,0.017,30,2.0*10**-8]
- same as run 1, but now we're using the Powell method
- optimization quit @ iteration 14, same thing as the Powell run from yesterday, it quit after trying a=1.017
    - need to look more around @ Powell method to see why the optimization attempts such drastic jumps

## Run 3
- initial parameters:
    -  params0=[36000.0,0.017,30,2.0*10**-8]
- Same as run1, but I've moved tlist back to 18ns - let's see what happens.
- I quit the run early during the meeting with Tyler, then worked on sfepy for the rest of the day, so I will redo this run tomorrow. It will be slightly different from today's run however, because of how I have written sfepy