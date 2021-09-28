# A Description of the optimization runs carried out on September 27, 2021

## Run 1  
- Parameters
    - params0=[24000, 0.028, 30.0, 1.8*10**-8]
- differs from the last run on the previous day, because I changed tlist to be well before the start of the data

### Notes
- it seems to be wayyyy off to the left. I think this has something to do with the tlist being corrdcted to where it was before


## Run 2
- Parameters
    - params0=[24000, 0.028, 30.0, 1.8*10**-8]
- same as the previous run, but now I am using the Powell method
-powell method finds some sort of optimization, then quits abruptly
- didn't make notes of any of the fits from this run as it ran
- the parameters didn't seem to change much
    - only parameter it changed was peak_temp, then it spiked the a-value and the optimization quit. Current hypothesis is that it quit

## Run 3
- Parameters
- Using the new starting values Tyler sent me over slack

- Notes:
    - continues to optimize towards the same values, new starting times seemed to have a minimal effect on what the optimization algorithm does. Fits still appear to be far to the left, so I will have to retry these optimization techniques again, using the previous tlist