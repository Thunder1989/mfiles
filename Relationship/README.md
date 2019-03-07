## Generate Event Sequence for Time Series using the Latent Event Model
Just run:
```
[~, Z, ~, ~, ~] = gibbs_sgf_K(data, Ks, F_switch, debug)
```
#data: is your time series input as a col vector,  
#Ks: is the total number of states,  
#F_swtich and debug are debug parameters, and can be set to 1 and 0 for now,  
and Z is the returned state sequence.

You might want to use the ```gibbs_sgf_K``` function with multiple cores, where each core handles a seperate file,
and an example on how to do this is in:
```
event_keti.m
```
