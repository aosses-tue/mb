form Define pitch params
    sentence inputfile input.wav
    sentence outputfile output.txt
    positive timestep 0.01
    positive minf0 75
    positive numcand 15
    positive maxf0 1000 
    positive silence_thr 0.01
    positive voicing_thr 0.45
    positive octave_cost 0.1
    positive octave_jump_cost 0.35
    positive vUv_cost 0.14
endform

sep$ = "/"

soundId = Read from file... 'inputfile$'
pitchId = To Pitch (ac)...  timestep 75 15 no 0.01 0.45 0.1 0.35 0.14 1
#pitchId = To Pitch (ac)... timestep 75 numcand no silence_thr voicing_thr octave_cost octave_jump_cost vUv_cost 1

select pitchId
Write to short text file... 'outputfile$'
