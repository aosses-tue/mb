# Praat version: 5.3.53 for Windows
# Tested in August 2013
# Last updated: 19/05/2014

form exportar
	sentence wavfile D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-sec-Harm\Stage5\model-ac-2-dist-rev.wav
	text filepath D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-sec-Harm\Stage5\model-ac-2-dist-rev.praat.txt
	positive speechFrameLength 0.01
	positive pitchFloor 75
	positive numCandidates 15
	text veryAccurate no
	positive silenceThreshold 0.01
	positive voicingThreshold 0.45
	positive octaveCost 0.01
	positive octaveJumpCost 0.35
	positive vUvCost 0.14
	positive pitchCeil 400
endform

if fileReadable(filepath$)
	pause File 'filepath$' Exists! Delete and continue?
	filedelete 'filepath$'
endif

Read from file... 'wavfile$'
sound = selected ("Sound")
tmin = Get start time
tmax = Get end time
#To Pitch (ac)... speechFrameLength pitchFloor numCandidates 'veryAccurate$' silenceThreshold voicingThreshold octaveCost octaveJumpCost vUvCost pitchCeil
To Pitch (ac)... speechFrameLength pitchFloor numCandidates 'veryAccurate$' silenceThreshold voicingThreshold octaveCost octaveJumpCost vUvCost pitchCeil

Rename... pitch
select sound
timeStep = speechFrameLength
clearinfo
echo Here are the results:
for i to (tmax-tmin)/timeStep
	time = tmin + i * timeStep
	select Pitch pitch
	pitch = Get value at time... time Hertz Linear
	printline 'time:4' 'pitch:4'
	fileappend 'filepath$' 'time:4' 'pitch:4' 'newline$' 
endfor
