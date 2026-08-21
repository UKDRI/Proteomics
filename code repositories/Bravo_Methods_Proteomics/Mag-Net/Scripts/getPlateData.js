var myRound = function (tmp) { return Math.round(tmp*10)/10 }
var myRef = function (tmp) { return (tmp ==="column") ? 8 : 1 }

labwareDataRetrievedOK = false

reagentsMaxVolDB = myRound(plateInfo(reagentsLabware).maxVolume) 
wasteMaxVolDB = myRound(plateInfo(wasteLabware).maxVolume)
beadsMaxVolDB = myRound(plateInfo(beadsLabware).maxVolume)
sampleMaxVolDB = myRound(plateInfo(sampleLabware).maxVolume)
plasmaMaxVolDB = myRound(plateInfo(plasmaLabware).maxVolume)
destinationMaxVolDB = myRound(plateInfo(destinationLabware).maxVolume)

reagentsMaxVol = myRound(reagentsMaxVolDB / myRef(reagentsMaxVolRef))
wasteMaxVol = myRound( wasteMaxVolDB / myRef(wasteMaxVolRef) )
beadsMaxVol = myRound( beadsMaxVolDB / myRef(beadsMaxVolRef) ) 
sampleMaxVol = myRound( sampleMaxVolDB / myRef(sampleMaxVolRef) )  
plasmaMaxVol = myRound( plasmaMaxVolDB / myRef(plasmaMaxVolRef) )  
destinationMaxVol = myRound( destinationMaxVolDB / myRef(destinationMaxVolRef) )  

labwareDataRetrievedOK = true


