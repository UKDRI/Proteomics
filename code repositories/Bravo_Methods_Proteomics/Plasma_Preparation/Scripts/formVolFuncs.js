function formCalcVols () {
   var volArr = []
   var lostVol = 3000 // 3 ml added for filling the inv pyramids
   var volIncBeadsPlate = 1.15 // 15% more volume will be added to the volumes in the beads plate
   var volIncReagPlate = 1.10 // same
   var volPerCycle, volStart 

   formInfo = ""
   reagVolEtOH1 = "Empty"
   reagVolEtOH2 = "Empty"
   reagVolBeads = "Empty"


   // modify volumes according to the check boxes

   var _volBeads = (doBeads) ? volBeads : 0

   var _volEtOH = (doEtOH) ? volEtOH : 0

   //var _volAcN2 = (doAcN2) ? volAcN2 : 0


   // check addition to sample plate
   if (_volEtOH > sampleMaxVol) {
      formInfo = "Wash EtOH will exceed sample plate capacity" 
      return false
   }

	// increase all volumes in the reagent plate by 10%
	_volEtOH *= volIncReagPlate


// calculate volumes for EtOH
   if (_volEtOH > 0) {
      volPerCycle = _volEtOH * nCols * 8
      volArr = [lostVol,lostVol,lostVol] // the 3 columns start with lostVol
      vInd = 0
      cInd = 0
      vIndJustIncremented = true 
      do {
         if (volArr[vInd] + volPerCycle <= reagentsMaxVol*8) {
            volArr[vInd] += volPerCycle
            vIndJustIncremented = false
            cInd++ // increment cycle
         }
         else {
            if (vIndJustIncremented) {
               formInfo = "Wash EtOH will exceed reagents plate capacity"
               reagVolEtOH1 = "*Error"
               reagVolEtOH2 = "*Error"
               return false
            }
            else { 
               vInd++ // increment array position
               vIndJustIncremented = true
            }
         }
      }
      while (cInd < 1)
      reagVolEtOH1 = ( volArr[0] <= lostVol ) ? "Empty" : volArr[0]
      reagVolEtOH2 = ( volArr[1] <= lostVol ) ? "Empty" : volArr[1]
   }
   else {
      reagVolEtOH1 = "Empty"
      reagVolEtOH2 = "Empty"
      reagVolEtOH3 = "Empty"
   }


   
   // check beads plate
   // incrementing the volumes in the beads plate well by 15%
   tmp  = _volBeads * (nCols + 1)
   if (tmp > beadsMaxVol) {
      formInfo = "Beads volume will exceed beads plate capacity"
      reagVolBeads = "*Error"
      return false
   }
   else {
      reagVolBeads = (!tmp) ? "Empty" : tmp * volIncBeadsPlate
   }

}


print("formVolFuncs.js read")