function formCalcVols () {
   var volArr = []
   var lostVol = 3000 // 3 ml added for filling the inv pyramids
   var volIncBeadsPlate = 1.15 // 15% more volume will be added to the volumes in the beads plate
   var volIncReagPlate = 1.10 // same
   var volPerCycle, volStart 

   formInfo = ""
   reagVolEtOH1 = "Empty"
   reagVolEtOH2 = "Empty"
   reagVolAcN1 = "Empty"
   reagVolAcN2 = "Empty"
   reagVolABC = "Empty"
   reagVolTFA = "Empty"
   reagVolTFA1 = "Empty"
   reagVolBeads = "Empty"
   reagVolTrypsin = "Empty"
   reagVolReduction = "Empty"
   reagVolAlkylation = "Empty"

   // modify volumes according to the check boxes
   var _volReduction = (doReduction) ? volReduction : 0
   var _volAlkylation = (doAlkylation) ? volAlkylation : 0
   var _volBeads = 0//(doBeadsAcN) ? volBeads : 0
   var _volBeadsAcN = (doBeads) ? volBeadsAcN : 0
   var _volEtOH = (doEtOH) ? volEtOH : 0
   var _volAcN = (doAcN) ? volAcN : 0
   //var _volAcN2 = (doAcN2) ? volAcN2 : 0
   var _volABC = (doABCTrypsin) ? volABC : 0
   var _volTrypsin = (doABCTrypsin) ? volTrypsin : 0
   var _volTFA = (doAcidificationRecovery) ? volTFA : 0
   var _volTFA1 = (doAcidificationRecovery) ? volTFA1 : 0

   // check addition to sample plate
   if (_volEtOH > sampleMaxVol) {
      formInfo = "Wash EtOH will exceed sample plate capacity" 
      return false
   }
   else if (_volAcN > sampleMaxVol) {
      formInfo = "Wash AcN will exceed sample plate capacity"
      return false
   }    
   else if ((_volABC+_volTrypsin) > sampleMaxVol) {
      formInfo = "Volume of ABC + Trypsin will exceed sample plate capacity"
      return false
   }    
   else if (_volTFA > sampleMaxVol) {
      formInfo = "Volume of TFA will exceed sample plate capacity"
      return false
   }       
   else if (_volTFA1 > sampleMaxVol) {
      formInfo = "Volume of TFA 1% will exceed sample plate capacity"
      return false
   }    
	// increase all volumes in the reagent plate by 10%
	_volEtOH *= volIncReagPlate
	_volAcN *= volIncReagPlate
   //_volAcN2 *= volIncReagPlate
	_volABC *= volIncReagPlate
	_volTFA *= volIncReagPlate
   _volTFA1 *= volIncReagPlate
	_volBeadsAcN *= volIncReagPlate

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
      while (cInd < cyclesEtOH)
      reagVolEtOH1 = ( volArr[0] <= lostVol ) ? "Empty" : volArr[0]
      reagVolEtOH2 = ( volArr[1] <= lostVol ) ? "Empty" : volArr[1]
   }
   else {
      reagVolEtOH1 = "Empty"
      reagVolEtOH2 = "Empty"
      reagVolEtOH3 = "Empty"
   }

 tmp =  _volAcN * cyclesAcN * nCols * 8 + lostVol
   if (tmp > reagentsMaxVol*8) {
      formInfo = "Total AcN will exceed reagents plate capacity"
      reagVolAcN = "*Error"
      return false
   }
   else {
      reagVolAcN = (tmp <= lostVol) ? "Empty" : tmp
   }

tmp =  _volBeadsAcN*8*nCols + lostVol
   if (tmp > reagentsMaxVol*8) {
      formInfo = "Total AcN will exceed reagents plate capacity"
      reagVolAcN2 = "*Error"
      return false
   }
   else {
      reagVolAcN2 = (tmp <= lostVol) ? "Empty" : tmp
   }
   
   tmp = _volABC * nCols * 8 + lostVol
    if (tmp > reagentsMaxVol*8) {
      formInfo = "Total ABC will exceed reagents plate capacity"
      reagVolABC = "*Error"
      return false
   }
   else {
      reagVolABC = (tmp <= lostVol) ? "Empty" : tmp
   }

    tmp = _volTFA * nCols * 8 + lostVol
    if (tmp > reagentsMaxVol*8) {
      formInfo = "Total TFA will exceed reagents plate capacity"
      reagVolTFA = "*Error"
      return false
   }
   else {
      reagVolTFA = (tmp <= lostVol) ? "Empty" : tmp
   }

tmp = _volTFA1 * nCols * 8 + lostVol
   if (tmp > reagentsMaxVol*8) {
     formInfo = "Total TFA 1% will exceed reagents plate capacity"
     reagVolTFA1 = "*Error"
     return false
  }
  else {
     reagVolTFA1 = (tmp <= lostVol) ? "Empty" : tmp
  }


   // check beads plate
   // incrementing the volumes in the beads plate well by 15%
   //tmp  = _volBeads * (nCols + 1)
   //if (tmp > beadsMaxVol) {
   //   formInfo = "Beads volume will exceed beads plate capacity"
   //   reagVolBeads = "*Error"
   //   return false
   //}
   //else {
   //   reagVolBeads = (!tmp) ? "Empty" : tmp * volIncBeadsPlate
   //}

  tmp  = _volTrypsin * nCols
   if (tmp > beadsMaxVol) {
      formInfo = "Trypsin volume will exceed beads plate capacity"
      reagVolTrypsin = "*Error"
      return false
   }
   else {
      reagVolTrypsin = (!tmp) ? "Empty" : Math.ceil(tmp * volIncBeadsPlate + 5)
   }

  tmp  = _volReduction * nCols
   if (tmp > beadsMaxVol) {
      formInfo = "Reductant volume will exceed beads plate capacity"
      reagVolReduction = "*Error"
      return false
   }
   else {
      reagVolReduction = (!tmp) ? "Empty" :  Math.ceil(tmp * volIncBeadsPlate + 10)
   }

  tmp  = _volAlkylation * nCols
   if (tmp > beadsMaxVol) {
      formInfo = "Alkylant volume will exceed beads plate capacity"
      reagVolAlkylation = "*Error"
      return false
   }
   else {
      reagVolAlkylation = (!tmp) ? "Empty" :  Math.ceil(tmp * volIncBeadsPlate + 10)
   }
   
	// if you get to here everything is ok!
	return true
}

function formSetAcN70 () {
	formInfo = ""

   //if (doBeadsAcN && volBeads <=0) {
   //   formInfo = "No beads? No party! Disabling beads + AcN addition."
   //   doBeadsAcN = false
	//   return false
   //}
   
   if (!doReduction || !doAlkylation || !doEtOH) {
	   formInfo = "Warning: addtion of some key compounds is disabled"
   }

   var _volReduction = (doReduction) ? volReduction : 0
   var _volAlkylation = (doAlkylation) ? volAlkylation : 0
   var _volEtOH = (doEtOH) ? volEtOH : 0
   
   var tmp = _volReduction + _volAlkylation + _volEtOH
   var volAfterAcN = tmp + tmp*(0.7/0.3)

   if (volAfterAcN  > sampleMaxVol) {
      volBeadsAcN = "!!!!!"
      formInfo = "Volume of AcN would exceed sample plate capacity! Vol = " + volAfterAcN.toFixed(2)
	   return false 
   }
   else {
      volBeadsAcN = tmp*(0.7/0.3)
	   return true 
   }
 }


print("formVolFuncs.js read")