function calcLoopsVols (wellVol, totCols, maxVolTips, solvent) {
   // wellVol: volume needed in each well
   // totCols: total columsn to process
   // maxVolTips: maximum volume one can aspirate with the tips
   // solvent: flag. If "true" then always dispense to all columns in several loops
   // 
   // To recap: if solvent === false then 
   // the routine calculate the loops based on then max vol of compound transferred in each step
   // if solvent  === true then 
   // it is assumed that one wants to dispense the solvent evenly across the whole number of cols
   //
   if (!wellVol || !totCols || !maxVolTips) {
      print("ERROR: no volume, no. of cols or max tip vol provided! Please abort!")
      err=true
      return
   }
   if (solvent === undefined) {
      print("ERROR: no solvent flag provided! Please Abort!")
      err = true
      return
   }
   var fillCols, num
   var o = {}
   var approx = function (num) {
		return (Math.round(num*100) /100)
   }
   // find how many columns can be filled with one single asp
   fillCols = Math.floor(maxVolTips/wellVol)
   if (fillCols >= totCols && !solvent) {
      o.loops = 1
      o.cols = totCols
      o.aspvol = approx(wellVol * totCols)
      o.dispvol = approx(wellVol)
      o.solventlike = false
   }
   else if (fillCols >= 1 && fillCols < totCols && !solvent) {
      // find the divisor 
      while (totCols%fillCols) fillCols-- 
      print("Number of columns to be filled with the same ASP = " + fillCols)
      o.loops = totCols / fillCols
      o.cols = fillCols
      o.aspvol = approx(wellVol * fillCols)
      o.dispvol = approx(wellVol)
	   o.solventlike = false
   }
   else if (fillCols < 1 || solvent) {
	   // the case fillcols < 1 is treated in the same way as when we dispense a solvent;	  
	   // this is useful for minimizing the total number of dispenses. 
	   // For example if wellVol is 80, totCols is 12 and tips contain 70 uL
	   // the number of loops will be 14 for solventlike dispense vs 24 when 40uL 
	   // are dispensed twice per well. 
      o.loops = Math.ceil(wellVol * totCols / maxVolTips)
      o.cols = totCols
      o.aspvol = approx(wellVol * totCols / o.loops)
      o.dispvol = approx(o.aspvol / totCols)
      o.solventlike = true 
   }
   return(o)
}

function calcLoopsVolsPrint (o) {
   print("loops = " + o.loops)
   print("cols = " + o.cols)
   print("aspvol = " + o.aspvol)
   print("dispvol = " + o.dispvol)
   print("solventlike = " + o.solventlike) 
}

function printTips(o,p) {
   print("task: " + o + " wellselection: " + p)
}

function printHeadmode (o,p) {
   print("task: " + o + " headmode: " + p)
}

function printArrayOfArray(a) {
var out = ""
   for (var i=0; i< a.length; i++) {
      out += "["+a[i][0]+","+a[i][1]+"]"
      if (i !== a.length-1) out += ","
   }
   return(out)
}


function colsFormUpdate () {
   // using global statusColsString
   var colsArr = statusColsString.split(",")
   
   // tentative for 384
   tipsGroup = tipsGroup || 1
   print("tipsGroup = " + tipsGroup)
   // end tentative code 
   
   for (i = 0; i < 12; i++) {
	   tipsStep = 12 * (tipsGroup-1)
      // when you need to input a new tips state during protocol
      // don't update selectors
      if (!tipsStatusChangeEnabled) eval("sCol"+(i+1)+ " = colsArr[i+tipsStep]")
      if (colsArr[i+tipsStep] === "E") {
         eval("iCol"+(i+1)+" = yellowImage")
      }
      else if (colsArr[i+tipsStep] === "C") {
         eval("iCol"+(i+1)+" = greenImage")
      }
      else if (colsArr[i+tipsStep] === "U") {
         eval("iCol"+(i+1)+" = redImage")
      }
   }
}


function myRound2 (n) {
   return Math.round(n*100)/100
}

function myRound1 (n) {
   return Math.round(n*10)/10
}


	    
