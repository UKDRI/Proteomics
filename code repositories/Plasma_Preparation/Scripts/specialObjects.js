tipsLoc3 = {
   statusCols: [],
   updateCols: function() {
      print("received global statusColsString = " + statusColsString)
	  // making sure that the resulting array is for 96 tip boxes
	  // using slice for keeping only the first 12 elements
      this.statusCols = statusColsString.split(",").slice(0,12)
      // need to call colsFormUpdate function 
      // which works in the global context
      colsFormUpdate() 
   },
   checkCols: function(availCols) {
      var cleanFound, emptyFound, tmpString, fullString, tmpIndex
      var cleanCount = 0
      fullString = this.statusCols.join("")

      // 1 - how many clean columns are available on the left and are they >= availCols? 
      tmpString = "C"
      while (!fullString.indexOf(tmpString)) {
         cleanCount++
         tmpString +="C"
      }
      cleanFound = (cleanCount >= availCols) ? true : false

      // 2 - are there any "availCols" empty cols anywhere?  
      tmpString = ""
      for (var i = 0; i < availCols; i++) tmpString += "E"  
      // lastIndexOf() starts for the end of fullString
      tmpIndex = fullString.lastIndexOf(tmpString)
      emptyFound = (tmpIndex >=0) ? true : false

      if (cleanFound && emptyFound) {
         return {
            tipsOK: true,
            tipsOff: tmpIndex+1,
            Mcols: cleanCount - availCols,
            Ncols: availCols
         }
      }
      else {
		  alert() // play a sound if tips status is not valid
         return {
            tipsOK: false
         }
      }   
   },
   shiftCols: function shiftCols (M,N) {
      print("Before shiftCols: " + this.statusCols.join(",")) 
      // shift M cols by N and set the remaining N to empty
      for (var i = 0; i < M+N; i++) {
         this.statusCols[i] = (i>M-1) ? "E" : this.statusCols[i+N]
      }
      print("After shiftCols: " + this.statusCols.join(",")) 
      //write the status to the interface
      statusColsString = this.statusCols.join(",")
      colsFormUpdate() 
   },
   tipsOnUpdate: function (N) {   
      print("Before tipsOnUpdate: " + this.statusCols.join(",")) 
      for (i=0; i < N; i++) {
         this.statusCols[i] = "E"
      }
      print("After tipsOnUpdate: " + this.statusCols.join(",")) 
      //write the status to the interface
      statusColsString = this.statusCols.join(",")
      colsFormUpdate() 
   },
   tipsOffUpdate: function (P,N) {
      print("Before tipsOffUpdate: " + this.statusCols.join(",")) 
      for (i = P-1; i < P-1+N; i++) {
         this.statusCols[i] = "U"
      }
      print("After tipsOffUpdate: " + this.statusCols.join(",")) 
      //write the status to the interface
      statusColsString = this.statusCols.join(",")
      colsFormUpdate() 
   },
   tipsOffCleanUpdate: function (P,N) {
      print("Before tipsOffCleanUpdate: " + this.statusCols.join(",")) 
      for (i = P-1; i < P-1+N; i++) {
         this.statusCols[i] = "C"
      }
      print("After tipsOffCleanUpdate: " + this.statusCols.join(",")) 
      //write the status to the interface
      statusColsString = this.statusCols.join(",")
      colsFormUpdate() 
   }      
}



// factory function used for creating reagLoc5 
// that in turn calculates  "distance from flat bottom on Agilent 6 col reservoir 201284-100"
// using arrays to make it more general 
 
function createReagLoc5 ()  {
   // hide all the following in a closure
   // all the reag* taken from global context as they come from the form
   var rg = {
      volAcN : [(typeof reagVolAcN !== "undefined") ? reagVolAcN : 0],
      volAcN2 : [(typeof reagVolAcN2 !== "undefined") ? reagVolAcN2 : 0],
      volEtOH : [
         (typeof reagVolEtOH1 !== "undefined") ? reagVolEtOH1 : 0,
         (typeof reagVolEtOH2 !== "undefined") ? reagVolEtOH2 : 0
      ],
      volABC : [(typeof reagVolABC !== "undefined") ? reagVolABC : 0],
      volTFA : [(typeof reagVolTFA !== "undefined") ? reagVolTFA : 0],
      volTFA1 : [(typeof reagVolTFA1 !== "undefined") ? reagVolTFA1 : 0]
   }
   var ct = {
      iAcN : 0,
      iAcN2 : 0,
      iEtOH : 0,
      iABC : 0,
      iTFA : 0,
      iTFA1 : 0
   }
   // return an object containing the getStuff() method
   return {
      getStuff: function (reag,volAsp,debug,refillAfterCalc) {
		  // refillAfterCalc (bool) is used to reinstate the amount of volume taken during precoating
         var area = 1111 // base area of a column
         var lostVol = 3000 // volume "lost" in the base pyramids
		  var distanceFineTuning = -2 // mm
         var debug = debug || true //forcing debug 
		  var refillAfterCalc = refillAfterCalc || false // if not passed don't use it
         var dist, wellsel
         // remember that wellselection refers to 96-well plate
         var ws = { EtOH:[1], AcN2:[3],AcN:[5], ABC:[7], TFA:[9], TFA1:[11] } 
         var reagS = "vol"+reag
         var iS = "i"+reag
         volAsp = Math.round(volAsp*10)/10
         if (!rg.hasOwnProperty(reagS)) {
            print("reagLog5: property "+ reagS +" does not exist")
            task.pause()
            return -1 
         }
         var selReag = rg[reagS] // REMINDER: passed by reference
         var i
         do {
            // get current status of the counter
            i = ct[iS] // REMINDER: passed by value
            if (debug) { 
               print("selReag = " + selReag)
               print("i = " + i)
            }     
            if (debug) {
              print("volAsp = " + volAsp)
              print("selReag[i] = " + selReag[i])
            }			  
			  // using Math.ceil() so that eg 1999.34 becomes 2000
			  // and it does not stop because of roundings in maths
            if ( Math.ceil(selReag[i]-volAsp) >= lostVol ) {
              if (debug) print("Volume before aspirating " + reag+" = " + selReag[i])
              selReag[i] -= volAsp
              if (debug) print("Volume after aspirating " + reag+" = " + selReag[i])
              dist = distanceFineTuning + (selReag[i] - lostVol) / area
              dist = (dist > 0) ? dist : 0 
              wellsel = [[ 1, ws[reag][i] ]]
              if (refillAfterCalc) {
				   if (debug) print("Refilling " + reag + " by " + volAsp + " uL")
				   selReag[i] += volAsp // patches, patches, everywhere...
			    } 
              if (debug) print ("Aspirating " + reag + " at " + Math.round(dist*10)/10 + " mm over flat bottom with wellsel = " + wellsel)
              // showing current status in the form (I know it is ugly...)
              reagVolAcN = rg.volAcN[0]
              reagVolAcN2 = rg.volAcN2[0]
              reagVolEtOH1 = rg.volEtOH[0]
              reagVolEtOH2 = rg.volEtOH[1]
              reagVolABC = rg.volABC[0]
              reagVolTFA = rg.volTFA[0]
              reagVolTFA1 = rg.volTFA1[0]
              return { dist: dist, wellsel: wellsel }
            }
            else {
               ct[iS]++
            }
         }
         while (ct[iS] < rg[reagS].length)
         print("reagLoc5: if you see this something weird has happened!") 
         task.pause()
     }
   }
}

// now instantiate the object
var reagLoc5 = createReagLoc5 ()

