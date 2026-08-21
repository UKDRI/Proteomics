function colsFormUpdate2 () {
	print(formFilePath)
   formSlashIndex = formFilePath.lastIndexOf("/")
   formPath = formFilePath.slice(0,formSlashIndex) + "/"
   print(formPath)
   
   // tentative for 384
   tipsGroup = tipsGroup || 1
   print("tipsGroup = " + tipsGroup)
   // end tentative code 
   

   blueImage =  formPath + "blueball.png"
   greenImage = formPath + "greenball.png"
   greyImage = formPath + "greyball.png"
   redImage = formPath + "redball.png"
   yellowImage = formPath + "yellowball.png"
   statusCols = statusColsString.split(",")

   // note I now use i + (tipsGroup-1) * 12 as index for statusCols
   // in a 384 tips box I have 4 groups
   // the first two are for tips on and the second 2 for tips off
   for (i = 0; i < 12; i++) {
	   tipsStep = 12 * (tipsGroup-1)
      if (statusCols[i+tipsStep] === "E") {
         eval("iCol"+(i+1)+" = yellowImage")
         eval("sCol"+(i+1)+" = 'E'")
      }
      else if (statusCols[i+tipsStep] === "C") {
			eval("iCol"+(i+1)+" = greenImage")
			eval("sCol"+(i+1)+" = 'C'")
      }
      else if (statusCols[i+tipsStep] === "U") {
			eval("iCol"+(i+1)+" = redImage")
			eval("sCol"+(i+1)+" = 'U'")
      }
   }
}


print("formTipsStuff read")
