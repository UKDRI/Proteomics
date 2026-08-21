function plateInfo (plateName) {
   var baseC = ["","Microplate","Filter plate","Reservoir","Tip Wash Station","Pin tool","Tip box","Lid","Tip trash bin","AM cartridge rack"] 
   var wellG = ["","Round","Square"]
   var wellB = ["","Rounded","Flat","V-Shaped"]
   var myKey, f, reg
   var vworksPath = "C:\\Program Files (x86)\\Agilent Technologies\\VWorks\\VWorks.exe"

   var f = new File()
   if ( f.Exists(vworksPath) ) {
      // 64 bit vworks
      myKey = "SOFTWARE\\Wow6432Node\\Velocity11\\Shared\\Labware\\Labware_Entries\\" + plateName
   }
   else {
      // 32 bit vworks
	   myKey = "SOFTWARE\\Velocity11\\Shared\\Labware\\Labware_Entries\\" + plateName
   }

   // create registry object 
   var reg = new Registry () 

   // now return an object with all the required info
   return { name: reg.Read(myKey,"NAME"),
                wells: reg.Read(myKey,"NUMBER_OF_WELLS"), 
                maxVolume: reg.Read(myKey,"WELL_TIP_VOLUME"),
                labwareType: baseC[reg.Read(myKey,"BASE_CLASS")],
                wellDepth: reg.Read(myKey,"WELL_DEPTH"),
                wellDiameter: reg.Read(myKey,"WELL_DIAMETER"),
                wellGeometry: wellG[reg.Read(myKey,"WELL_GEOMETRY")],
                wellBottom: wellB[reg.Read(myKey,"WELL_BOTTOM_SHAPE")] }
}


   