function listMethods (rootPath) {
	var myPath = rootPath + "/Methods/"
	var myScriptPath = rootPath + "/Scripts/"
	var myFile = myPath + "methodsList"
	var fileList
	
	// var dosString = "\"cd " + myPath + "&& dir \/a-d \/b *.js > methodsList \""
	// var dosCommand =  "cmd \/c " + dosString 
	
	// using now batch file in case one has a bugged command line
	var dosCommand = "\"" + myScriptPath + "createMethodsList.bat" + "\""  

	// with "true" the filelist is correctly generated 
	// and printed to the log
	// however it is never shown in the droplist
	// without "true" (default = "false") one needs to 
	// click on "read list" twice... this is a known bug!!!

	run(dosCommand)
	
	var f = new File()	
	if (f.Exists(myFile)) {
		f.Open(myFile)
		fileList = f.Read().replace(/\n/g,";")
		print(fileList)
		f.Close()
		return (fileList || "No methods detected")
	}
	else {
		return "Error: Cannot find methodsList"
	}
}

function deleteMethod (rootPath, myMethod, doConfirm) {
	var myPath = rootPath + "/Methods/"
	var myFile = myPath + myMethod
	var dosString = "\"cd " + myPath + " && del \"" + myMethod + "\" \""
	var dosCommand =  "cmd \/c " + dosString
	print(dosCommand)
	if (!doConfirm) {
		return {status:false, info:"Can't delete: confirm box not checked" }
	}
	var f = new File()
	if (f.Exists(myFile)) {
		run(dosCommand)
		return {status: true, info:"Method " + myMethod + " deleted "} 
	}
	else {
		return {status:false, info:"Can't find method " + myMethod}
	}
}

function saveMethod (rootPath, myMethod, overwrite) {
	var myPath = rootPath + "/Methods/"
	var myScriptPath = rootPath + "/Scripts/"
	  
	// if the Methods directory does not exist, create it
	// using a batch file just in case someone's command line is bugged
	var dosCommand = "\"" + myScriptPath + "createMethodsDir.bat" + "\""
	run(dosCommand)
	
	var myFile = myPath + myMethod	+ ".js"
	
	// WARNING - no tips status written to the method
	var protVars = [
		"colTrypsin",
		"colBeads",
		"colAlkylation",
		"colReduction",
		"wasteTTForm",
		"tips6filtered",
		"tips3filtered",
		"tips6unfiltered",
		"tips3unfiltered",
		"sampleMaxVolRef",	
		"beadsMaxVolRef",
		"wasteMaxVolRef",
		"sampleLabware",
		"beadsLabware",
		"destinationLabware",
		"wasteLabware",
		"cyclesAcN",
		"cyclesEtOH",
		"volStart",
		"volSampTransf",
		"doRemoveBeadsAcN",
		"timeBeads",
		"timeABC",
		"timeAlkylation",
		"timeReduction",
		"offTempABC",
		"offTempAlkylation",
		"offTempReduction",
		"tempABC",
		"tempAlkylation",
		"tempReduction",
		"doAcidificationRecovery",
		"nCols",
		"doABCTrypsin",
		"volReduction",
		"doAcN",
		"volAlkylation",
		"volBeads",
		"doEtOH",
		"volBeadsAcN",
		"volEtOH",
		"doBeadsAcN",
		"volAcN",
		"doReduction",
		"doAlkylation",
		"volTFA",
		"stockTFA",
		"finalPercTFA",
		"volABC",
		"volTrypsin",
		"testMode",
		"doFinalShakeDestPlate"
	]
	
	var f = new File()
	var tmp, q 

	if ( (f.Exists(myFile) && overwrite) || !f.Exists(myFile) ) {
		f.Open(myFile,1,1)
		for (var i = 0; i < protVars.length ; i++) {
			tmp = this[protVars[i]]
			q = (typeof tmp === "string") ? "\"" : ""
			f.Write(protVars[i] + " = " + q + tmp + q + "\n")
		}
		f.Close()
		return {status: true, info: "Method " + myMethod + " saved "} 
	}
	else if (f.Exists(myFile)) {
		return {status: false, info: "Can't save: method " + myMethod + " exists "} 	
	}
}

function loadMethod (rootPath, myMethod) {
	var myPath = rootPath + "/Methods/"
	var myFile = myPath + myMethod
	var f = new File()
	if (f.Exists(myFile)) {
		open(myFile)
		return {status: true, info:"Method " + myMethod + " loaded "} 
	}
	else {
		return {status:false, info:"Can't find method " + myMethod}
	}
}
	