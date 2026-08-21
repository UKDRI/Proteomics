function alert(){
	//var soundCommand = pow = "powershell -c (New-Object Media.SoundPlayer 'c:\\vworks workspace\\SP3\\Scripts\\Sounds\\alert.wav').PlaySync()"
	//var soundCommand = "C:\\VWorks Workspace\\Sounds\\alert.exe"
	var f = new File()
	var myFile1 = "C:\\VWorks Workspace\\Sounds\\alert.exe"
	var myFile2 = "C:\\VWorks Workspace\\Sounds\\alert.wav"
	if ( f.Exists(myFile1) && f.Exists(myFile2) ) {
		//run("cmd /C \"C:\\VWorks Workspace\\Sounds\\alert.exe\"")
		run("\"C:\\VWorks Workspace\\Sounds\\alert.exe\"")
   }
	else {
		print("Either program alert.exe or file alert.wav not found. Please check the sound directory")
	}
}	