function beadsPlateColsCheck() {
	// check for duplicate values in a numeric array
	// there must not be two duplicate columns in the beads plate
	// only valid for numbers!
	var myCols = [colReduction,colAlkylation,colBeads,colTrypsin]
	var myArr = []
	for (var i=0 ; i<myCols.length ; i++) myArr[myCols[i]] = "F"
	print("string after myArr join = " + myArr.join(""))
	return (myCols.length === myArr.join("").length) ? true : false
}
	
