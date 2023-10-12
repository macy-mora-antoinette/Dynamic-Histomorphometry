/*
-------------------------------------------------------------------------------------
 ***************************   DYNAMIC HISTOMORPHOMETRY   ***************************
-------------------------------------------------------------------------------------
This code contains all of the different macros associated with Imagej that are needed 
to help analyze images of mouse bones for dynamic histomorphometry. There are 6 parts.
Note that each part is independent of the other and could be copy and pasted into an 
individual .ijm file of it's own. For editing/updating, it's recommended to work with
each code in it's own separate file.

To test these macros, open ImageJ, select Plugins>Macros>Install Macros.
Select this file. Now, when you open Plugins>Macros, you should be able
to see the different functions available below:
*/

/*
-------------------------------------------------------------------------------------
 ************************************   PART 1   ************************************
-------------------------------------------------------------------------------------

OPEN IMAGES: This macro lets you open all the TIF files in a single directory at once
*/

macro "Open Images" {
	// Select the directory  where all of the images are saved;
	dir=getDirectory("Choose a Directory");  
	
	// acquire a list of all of the images in the directory selected above;
	arr = getFileList(dir);
	list = Array.filter(arr, "tif")
	
	// Open all of the images
	for( k = 0; k < list.length; k++ ) { 
		open(dir + list[k]); 
		setSlice(3);
		//run("Enhance Contrast", "saturation=0.35");
	}
}

/*
-------------------------------------------------------------------------------------
 ************************************   PART 2   ************************************
-------------------------------------------------------------------------------------

MAKE COMPOSITE IMAGES: In histomorphometry, bones are labeled with
flourophores of different colors, which you image using different filters under a 
microscope. We assume here you are using the filter Texas Red (TXR) and FITCI (GFP).
If you have these images saved in separate folders (TXR and GFP respectively), then 
this will create the necessary composites you need combining both fluorescent labels.

*/
macro "Make Composites" {

// Select the directory (folder) where all of the image folders (samples) are saved;
dir=getDirectory("Choose a Directory");  

// acquire a list of all of the image folders in the directory selected above;
list = getFileList(dir);
main_list = Array.filter(list, "/");

// For loop to iterate over the entire length of the data list;
for( k = 0; k < main_list.length; k++ ) { 

	name = main_list[k];
	
	//get list of samples within GFP and Texas Red (TXR) folder ;
	file = dir + name;
	sample_file = getFileList(file);

	//open Overlay;
	open(file + sample_file[2]);
	title_overlay = getTitle();
	
	//don't know why third channel is saved but it's empty
	Stack.setChannel(3);
	run("Delete Slice", "delete=channel");
	
	//separate red and green channels
	run("Split Channels");
	red_title = "C1-" + title_overlay;
	green_title = "C2-" + title_overlay;
	selectWindow(red_title);
	run("RGB Color");
	selectWindow(green_title);
	run("RGB Color");
	
	//create new composite TIFF
	imageCalculator("Add create", red_title, green_title);
	run("Images to Stack", "use");
	title_Composite = substring(main_list[k], 0, name.length-1);
	saveAs("Tiff", dir + title_Composite  + ".tif");
	run("Close All");
}

}

/*
-------------------------------------------------------------------------------------
 ************************************   PART 3   ************************************
-------------------------------------------------------------------------------------

MEASURE LABELED SURFACES: This code segments the perimeters using multipoints. 
It is split into two parts.

PART 1: Functions that are used to segment a polygon using specified points 
along its coordinates

PART 2: For loops to iterate over directories containing the images and ROI sets. PART 2
uses the functions created in PART 1. The single and double labeled surfaces are measured
and saved in a csv.

*/

macro "Label Surfaces" {
	
	// --------------------------------------//
	// *************** PART 1 *************** //
	// --------------------------------------//
	
	function get_shortest_distance(xperim, yperim, xpoints, ypoints) {
		
		//check all perimeter coordinates to see which ones are closest to the points
		idx = newArray(xpoints.length);
		for (i = 0; i < xpoints.length; i++) {
			dists = newArray(0);
			for (j = 0; j < xperim.length; j++) {
				//compute distance of each point
				dist = sqrt(pow(xpoints[i] - xperim[j], 2) + pow(ypoints[i] - yperim[j], 2)); 
				dists = Array.concat(dists, dist); //add to distance array
			}
			
			//Get the minimum distance from the distance array
			jMinAr = Array.findMinima(dists, 0);
			idx[i] = jMinAr[0];
		}
		Array.sort(idx);
		return idx;
	}
	
	function arr_contains( array, vals ) {
		values = newArray(0);
		values = Array.concat(values, vals);
	    for (i=0; i<array.length; i++) {
	    	for (j=0; j<values.length; j++) { 
		        if ( array[i] == values[j] ) { return true; }
		    }
		}
		return false;
	}
	
	
	function arr_sum( array ) {
	    sum = 0;
	    for (i=0; i<array.length; i++) {
			sum = sum + array[i];
		}
		return sum;
	}
	
	function segment_perimeter(xPerim, yPerim, xPoints, yPoints, name, color) {
		/*
		 * this function will create ROIs that label either red or green
		 * 
		Inputs:
		xPerim and yPerim are an array of the coordinates from a polygon in the ROI manager.
		xPoints and yPoints are an array of the coordinates of multiple points along the 
		perimeter of the polygon.
		Name is a string used to label each segment uniformly. Color is string.
		
		Outputs:
		all_array are all the indices to label the red or green segments
		*/
		
		
		//get labeling coordinates
		labelname = name + color + "label";
		RoiManager.selectByName(labelname);
		getSelectionCoordinates(X, Y);
		
		//remove the coordinates already in the xPoints and yPoints
		X = Array.slice(X, xPoints.length, X.length);
		Y = Array.slice(Y, xPoints.length, Y.length);
		
		//get index of the points on perim using shortest distance
		idx_label = get_shortest_distance(xPerim, yPerim, X, Y);
		idx = get_shortest_distance(xPerim, yPerim, xPoints, yPoints);
		
		//Check which segments the labels are located
		all_array = newArray(0);
		perim_idx = Array.getSequence(xPerim.length);
		for (t = 0; t < idx.length; t++) {
	
			if (idx.length < t+2 ) { //indices are sorted so last index indicated end of perim
				seg1 = Array.slice(perim_idx, 0, idx[0]);
				seg2 = Array.slice(perim_idx, idx[t],perim_idx.length-1);
				seg = Array.concat(seg2, seg1);
			} else {
				seg = Array.slice(perim_idx, idx[t],idx[t+1]);
			}
			//rename based on labeling 
			if (arr_contains(seg,idx_label)) {
				
				//get all indices will labels
				all_array = Array.concat(all_array, seg);
				//add line with labeling
				segx = newArray(0);
				segy = newArray(0);
				for (z = 0; z < seg.length; z++) {
					segx = Array.concat(segx, xPerim[seg[z]]);
					segy = Array.concat(segy, yPerim[seg[z]]);
				}
				makeSelection("freeline", segx, segy);
				roiManager("Add");
				roiManager("select", roiManager("count")-1);
			    roiManager("rename", name + color + "-segment" + toString(t));
			}
		
		}
		return all_array;		
		
	}
	
	
	function double_or_single(XPerim, YPerim, idx_red, idx_green, name) {
		/*
		 * this function will create ROIs that label either double or single labels
		 * 
		Inputs:
		xPerim and yPerim are an array of the coordinates from a polygon in the ROI manager.
		idx_red is the output from segment_perimeter
		idx_green is the output from segment_perimeter
		Name is a string used to label each segment uniformly. Color is string.
		
		Outputs:
		This function will simply add to the ROI manager. No variable outputs created.
		*/
		
		lab_count = 0;
		num = XPerim.length; 
		for(f = 0; f < num; f++ ) {
		
		//DOUBLE LABEL
		if (arr_contains(idx_red, f) & arr_contains(idx_green, f) & (f < num) ) {
			segx = newArray(0); segy = newArray(0);
		  	while ( arr_contains(idx_red, f) & arr_contains(idx_green, f) & (f < num)) {
		  		segx = Array.concat(segx, XPerim[f]); segy = Array.concat(segy, YPerim[f]);
		  		f++;
		  		
		  	} 
			if (segx.length > 1) {
			  	makeSelection("freeline", segx, segy);
				roiManager("Add");
				roiManager("select", roiManager("count")-1);
			    roiManager("rename", name + "-dL" + toString(lab_count));
			    roiManager("deselect");
			    lab_count++;
		    }
		}
				
		//SINGLE LABELS
		else if ( ( (arr_contains(idx_red, f) & (!arr_contains(idx_green, f)) ) || ((!arr_contains(idx_red, f)) & arr_contains(idx_green, f)) ) & (f < num) ) {
			segx = newArray(0); segy = newArray(0);
			while ( ( (arr_contains(idx_red, f) & (!arr_contains(idx_green, f)) ) || ((!arr_contains(idx_red, f)) & arr_contains(idx_green, f)) ) & (f < num) ) {
		  		segx = Array.concat(segx, XPerim[f]); segy = Array.concat(segy, YPerim[f]);
		  		f++;
		  	} 
		  	
		  	//print("x and y coordinates");
		  	if (segx.length > 1) {
			  	makeSelection("freeline", segx, segy);
				roiManager("Add");
				roiManager("select", roiManager("count")-1);
			    roiManager("rename", name + "-sL" + toString(lab_count));
			    roiManager("deselect");
			    lab_count++;
		    	}
			}					
		
		} 
	}
	
	// --------------------------------------//
	// *************** PART 2 *************** //
	// --------------------------------------//
	
	// Select the folder where all of the image folders are saved;
	main_dir=getDirectory("Select the folder where all of the image folders are saved");
	
	// acquire a list of all of the images in the directory selected above;
	main_arr = getFileList(main_dir);
	main_list = Array.filter(main_arr, "/");
	
	for( m = 0; m < main_list.length; m++ ) { 
	
		dir = main_dir + main_list[m];	
		arr = getFileList(dir);
		list = Array.filter(arr, "_label.zip");
		
		for( l = 0; l < list.length; l++ ) {	
	
			
			//Open TIFF and RoiSet
			end = lengthOf(list[l])-lengthOf("_label.zip"); //Removes "zip" from name
			title = substring(list[l], 0, end); //get the name of the tif from RoiSet
			open(dir + title + ".tif"); //open TIF
			roiManager("Open", dir + list[l]); //Open RoiSet
			print(list[l]);
			//iterate over ROI set
			count = RoiManager.size;
			for( p = 0; p < count; p++ ) {
				
				roiName = RoiManager.getName(p);
				print(roiName);
				if (roiName.contains("area")) { //we are looking for ROIs labeled psarea or esarea
				
					//extract name
					end = lengthOf(roiName)-lengthOf("area"); //remove "area" from name
					name = substring(roiName, 0, end); //name should be ps or es
					
					//get perimeters coordinates and smooth (spline)
					RoiManager.selectByName(roiName);
					run("Fit Spline"); 
					getSelectionCoordinates(XPerim, YPerim);
					
					// add smooth perimeter
					roiManager('add'); //add fitted spline as new perimeter
					roiManager('select', RoiManager.size - 1);
					roiManager("rename", name + "-Perim");
					roiManager("deselect");
					
					//get coordinates for RED segments
					redName = name + "-red";
					if (RoiManager.getIndex(redName) == -1) { //check that red label exists
						idx_red = newArray(0);
						idx_red = Array.concat(idx_red, -1);
					} else {
						RoiManager.selectByName(redName);
						getSelectionCoordinates(XPointsred, YPointsred);
						idx_red = segment_perimeter(XPerim, YPerim, XPointsred, YPointsred, name, "-red");
						roiManager("deselect");
					}
										
					//get coordinates for GREEN segments
					greenName = name + "-green";
					if (RoiManager.getIndex(greenName) == -1) {//check that green label exists
						idx_green = newArray(0);
						idx_green = Array.concat(idx_green, -1);
					} else {
						RoiManager.selectByName(greenName);
						getSelectionCoordinates(XPointsgreen, YPointsgreen);
						idx_green = segment_perimeter(XPerim, YPerim, XPointsgreen, YPointsgreen, name, "-green");				
						roiManager("deselect");
					}
					
					//label red/green segments as either double or single labels
					double_or_single(XPerim, YPerim, idx_red, idx_green, name);
				}
			}	
		roiManager("save", dir + title + "_label-Measured" + ".zip");
		roiManager("deselect");
		roiManager("delete");
		run("Close All");
		}

	}
}


macro "Measure Labels" {
	
	function arr_sum( array ) {
	    sum = 0;
	    for (i=0; i<array.length; i++) {
			sum = sum + array[i];
		}
		return sum;
	}
	
	results_name = getString("Name your results file:", "Single_and_Double_label_results");
	
	// Select the folder where all of the image folders are saved;
	main_dir=getDirectory("Select the folder where all of the image folders are saved");
	
	// acquire a list of all of the images in the directory selected above;
	main_arr = getFileList(main_dir);
	main_list = Array.filter(main_arr, "/");
	
	//to store table values
	titles = newArray(0);
	Ps_Pm_all = newArray(0);
	Ps_Ar_all = newArray(0);
	Ps_sL_all = newArray(0);
	Ps_dL_all = newArray(0);
	Es_Pm_all = newArray(0);
	Es_Ar_all = newArray(0);
	Es_sL_all = newArray(0);
	Es_dL_all = newArray(0);
	
	for( m = 0; m < main_list.length; m++ ) { 

		dir = main_dir + main_list[m];	
		arr = getFileList(dir);
		list = Array.filter(arr, "_label-Measured.zip");
		
		for( l = 0; l < list.length; l++ ) {	
	
			
			//Open TIFF and RoiSet
			end = lengthOf(list[l])-lengthOf("_label-Measured.zip"); //Removes "zip" from name
			title = substring(list[l], 0, end); //get the name of the tif from RoiSet
			titles = Array.concat(titles, title);
			open(dir + title + ".tif"); //open TIF
			//run("Set Scale...", "distance=1 known=1.8872 unit=micron"); //keyence for 4x scale
			run("Set Scale...", "distance=2.03 known=1 unit=micron");//set scale from microsope
			roiManager("Open", dir + list[l]); //Open RoiSet
			
			//to store lengths that will be used to compute table values
			Ps_Pm = newArray(0);
			Ps_Ar = newArray(0);
			Ps_sL = newArray(0);
			Ps_dL = newArray(0);
			Es_Pm = newArray(0);
			Es_sL = newArray(0);
			Es_dL = newArray(0);
			Es_Ar = newArray(0);
			
			//iterate over ROIs to get table
			count = RoiManager.size;
			for( c = 0; c < count; c++ ) {
				roiName = RoiManager.getName(c);
				roiManager("select", c);
						
				if (roiName.contains("Perim")) {
					
					len = getValue("Perim.");
					ar = getValue("Area");
					
					
					//get total perimeter length for ps and es
					if (roiName.contains("Ps") || roiName.contains("ps")) {Ps_Pm = Array.concat(Ps_Pm, len); Ps_Ar = Array.concat(Ps_Ar, ar);}
					else if (roiName.contains("Es") || roiName.contains("es")) {Es_Pm = Array.concat(Es_Pm, len);Es_Ar = Array.concat(Es_Ar, ar); }
					
					//get length of double and single labels
				} else if (roiName.contains("sL") || roiName.contains("dL")) {
					len = getValue("Length");
					//for ps
					if (roiName.contains("Ps") || roiName.contains("ps")) {
						if (roiName.contains("sL")) {Ps_sL = Array.concat(Ps_sL, len); }
						else if (roiName.contains("dL")) {Ps_dL = Array.concat(Ps_dL, len); }				
					}
					//for es
					else if (roiName.contains("Es") || roiName.contains("es")) {
						if (roiName.contains("sL")) {Es_sL = Array.concat(Es_sL, len); }
						else if (roiName.contains("dL")) {Es_dL = Array.concat(Es_dL, len); }	
					}
				} 				
			} 				
	
			// compute metrics as raw sums and percents
			Ps_Pm_sum = arr_sum(Ps_Pm);
			Ps_Ar_sum = arr_sum(Ps_Ar);
			Ps_sL_sum = arr_sum(Ps_sL);
			Ps_dL_sum = arr_sum(Ps_dL);
			Es_Pm_sum = arr_sum(Es_Pm);
			Es_Ar_sum = arr_sum(Es_Ar);
			Es_sL_sum = arr_sum(Es_sL);
			Es_dL_sum = arr_sum(Es_dL);
			
			Ps_sL_percent = Ps_sL_sum / Ps_Pm_sum;
			Ps_dL_percent = Ps_dL_sum / Ps_Pm_sum;
			Es_sL_percent = Es_sL_sum / Es_Pm_sum;
			Es_dL_percent = Es_dL_sum / Es_Pm_sum;
			
			// Manually add all metrics to their own array (ps and es both get: area, sl %, dl %, perimeter,)
			Ps_Pm_all = Array.concat(Ps_Pm_all, Ps_Pm_sum);
			Ps_Ar_all = Array.concat(Ps_Ar_all, Ps_Ar_sum);
			Ps_sL_all = Array.concat(Ps_sL_all, Ps_sL_percent);
			Ps_dL_all = Array.concat(Ps_dL_all, Ps_dL_percent);
			Es_Pm_all = Array.concat(Es_Pm_all, Es_Pm_sum);
			Es_Ar_all = Array.concat(Es_Ar_all, Es_Ar_sum);
			Es_sL_all = Array.concat(Es_sL_all, Es_sL_percent);
			Es_dL_all = Array.concat(Es_dL_all, Es_dL_percent);	
			
			roiManager("deselect");
			roiManager("delete");
			run("Close All");		
		}	
	}
	
	//table of results
	Table.create("Single_and_Double_labels");
	Table.setColumn("Names", titles);
	Table.setColumn("Bone_Ar", Ps_Ar_all);
	Table.setColumn("Marrow_Ar", Es_Ar_all);
	Table.setColumn("Ps_Pm", Ps_Pm_all);
	Table.setColumn("Ps_sL", Ps_sL_all);
	Table.setColumn("Ps_dL", Ps_dL_all);
	Table.setColumn("Es_Pm", Es_Pm_all);
	Table.setColumn("Es_sL", Es_sL_all);
	Table.setColumn("Es_dL", Es_dL_all);
	Table.save(main_dir + results_name + ".csv");	
}

/*
-------------------------------------------------------------------------------------
 ************************************   PART 4   ************************************
-------------------------------------------------------------------------------------

MEASURE INTERLABEL THICKNESS: In order for this to work, make sure that each pair of
lines to be measured come one right after the other (i.e. right next to each other).
This code assumes all even indices are line 1 and all odd indices are the partnering
line 2
*/

macro "Measure Interlabel Thickness" {
	
	results_name = getString("Name your results file:", "Interlabel_results");
	
	// Select the directory (folder) where all of the image folders (samples) are saved;
	main_dir=getDirectory("Choose a Directory");
	
	// acquire a list of all of the images in the directory selected above;
	main_arr = getFileList(main_dir);
	main_list = Array.filter(main_arr, "/");
	
	sample_names = newArray(0);
	Ps_means = newArray(0);
	Es_means = newArray(0);
	for( m = 0; m < main_list.length; m++ ) { 
		dir = main_dir + main_list[m];
		print(main_list[m]);
		arr = getFileList(dir);
		list = Array.filter(arr, "_interlabel.zip");
		
		
		// iterate over each RoiSet
		for( l = 0; l < list.length; l++ ) { 
			
			//Open TIFF and RoiSet
			end = lengthOf(list[l])-15; //Removes "zip" from name
			title = substring(list[l], 0, end); //get the name based on the ROIset
			print(title);
			open(dir + title + ".tif"); //open TIF
			
			//set scale (known from microsope)
			//run("Set Scale...", "distance=1 known=1.8872 unit=micron"); //keyence for 4x scale
			run("Set Scale...", "distance=2.03 known=1 unit=micron");//set scale from microsope
			roiManager("Open", dir + list[l]); //Open RoiSet
			sample_names = Array.concat(sample_names, title);
		
			//Iterate over ROI manager
			count = RoiManager.size/2;
			ps_count = 1;
			es_count = 1;
			
			Ps_IL = newArray(0);
			Es_IL = newArray(0);
			for(r = 0; r < count; r++){
				
				//Get first line
				roiManager('select', 2*r);
				name = Roi.getName; //used later to check if Ps or Es (group 4 or 5)
				run("Fit Spline"); 
				getSelectionCoordinates(x1, y1);
				length1 = getValue('Length');
				roiManager('add'); //add fitted spline as new perimeter
				roiManager('select', RoiManager.size - 1);
				roiManager("rename", name + "-spline");
				roiManager("deselect");
				
				
				//Get second line
				roiManager('select', 2*r+1);
				name2 = Roi.getName;
				run("Fit Spline"); 
				getSelectionCoordinates(x2, y2);
				length2 = getValue('Length');
				roiManager('add'); //add fitted spline as new perimeter
				roiManager('select', RoiManager.size - 1);
				roiManager("rename", name2 + "-spline");
				roiManager("deselect");
				
				//deselect both lines once you have the coordinate values
				roiManager('deselect');
				
				// Identify which line is shorter (testLine) and longer (refLine);
				if(x1.length > x2.length){
				  refLineX = x1; refLineY = y1;
				  testLineX = x2; testLineY = y2;
				}else{
				  refLineX = x2; refLineY = y2;
				  testLineX = x1; testLineY = y1;
				}
				
				//Get indices for testing, use every other index on the testLineX
				p = Array.sequence( Math.round(testLineX.length/2-2) );   
				for(k = 0; k < p.length; k++){
					idx = p[k]*2;
					
					//check every point on the longer line to determine which points 
					//create the shortest path to the testline
					dists = newArray(0);
					for(j = 0; j < refLineX.length; j++){
					//compute distance of each point
						dist = sqrt(pow(testLineX[idx] - refLineX[j], 2) + pow(testLineY[idx] - refLineY[j], 2)); 
						dists = Array.concat(dists, dist); //add to distance array
					}
				
					//Get the minimum distance from the distance array
					jMinAr = Array.findMinima(dists, 0);
					jMin = jMinAr[0];
				    
				    //add line to ROI manager and add to Ps or Es array
				    //make new line with this distance
				    makeLine(refLineX[jMin], refLineY[jMin], testLineX[idx], testLineY[idx]); 
				    roiManager('add'); //add to the roiManager
				    roiManager("select", roiManager("count")-1);
				    IL_length = getValue('Length');
				    
				    //label names as Ps or Es
				    if (name.contains("Ps") | name.contains("ps")) { 
					    namePs = "Ps-IL-" + toString(ps_count);
					    roiManager("rename", namePs);
					    Ps_IL = Array.concat(Ps_IL,IL_length);
					    ps_count++;
				    } 
					if (name.contains("Es") | name.contains("es")) {
					    nameEs = "Es-IL-" + toString(es_count);
					    roiManager("rename", nameEs);
					    Es_IL = Array.concat(Es_IL,IL_length);
					    es_count++;
					} 
					
				}
				
			}		
			//Save ROIs
			roiManager("save", dir +  title + "-Interlabel-Measured.zip");		
			//get the average Ps and Es interlabel thickness
			Array.getStatistics(Ps_IL, min, max, mean_Ps);
			Array.getStatistics(Es_IL, min, max, mean_Es);
			Ps_means = Array.concat(Ps_means, mean_Ps);
			Es_means = Array.concat(Es_means, mean_Es);
			roiManager("deselect");
			roiManager("delete");
			run("Close All");
		}
		
	}
	
	//make table of results
	Table.create("Interlabel-thickness");
	Table.setColumn("Names", sample_names);
	Table.setColumn("Ps_IL", Ps_means);
	Table.setColumn("Es_IL", Es_means);
	Table.save(main_dir + results_name + ".csv")
}

/*
-------------------------------------------------------------------------------------
 ************************************   PART 5   ************************************
-------------------------------------------------------------------------------------

*/

macro "Save Areas [9]" {
	
		
	//dir=getDirectory("Choose a Directory");
	dir = "D:\\Macy Castaneda\\Dynamic Histomorphometry\\2_15_2023\\";
	
	//Open TIFF and RoiSet
	title_image = getTitle();
	title = substring(title_image, 0, lengthOf(title_image) - lengthOf(".tif")); //get the name of the tif from RoiSet
	list = roiManager("count");
		
	//determine periosteum and endosteum using size of areas
	areas = newArray(0);
	for( j = 0; j < roiManager("count"); j++ ) { 
		roiManager("select", j);
		area = getValue("Area");
		areas = Array.concat(areas,area);
		}
	
	//Get the Periosteum (Ps)
	areas = Array.sort(areas); //sort in ascending order
	roiManager("Deselect");
	roiManager("select", (areas.length-1)); //max is last element
	roiManager("rename", "psarea") //periosteum has the largest area
	
	//Get the Endosteum (Es)
	val = 0;
	for( b = 0; b < areas.length-1; b++ ) { //Uncommon but there may be more than one endosteal surface
		roiManager("Deselect");
		roiManager("select", b); //all other elements
		if (b==0) {
			name = "esarea";
		} else {
			name = "es2area";
			}
		roiManager("rename", name) //Endosteum is the smaller area
		}
	
	//save ROIs
	roiManager("save", dir + title + "_label.zip");
	roiManager("deselect");
	roiManager("delete");


}

/*	
-------------------------------------------------------------------------------------
 ************************************   PART 7   ************************************
-------------------------------------------------------------------------------------
*/
macro "Save interlabel ROIs" {
	
		
	//dir=getDirectory("Choose a Directory");
	dir = "D:\\Macy Castaneda\\Dynamic Histomorphometry\\2_16_2023\\";
	
	//Open TIFF and RoiSet
	title_image = getTitle();
	title = substring(title_image, 0, lengthOf(title_image) - lengthOf(".tif")); //get the name of the tif from RoiSet
	
	//save ROIs
	roiManager("save", dir + title + "_interlabel.zip");
	roiManager("deselect");
	roiManager("delete");


}

/*		
-------------------------------------------------------------------------------------
 ************************************   PART 8   ************************************
-------------------------------------------------------------------------------------

OPEN IMAGES THRESHOLDED: This macro lets you open all the TIF files in a single
directory at once and thresholds the green and red images for easier identification
of fluorescent labeling
*/

macro "Open Threshold Images" {
	// Select the directory  where all of the images are saved;
	dir=getDirectory("Choose a Directory");  
	
	// acquire a list of all of the images in the directory selected above;
	arr = getFileList(dir);
	list = Array.filter(arr, "tif")
	
	// Open all of the images
	for( k = 0; k < list.length; k++ ) { 
		open(dir + list[k]); 
		
		//Open TIFF and RoiSet
		end = lengthOf(list[k])-lengthOf(".tif"); //Removes "zip" from name
		title = substring(list[k], 0, end); //get the name of the tif from RoiSet
		roiManager("Open", dir + title + "_label.zip"); //open TIF

		
		n = roiManager("count");
		for (i = 0; i < n; i++) {

			roiManager("select", i);
			name = RoiManager.getName(i);
		
			if (name.contains("psarea") | name.contains("Psarea") ) {
			    run("Enlarge...", "enlarge=50");
			    roiManager("add");
				roiManager("select", i);
				setSlice(1);
				
				mu = getValue("Mean");
				std = getValue("StdDev");
				max = getValue("Max");
				thresh = round(mu + 3*std);
				setMinAndMax(thresh, max);
				
				setSlice(2);
				mu2 = getValue("Mean");
				std2 = getValue("StdDev");
				max2 = getValue("Max");
				thresh2 = round(mu2 + 3*std2);
				setMinAndMax(thresh2, max2);
				
				setSlice(3);
				}
		}
		roiManager("deselect");	
		roiManager("delete");
	
	} 
}	


macro "Open ROIS [0]" {
	//Open TIFF and RoiSet
	title = getTitle();
	end = lengthOf(title)-lengthOf(".tif");; //Removes "zip" from name
	name = substring(title, 0, end); //get the name based on the ROIset
	roiManager("open", name + "_interlabel.zip"); //open TIF
}

macro "Open ROIS [1]" {
	//Open TIFF and RoiSet
	title = getTitle();
	end = lengthOf(title)-lengthOf(".tif");; //Removes "zip" from name
	name = substring(title, 0, end); //get the name based on the ROIset
	roiManager("open", name + "_label.zip"); //open TIF
}

	
/*
-------------------------------------------------------------------------------------
 ************************************   PART 9   ************************************
-------------------------------------------------------------------------------------
*/

macro "Label ps [2]" {
	idx = roiManager("index");
	roiManager("select", idx);
	roiManager("rename", "ps");
	}
	
macro "Label es [3]" {
	idx = roiManager("index");
	roiManager("select", idx);
	roiManager("rename", "es");
	}
	
macro "Label red [4]" {
	idx = roiManager("index");
	roiManager("select", idx);
	roi_name = Roi.getName;
	items = newArray("-red", "-green")
	
	//if it already has a label, replace it. Don't duplicate the label
	for( g = 0; g < items.length; g++ ) {
		if (roi_name.contains(items[g])) {
			end = lengthOf(roi_name)-lengthOf(items[g]);
			roi_name = substring(roi_name, 0, end); 
		}
	}
	//rename the ROI
	new_name = roi_name + "-red";
	roiManager("rename", new_name);
	}
	
macro "Label green [5]" {
	idx = roiManager("index");
	roiManager("select", idx);
	roi_name = Roi.getName;
	items = newArray("-red", "-green");
	
	//if it already has a label, replace it. Don't duplicate the label
	for( g = 0; g < items.length; g++ ) {
		if (roi_name.contains(items[g])) {
			end = lengthOf(roi_name)-lengthOf(items[g]);
			roi_name = substring(roi_name, 0, end); 
		}
	}
	//rename the ROI
	new_name = roi_name + "-green";
	roiManager("rename", new_name);
}

macro "Label label [6]" {
	idx = roiManager("index");
	roiManager("select", idx);
	roi_name = Roi.getName;

	
	//rename the ROI
	new_name = roi_name + "label";
	roiManager("rename", new_name);
}
	
