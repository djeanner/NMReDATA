/******************************************************************************************************
Generate .nmredata.sdf from Mnova 11
by Damien Jeannerat, University of Geneva
http://nmredata.org
Starting program: assignmenetReport.qs for Mnova V11.0.1-17801
*****************************************************************************************************/

/*globals Str, Molecule, Application, NMRAssignments, NMRSpectrum, settings, print, MessageBox, MnUi, FileDialog, Dir, File, TextStream, AssignmentReporter*/
/*jslint plusplus: true, indent: 4*/

function assignmentReport() {
	'use strict';
	var mol, assign, spectra, foundMulti, specIndex, spectrum, multi, standardReporter, correlationsReporter,
		diag, drawnItems, tableText, pageItem, i, standardTable, correlationsTable, table, cols, j, width,
		addNumberOfNuclides, addMultiplicity, parameters, correlationsTableStart, assignmentsArray, format,
		lines = 25,
		dw = Application.mainWindow.activeDocument,
		clipBoardKey = "Correlation Reporter/ClipBoard",
		correlations2DKey  = "Correlation Reporter/2D Correlations",
		orderKey = "Correlation Reporter/Order by shift",
		decimalsForProtonKey = "Correlation Reporter/Number of decimals for proton",
		decimalsForCarbonKey = "Correlation Reporter/Number of decimals for carbon and x-nuclei",
		showShiftKey = "Correlation Reporter/Show Shift",
		includeMultiplicityKey = "Correlation Reporter/Include Multiplicity",
		addNumberOfNuclidesKey = "Correlation Reporter/Add number of nuclides",
		exportToFileKey = "Correlation Reporter/Export to File",
		exportingFormatKey = "Correlation Reporter/Exporting Format",
		dropLinesWithoutCorrelationKey = "Correlation Reporter/Drop Lines Without Correlation",
		formatKey = "Correlation Reporter/Format",
		showDeltaForCarbonKey = "Correlation Reporter/Show Delta for carbon",
		reportTxtFileKey = "Correlation Reporter/Report Txt File",
 		FileNameNmredata  = "out_export_nmredata.sdf.txt",
		reportHTMLFileKey = "Correlation Reporter/Report HTML File";

	function getActiveMolecule(aDocWin, aMolPlugin) {

		var mol = aMolPlugin.activeMolecule();
		if (mol.isValid()) {
			return mol;
		}

		if (aDocWin.itemCount("Molecule") === 1) {
			mol = new Molecule(aDocWin.item(0, "Molecule"));
			return mol;
		}
		return undefined;
	}

	function createHTMLReport(table, index, lines) {

		var i, j, output;

		function fillVoids(cell) {
			if (cell === undefined || cell === "") {
				return "-";
			}
			return cell;
		}


		output = "<font style=\"font-size: 8pt; font-family: Arial; color: black\">";
		output += "<html><head>";
		output += "<title>Correlations Table</title>";
		output += "</head><body>";
		output += '\n<table border="1" cellSpacing="0" cellPadding="4" width="100%">';
		output += '\n<tr style="background-color:silver">';


		for (i = 0; i < table[0].length; i++) {
			output += '<td><b>' + table[0][i] + '</b></td>';
		}
		output += '</tr>';

		if (index === undefined) {
			for (i =  1; i < table.length; i++) {
				output += '<tr>';
				for (j = 0; j < table[i].length; j++) {
					output += '<td>' + fillVoids(table[i][j]) + '</td>';
				}
				output += '</tr>';
			}
		} else {

			if (index !== 0) {
				index = index * lines;
			}

			for (i = index + 1; (i < table.length && (i <= (index + lines))); i++) {
				output += '<tr>';
				for (j = 0; j < table[i].length; j++) {
					output += '<td>' + fillVoids(table[i][j]) + '</td>';
				}
				output += '</tr>';
			}
		}
		output += '</table>';
		output += '</body></html>';

		return output;
	}



	function exportToFile(aFormat) {
//	function exportToFile(parameters) {

		function formatHeader(aHeader) {
			var re, header = aHeader.toString();
			re = new RegExp("<(.*?)>|[&;]", 'g');
			header = header.replace(re, '');
			return header;
		}

		var i, file, stream, dataFile;

//		if (aFormat) {
//			dataFile = FileDialog.getSaveFileName("*.txt", "Save report in .sdf", settings.value(reportTxtFileKey, Dir.home()));
			dataFile = Dir.home() + "/Mnova_table_of_correlations.sdf";
//		} else {
//			dataFile = FileDialog.getSaveFileName("*.html", "Save report in HTML format", settings.value(reportHTMLFileKey, Dir.home()));
//		}

		if (dataFile !== "") {

			file = new File(dataFile);

//			if (aFormat) {
				settings.setValue(reportTxtFileKey, dataFile);
//			} else {
//				settings.setValue(reportHTMLFileKey, dataFile);
//			}

			file.open(File.WriteOnly);
			stream = new TextStream(file);
			out_mol = mol.getMolfile();
			stream.writeln(out_mol);
			stream.writeln(">  <Mestre_correlation_table>");

//			if (aFormat) {
				for (i = 0; i < table.length; i++) {
////				for (i = 1; i < table.length; i++) {
					if (i === 0) {
						stream.writeln(formatHeader(table[i].join("\t")));
					} else {
						stream.writeln(table[i].join("\t"));
					}
					stream.flush();
				}
//			} else {
//				tableText = createHTMLReport(table);
//				stream.write(tableText);
//				stream.flush();
//			}
			stream.writeln("");
			stream.writeln("");
				stream.flush();
			file.close();
		}
	}



	function getCorrelationsArray() {

		return ["HSQC", "HMBC", "H2BC", "COSY", "NOESY", "TOCSY", "ROESY"];
	}

	function getCorrelationsDescriptions() {
		return [assign.realAssignedExp("HSQC"), "HMBC", "H2BC", assign.realAssignedExp("COSY"), assign.realAssignedExp("NOESY"), assign.realAssignedExp("TOCSY"), "ROESY"];
	}

	function getAssignmentsArray() {
		var headerArray = [];
		headerArray.push("No");

		if (!addNumberOfNuclides && !addMultiplicity) {
			headerArray.push("&delta; <sub>H</sub>");
		} else if (addNumberOfNuclides && !addMultiplicity) {
			headerArray.push("&delta; <sub>H</sub> (nH)");
		} else if (!addNumberOfNuclides && addMultiplicity) {
			headerArray.push("&delta; <sub>H</sub> (Mul, <i>J</i>)");
		} else if (addNumberOfNuclides && addMultiplicity) {
			headerArray.push("&delta; <sub>H</sub> (Mul, <i>J</i>, nH)");
		}
		return headerArray;
	}

	function getAssignmentsDescriptions() {
		var headerArray = [];
		headerArray.push("No");

		if (!addNumberOfNuclides && !addMultiplicity) {
			headerArray.push("&delta; <sub>H</sub>");
		} else if (addNumberOfNuclides && !addMultiplicity) {
			headerArray.push("&delta; <sub>H</sub> (nH)");
		} else if (!addNumberOfNuclides && addMultiplicity) {
			headerArray.push("&delta; <sub>H</sub> (Multiplicity, <i>J</i>)");
		} else if (addNumberOfNuclides && addMultiplicity) {
			headerArray.push("&delta; <sub>H</sub> (Multiplicity, <i>J</i>, nH)");
		}
		return headerArray;
	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (dw === undefined || Application.molecule === undefined) {
		return;
	}

	mol = getActiveMolecule(dw, Application.molecule);
	if (mol === undefined || !mol.isValid()) {
		MessageBox.critical("Invalid Molecule");
		return;
	}

	assign = new NMRAssignments(mol);
	spectra = dw.itemCount("NMR Spectrum");
	foundMulti = false;
	specIndex = 0;

	while (!foundMulti  && specIndex < spectra) {
		spectrum = new NMRSpectrum(dw.item(specIndex, "NMR Spectrum"));

		if ((spectrum.nucleus() === "1H") && (spectrum.dimCount === 1)) {
			foundMulti = true;
			multi = spectrum.multiplets();
		} else {
			specIndex++;
		}
	}
	if (!foundMulti) {
		MessageBox.critical("Invalid Spectrum");
		return;
	}
/*
	diag = Application.loadUiFile("ricares:assignmentReport.ui");

	diag.widgets.gb2DCorrelations.checked = settings.value(correlations2DKey, true);
	diag.widgets.sbDecimalsForProton.value = settings.value(decimalsForProtonKey, 2);
	diag.widgets.sbDecimalsForCarbon.value = settings.value(decimalsForCarbonKey, 1);
	diag.widgets.ckOrder.checked = settings.value(orderKey, true);
	diag.widgets.gbShowDeltaForCarbon.checked = settings.value(showShiftKey, true);
	diag.widgets.gbExportToFile.checked = settings.value(exportToFileKey, false);
	diag.widgets.rbText.checked = settings.value(exportingFormatKey, true);
	diag.widgets.rbHTML.checked = !diag.widgets.rbText.checked;
	diag.widgets.ckClipBoard.checked = settings.value(clipBoardKey, true);
	diag.widgets.ckIncludeMultiplicity.checked = settings.value(includeMultiplicityKey, false);
	diag.widgets.ckAddNumberOfNuclides.checked = settings.value(addNumberOfNuclidesKey, false);
	diag.widgets.ckDropLinesWithoutCorrelation.checked = settings.value(dropLinesWithoutCorrelationKey, false);
	diag.widgets.ckShowDeltaForCarbon.checked = settings.value(showDeltaForCarbonKey, true);
	format =  settings.value(formatKey, 1);


	if (format === 0) {
		diag.widgets.rbN.checked = true;
		diag.widgets.rbDeltaN.checked = false;
		diag.widgets.rbCnDelta.checked = false;
	} else if (format === 1) {
		diag.widgets.rbN.checked = false;
		diag.widgets.rbDeltaN.checked = true;
		diag.widgets.rbCnDelta.checked = false;
	} else if (format === 2) {
		diag.widgets.rbN.checked = false;
		diag.widgets.rbDeltaN.checked = false;
		diag.widgets.rbCnDelta.checked = true;
	}
*/
	drawnItems = [];

//	if (diag.exec()) {

		settings.setValue(correlations2DKey, true);// diag.widgets.gb2DCorrelations.checked);
		settings.setValue(decimalsForProtonKey, 4);//diag.widgets.sbDecimalsForProton.value);
		settings.setValue(decimalsForCarbonKey, 4);//diag.widgets.sbDecimalsForCarbon.value);
		settings.setValue(orderKey, true);//diag.widgets.ckOrder.checked);
		settings.setValue(showShiftKey, true);//diag.widgets.gbShowDeltaForCarbon.checked);
		settings.setValue(exportToFileKey, true);//diag.widgets.gbExportToFile.checked);
		settings.setValue(exportingFormatKey, true);//diag.widgets.rbText.checked);
		settings.setValue(clipBoardKey, true);//diag.widgets.ckClipBoard.checked);
		settings.setValue(includeMultiplicityKey, true);//diag.widgets.ckIncludeMultiplicity.checked);
		settings.setValue(addNumberOfNuclidesKey, true);//diag.widgets.ckAddNumberOfNuclides.checked);
		settings.setValue(dropLinesWithoutCorrelationKey, true);//diag.widgets.ckDropLinesWithoutCorrelation.checked);
		settings.setValue(showDeltaForCarbonKey, true);//diag.widgets.ckShowDeltaForCarbon.checked);

		/*if (diag.widgets.rbN.checked) {
			settings.setValue(formatKey, 0);
		} else if (diag.widgets.rbDeltaN.checked) {
			settings.setValue(formatKey, 1);
		} else if (diag.widgets.rbCnDelta.checked) {
			settings.setValue(formatKey, 2);
		}*/
		settings.setValue(formatKey, 0);

		addNumberOfNuclides = true;//diag.widgets.ckAddNumberOfNuclides.checked;
		addMultiplicity =  true;//diag.widgets.ckIncludeMultiplicity.checked;

		assignmentsArray = getAssignmentsArray();
		correlationsTableStart = assignmentsArray.length;

		standardReporter = new AssignmentReporter(assignmentsArray, "Main", getAssignmentsDescriptions(), "Correlation Reporter/H&C");
		correlationsReporter = new AssignmentReporter(getCorrelationsArray(), "2D Correlations", getCorrelationsDescriptions(), "Correlation Reporter/2D");

		parameters = {};
		parameters.protonDecimals = 4;//diag.widgets.sbDecimalsForProton.value;
		parameters.carbonDecimals = 4;//diag.widgets.sbDecimalsForCarbon.value;
		parameters.assignmentObject = assign;
		parameters.molecule = mol;
		parameters.reporter = standardReporter;
		parameters.multi = multi;
		parameters.addNumberOfNuclides = addNumberOfNuclides;
		parameters.addMultiplicity = addMultiplicity;
		parameters.showDeltaForCarbon = true;//diag.widgets.ckShowDeltaForCarbon.checked;
		parameters.FileNameNmredata = FileNameNmredata; 

		standardTable = AssignmentReporter.assignmentReport(parameters);

//		if (diag.widgets.gb2DCorrelations.checked) {
			parameters.reporter = correlationsReporter;

//			if (diag.widgets.rbN.checked) {
//				parameters.format = 0;
//			} else if (diag.widgets.rbDeltaN.checked) {
//				parameters.format = 1;
//			} else if (diag.widgets.rbCnDelta.checked) {
				parameters.format = 2;
//			}

			correlationsTable = AssignmentReporter.assignmentReportWithCorrelations(parameters);
//		}

//		if (diag.widgets.gbShowDeltaForCarbon.checked) {
			correlationsTableStart++;
//		}

		table = AssignmentReporter.getFinalTable(standardTable, correlationsTable);


//		if (diag.widgets.ckOrder.checked) {
			table = AssignmentReporter.getOrderedTable(table);
//		}

//		if (diag.widgets.ckShowDeltaForCarbon.checked) {
  // back//			table = AssignmentReporter.removeVoidAssignmentsRows(table, standardReporter.xNuclidesIndex);

//		} else {
			table = AssignmentReporter.removeVoidAssignmentsRows(table, 1);
//		}

//		if (diag.widgets.ckDropLinesWithoutCorrelation.checked && diag.widgets.gb2DCorrelations.checked) {
			table = AssignmentReporter.removeVoidCorrelationsRows(table, standardReporter.xNuclidesIndex);
//		}

//		table = AssignmentReporter.removeVoidColumns(table);


//		if (diag.widgets.ckClipBoard.checked) {
//
//			tableText = createHTMLReport(table);
//			pageItem = Application.draw.text(tableText, "Report Special", "Assignments Proton", true);
//			drawnItems.push(pageItem);
//
//			for (i = 1; i < drawnItems.length; i++) {
//				drawnItems[i].top = drawnItems[i - 1].top;
//				drawnItems[i].left = drawnItems[i - 1].right;
//			}
//			settings.setValue(clipBoardKey, true);
//			dw.setSelection(drawnItems);
//			Application.mainWindow.doAction("action_Edit_Copy");
//			Application.mainWindow.activeDocument.curPage().deleteItems(drawnItems);
//			dw.update();
//
//		} else {
//
//			cols = Math.ceil((table.length - 1) / lines);
//			for (i = 0; i < cols; i++) {
//
//				tableText = createHTMLReport(table, i, lines);
//				pageItem = Application.draw.text(tableText, "Report Special", "Assignments Proton", true);
//				drawnItems.push(pageItem);
//			}
//			for (j = 1; j < drawnItems.length; j++) {
////				drawnItems[j].top = drawnItems[j - 1].top;
//				width = drawnItems[j].width;
//				drawnItems[j].right = drawnItems[j - 1].right + width;
//				drawnItems[j].left = drawnItems[j - 1].right;
//				drawnItems[j].update();
//				dw.update();
//			}
//			settings.setValue(clipBoardKey, false);
//			Application.mainWindow.activeDocument.curPage().update();
//		}


//		if (diag.widgets.gbExportToFile.checked) {
			exportToFile(true);//exportToFile(diag.widgets.rbText.checked);
//		}
//	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function AssignmentReporter(aArray, aDescription, aCorrelationDesc) {
	'use strict';
	this.fCorrelations = aArray;
	this.fDescription = aDescription;
	this.fCorrelationsDescription = aCorrelationDesc;
	this.xNuclides = {"C": 2};
	this.xNuclidesIndex = 2;
	this.nucleids = {};
}

AssignmentReporter.prototype = {};


AssignmentReporter.getXNucleus = function (aElement, aAssignmentReporter) {
	'use strict';
	if (aAssignmentReporter.xNuclides[aElement] === undefined && aElement !== "H") {
		aAssignmentReporter.xNuclidesIndex++;
		aAssignmentReporter.xNuclides[aElement] = aAssignmentReporter.xNuclidesIndex;
	}
	return aAssignmentReporter.xNuclides[aElement];

};


AssignmentReporter.assignmentReport = function (parameters) {
	'use strict';

	var i, j, at, noEqHs, hIndex, atomLabel, atomRow, h, shift, hIsHeavyIndex, skip, shiftH, shiftH0, shiftH1, element, shifts, atomNH,
		aAssignmentObject = parameters.assignmentObject,
		aMolecule = parameters.molecule,
		aAssignmentReporter = parameters.reporter,
		aMultiplets = parameters.multi,
		aAddNumberOfNuclides = parameters.addNumberOfNuclides,
		aAddMultiplicity = parameters.addMultiplicity,
		aDecimals = parameters.protonDecimals,
		aDecimalsForCarbon = parameters.carbonDecimals,
		aCarbonAssignments = parameters.showDeltaForCarbon,
		mol=parameters.molecule,
		aCount = aMolecule.atomCount,
		tableRows = {},
		implicitH,//DJ add
		tagc = "", tagh = "",//DJ add
		lich = [], licc = [],//chemical shifts H and C
		lith = [], litc = [],//text associated to shifts H and C
		FileNameNmredata = parameters.FileNameNmredata,
		value_c_shift,//add DJ
		label = "",//add DJ
		counth = 0, countc = 0,//add dj
		ii , //add dj
                separ = ", ", nmredataLine,//DJ add
		file, stream, dataFile, outmol, //DJ add
		headerRow = [];

/// dj add to dump in file
//		dataFile = "/Volumes/san256/users_for_mac_system_macPro/jeannerat/Dropbox/mrc_working_group/work_on_NMR_record_generation/dj_mnova_generate_nmredata_sdf_file/toto.txt";
//		dataFile = "/Users/djeanner/Dropbox/mrc_working_group/work_on_NMR_record_generation/dj_mnova_generate_nmredata_sdf_file/toto.txt";
// 		dataFile = FileDialog.getSaveFileName("*.txt", "Save report in .sdf", settings.value(reportTxtFileKey, Dir.home()));

		dataFile = Dir.home() + "/Assignment_" + FileNameNmredata;
		if (dataFile !== "") {
			file = new File(dataFile);
			//settings.setValue(reportTxtFileKey, dataFile);
			file.open(File.WriteOnly);
			stream = new TextStream(file);
			//stream.writeln(Dir.home());
			//stream.writeln(FileNameNmredata);
			out_mol = mol.getMolfile();
			stream.writeln(out_mol);
			stream.writeln(">  <NMREDATA_ASSIGNMENT>");
			stream.flush();
		}
/// 


	if (aAssignmentReporter !== undefined) {
		for (i = 0; i < aAssignmentReporter.fCorrelations.length; i++) {
			headerRow.push(aAssignmentReporter.fCorrelationsDescription[i]);
		}
	}


	for (at = 1; at <= aCount; at++) {// loop over atoms…
		noEqHs = aAssignmentObject.notEqHs(at);
		skip = true;
		hIsHeavyIndex = false;
		atomLabel = aMolecule.atom(at).number;
		element = aMolecule.atom(at).elementSymbol;
		atomNH = aMolecule.atom(at).nHAll;
/// dj add to dump in file
		if (dataFile !== "") {
			stream.write(";                                        atom number : ");
			stream.writeln(at);
			stream.flush();
		}
/// 
		if (noEqHs.length === 0  && element !== "H") {

			atomRow = [];
			atomRow[0] = AssignmentReporter.atomIndexToString(atomLabel, at);
			shifts = [];
			atomRow[1] = "";
			nmredataLine="";
			if (aCarbonAssignments) {// parameters.showDeltaForCarbon
				shift =  aAssignmentObject.chemShiftArr(at);

				if (shift) {
					if (shift[1]) {
						shiftH0 = Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimalsForCarbon);
						shiftH1 = Number((shift[1].max + shift[1].min) / 2).toFixed(aDecimalsForCarbon);
						atomRow[AssignmentReporter.getXNucleus(element, aAssignmentReporter)] = shiftH0 + "," + shiftH1;
						nmredataLine= "C" + AssignmentReporter.atomIndexToString(atomLabel, at) + separ + Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimalsForCarbon) + separ +at + "      ;K1";
					} else {
						atomRow[AssignmentReporter.getXNucleus(element, aAssignmentReporter)] = Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimalsForCarbon);
						nmredataLine= "C" + AssignmentReporter.atomIndexToString(atomLabel, at) + separ + Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimalsForCarbon) + separ + at + "      ;K2" ;
					}
				} else {
					nmredataLine=";" +  AssignmentReporter.atomIndexToString(atomLabel, at) + separ + at + "delete this line atom type : " + element + "; K3 ";
					if (element !== "C") {
						atomRow[AssignmentReporter.getXNucleus(element, aAssignmentReporter)] = "";
					} else {
						atomRow[AssignmentReporter.getXNucleus(element, aAssignmentReporter)] = "-";
					}
				}
			}
			tableRows[atomRow[0]] = atomRow;
	/// dj add to dump in file
		if (dataFile !== "") {
			stream.write(";                                       a:");
			stream.writeln(atomRow);
			stream.writeln(nmredataLine);
			stream.flush();
		}
	//
		} else {

			for (hIndex = 0; hIndex < noEqHs.length; hIndex++) {
				atomRow = [];
				atomRow[0] = AssignmentReporter.atomIndexToString(atomLabel, at);
				atomRow[1] = "";
				shifts = [];
				h = noEqHs[hIndex];
				if (h === 0) {
					hIsHeavyIndex = true;//H not attached to any C
				}
				shift =  aAssignmentObject.chemShiftArr(at, h);

				nmredataLine= "";
				if (shift) {
					if  (element !== "H") {
						implicitH="H";
					}else{
						implicitH="";
					}
					if (noEqHs.length > 1) {
						atomRow[0] = AssignmentReporter.atomIndexToString(atomLabel, at, h, true);
					label="H" + AssignmentReporter.atomIndexToString(atomLabel, at, h, true);
					nmredataLine = label + separ + Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimalsForCarbon) + separ + implicitH + at + " ;check explicit H2 is OK ; L3";
					} else if (noEqHs.length > 0) {
						atomRow[0] = AssignmentReporter.atomIndexToString(atomLabel, at, h, false);
					label=AssignmentReporter.atomIndexToString(atomLabel, at, h, false);
//					nmredataLine = nmredataLine + "\n" + label + separ + Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimalsForCarbon) + separ + "H" + at + " ;check explicit H is OK ; L2";
					label="H" + AssignmentReporter.atomIndexToString(atomLabel, at, h, false);
					nmredataLine = label + separ + Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimalsForCarbon) + separ + implicitH + at + " ;check explicit H is OK ; L2";
					}
					skip = false;

					if (shift[1]) {
						shiftH0 = Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimals);
						shiftH1 = Number((shift[1].max + shift[1].min) / 2).toFixed(aDecimals);
						shifts.push(shiftH0);
						shifts.push(shiftH1);

					} else {
						shiftH = Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimals);
						shifts.push(shiftH);
					}

					if (aAssignmentReporter.nucleids[atomRow[0] + "_" + shift]) {
						aAssignmentReporter.nucleids[atomRow[0] + "_" + shift]++;

					} else {
						if (noEqHs.length > 1) {
							aAssignmentReporter.nucleids[atomRow[0] + "_" + shift] = 1;
						} else {
							if (shift.length > 1 || hIsHeavyIndex) {
								aAssignmentReporter.nucleids[atomRow[0] + "_" + shift] = 1;
							} else {
								aAssignmentReporter.nucleids[atomRow[0] + "_" + shift] = atomNH;
							}
						}
					}
					atomNH = aAssignmentReporter.nucleids[atomRow[0] + "_" + shift];
					atomRow[1] = AssignmentReporter.findInformation(aDecimals, aMultiplets, shifts, atomNH, aAddNumberOfNuclides, aAddMultiplicity, label);
					lith[counth]= AssignmentReporter.findInformation(aDecimals, aMultiplets, shifts, atomNH, aAddNumberOfNuclides, aAddMultiplicity, label);
					lich[counth]= shift; 
					counth++;

				} else {
					atomRow[1] = "-";
				}

/// dj add to dump in file
		if (dataFile !== "") {
			stream.write(";                                       b:");
			stream.writeln(atomRow);
			stream.writeln(nmredataLine);
			stream.flush();
		}
	//

				if (aCarbonAssignments) {// parameters.showDeltaForCarbon
				nmredataLine ="";

					shift =  aAssignmentObject.chemShiftArr(at);
					if (!hIsHeavyIndex && shift && element !== "H") {
						skip = false;

						if (shift[1]) {
							shiftH0 = Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimalsForCarbon);
							shiftH1 = Number((shift[1].max + shift[1].min) / 2).toFixed(aDecimalsForCarbon);
							atomRow[AssignmentReporter.getXNucleus(element, aAssignmentReporter)] = shiftH0 + "," + shiftH1;
				nmredataLine =                       "C" + AssignmentReporter.atomIndexToString(atomLabel, at) + separ + shiftH0 + separ +  at + ";  LC1";
				nmredataLine = nmredataLine + "\n" + "C" + AssignmentReporter.atomIndexToString(atomLabel, at) + separ + shiftH1 + separ +  at + ";  LC2";
						} else {
							atomRow[AssignmentReporter.getXNucleus(element, aAssignmentReporter)] = Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimalsForCarbon);
					value_c_shift=Number((shift[0].max + shift[0].min) / 2).toFixed(aDecimalsForCarbon);
					nmredataLine = "C" + AssignmentReporter.atomIndexToString(atomLabel, at) + separ + value_c_shift + separ + at + "    ; LC min/max:" + shift[0].max + " " + shift[0].min;
					litc[countc]= value_c_shift; 
					licc[countc]= (shift[0].max + shift[0].min) / 2; 
					countc++;
						}
/// dj add to dump in file
if ( hIndex === 0 ) {// only write c with the first proton bound to it
		if (dataFile !== "") {
			stream.write(";                                       c:");
			stream.writeln(atomRow);
			stream.writeln(nmredataLine);
			stream.flush();
		}
}
	//
					} else {
						atomRow[AssignmentReporter.getXNucleus(element, aAssignmentReporter)] = "-";
					}
				}

				tableRows[atomRow[0]] = atomRow;
			}

		}
	}
/// dj add to dump in file
	if (dataFile !== "") {
			stream.writeln("");
			stream.writeln("");
			stream.writeln(">  <NMREDATA_1D_1H_NOTOK>");
			stream.writeln(";sometimes multiplet is not described");
			stream.writeln("Larmor= 111.1 ; .... ");
			//	stream.write(lich[ii]); could be used to sort....
	 		for (ii = 0; ii < counth; ii++) {
				stream.writeln(lith[ii]);
			}
			counth=0;
			stream.writeln("");
			stream.writeln("");
			stream.writeln(">  <NMREDATA_1D_13C_NOTOK>");
			stream.writeln("Larmor= 111.1 ; .... ");
			//	stream.write(licc[ii]); could be used to sort....
	 		for (ii = 0; ii < countc; ii++) {
				stream.writeln(litc[ii]);
			}
			counth=0;
			stream.writeln("");
			stream.writeln("");
			stream.flush();

		file.close();
	}
//////

	for (j in aAssignmentReporter.xNuclides) {
		if (aAssignmentReporter.xNuclides.hasOwnProperty(j)) {
			headerRow.push("&delta; <sub>" + j + "</sub>");
		}
	}
	tableRows.header = headerRow;
	return tableRows;
};

AssignmentReporter.assignmentReportWithCorrelations = function (parameters) {
	'use strict';

	var i, at, noEqHs, hIndex, atomRow, h, c, shift, atomLabel, element, correlations,
		aAssignmentObject = parameters.assignmentObject,
		aMolecule = parameters.molecule,
		aAssignmentReporter = parameters.reporter,
		aFormat = parameters.format,
		aProtonDecimals = parameters.protonDecimals,
		aCarbonDecimals = parameters.carbonDecimals,
		correlationsArray = [],
		aCount = aMolecule.atomCount,
		tableRows = {},
		file,dataFile,stream,
		FileNameNmredata = parameters.FileNameNmredata,
		headerRow = [],
		type,//add dj
		emptynmr = [],
		nmredata = [];


	function sortFunctionForFloats(a, b) {
		return parseFloat(a) - parseFloat(b);
	}

	tableRows.header = headerRow;
		dataFile = Dir.home() + "/2D_spectra" + FileNameNmredata;
//		dataFile = Dir.home() + "/Mnova_table_of_correlations.sdf";
		if (dataFile !== "") {
			file = new File(dataFile);
			file.open(File.WriteOnly);
			stream = new TextStream(file);
		}
/// 

	if (aAssignmentReporter !== undefined) {
		for (i = 0; i < aAssignmentReporter.fCorrelations.length; i++) {
			nmredata[i] = "";
			headerRow.push(aAssignmentReporter.fCorrelationsDescription[i]);
			emptynmr[i] = true; 
			type= aAssignmentReporter.fCorrelationsDescription[i];
			if (type === "HSQC") { nmredata[i] += ">  <NMREDATA_2D_13C_1J_1H>\n"; }
			else if (type === "HMBC") { nmredata[i] += ">  <NMREDATA_2D_13C_NJ_1H>\n"; }
			else if (type === "H2BC") { nmredata[i] += ">  <NMREDATA_2D_13C_2J_1H>\n"; }
			else if (type === "COSY") { nmredata[i] += ">  <NMREDATA_2D_1H_NJ_1H>\n"; }
			else if (type === "NOESY") { nmredata[i] += ">  <NMREDATA_2D_1H_D_1H>\n"; }
			else if (type === "TOCSY") { nmredata[i] += ">  <NMREDATA_2D_1H_TJ_1H>\n"; }
//			else { nmredata[i] += ";" + type + ";type of other spectra\n" + "**************************************************\n" + "**************************************************\n";}
			nmredata[i] += "CorrType=" + type + "\n";
			nmredata[i] += "Pulseprogram=" + "na ;to be filled....." + "\n";
			nmredata[i] += "Larmor=" + "111.1 ;to be filled....." + "\n";
		}
	}
	for (at = 1; at <= aCount; at++) {
		noEqHs = aAssignmentObject.notEqHs(at);
		atomLabel = aMolecule.atom(at).number;
		element = aMolecule.atom(at).elementSymbol;

		for (hIndex = 0; hIndex < noEqHs.length; hIndex++) {
			atomRow = [];
			atomRow.push("");
			h = noEqHs[hIndex];
			shift =  aAssignmentObject.chemShiftArr(at, h);

			if (aAssignmentReporter !== undefined) {
				for (c = 0; c < aAssignmentReporter.fCorrelations.length; c++) {
					correlations = AssignmentReporter.correlationToString(aAssignmentObject, aMolecule, aProtonDecimals, aCarbonDecimals, at, h, aAssignmentReporter.fCorrelations[c], aFormat);

					if (shift) {//add label in first column 
						if (noEqHs.length > 1) {
							atomRow[0] = AssignmentReporter.atomIndexToString(atomLabel, at, h, true);
						} else if (noEqHs.length > 0) {
							atomRow[0] = AssignmentReporter.atomIndexToString(atomLabel, at, h, false);
						}
					}

					if (correlations !== "") {

						emptynmr[c] = false; 
						correlationsArray = correlations.split(", ");
						correlationsArray.sort(sortFunctionForFloats);

						correlations = correlationsArray.toString();
						for (i=0 ; i < correlationsArray.length ; i ++) {
							nmredata[c] += correlationsArray[i];
							nmredata[c] += "/";
							nmredata[c] += "H" + atomRow[0] ;  
							nmredata[c] += "\n";
						}
					}

					atomRow.push(correlations);

				}
				tableRows[atomRow[0]] = atomRow;
			}
		}
	}
	if (dataFile !== "") {
		for (c = 0; c < aAssignmentReporter.fCorrelations.length; c++) {
			if ( emptynmr[c] === false ) {
			stream.writeln(nmredata[c]);
			stream.writeln();
			stream.flush();
			}
		}
		stream.flush();
		file.close;
	}
	return tableRows;
};



AssignmentReporter.removeVoidAssignmentsRows = function (table, lastxNuclidesIndex) {
	"use strict";
	var i, j,
		counter = {},
		newTable = [];

	for (i = 1; i < table.length; i++) {


		for (j = 1; j <= lastxNuclidesIndex; j++) {
			if (table[i][j] !== "") {
				counter[i] = true;
			}

		}
	}

	newTable.push(table[0]);
	for (i = 0; i < table.length; i++) {
		if (counter[i]) {
			newTable.push(table[i]);
		}
	}

	return newTable;
};


AssignmentReporter.removeVoidCorrelationsRows = function (table, startOfCorrelations) {
	"use strict";
	var i, j,
		counter = {},
		newTable = [];


	for (i = 1; i < table.length; i++) {
		for (j = parseInt(startOfCorrelations + 1, 10); j < table[i].length; j++) {
			if (table[i][j] !== "") {
				counter[i] = true;
			}
		}
	}

	newTable.push(table[0]);
	for (i = 0; i < table.length; i++) {
		if (counter[i]) {
			newTable.push(table[i]);
		}
	}

	return newTable;
};

AssignmentReporter.removeVoidColumns = function (table) {
	"use strict";

	var i, j, row,
		counter = [],
		newTable = [];

	for (i = 0; i < table[0].length; i++) {
		counter.push(0);
	}

	for (i = 1; i < table.length; i++) {
		for (j = 0; j < table[0].length; j++) {
			if (table[i][j] !== "" && table[i][j] !== undefined && table[i][j] !== "-") {
				counter[j]++;
			}
		}
	}

	for (i = 0; i < table.length; i++) {
		row = [];
		for (j = 0; j < counter.length; j++) {
			if (counter[j] !== 0) {
				row.push(table[i][j]);
			}
		}
		newTable.push(row);
	}
	return newTable;
};

AssignmentReporter.getFinalTable =  function (firstTable, secondTable) {
	"use strict";

	var i, j, joinedTable = [], aux = [];

	if (secondTable) {
		joinedTable.push(firstTable.header.concat(secondTable.header));
		for (i = 0; i < secondTable.header.length; i++) {
			aux.push("");
		}
	} else {
		joinedTable.push(firstTable.header);
	}

	for (i in firstTable) {
		if (firstTable.hasOwnProperty(i) && i !== "header") {
			for (j = 0; j < firstTable.header.length; j++) {
				if (firstTable[i][j] === undefined) {
					firstTable[i][j] = "";
				}
			}

			if (secondTable) {
				if (secondTable[i]) {
					joinedTable.push(firstTable[i].concat(secondTable[i].slice(1)));
				} else {
					joinedTable.push(firstTable[i].concat(aux));
				}
			} else {
				joinedTable.push(firstTable[i]);
			}
		}
	}

	return joinedTable;
};


AssignmentReporter.getOrderedTable =  function (tableRows) {
	'use strict';

	var i, k, a, b, temp,
		len = tableRows.length - 1;

	for (i = 1; i < len; i++) {
		for (k = 1; k < len; k++) {
			a = parseFloat(tableRows[k][1]);
			b = parseFloat(tableRows[k + 1][1]);
			if (isNaN(a)) {
				a = 0;
			}
			if (isNaN(b)) {
				b = 0;
			}
			if (a < b) {
				temp = tableRows[k + 1];
				tableRows[k + 1] = tableRows[k];
				tableRows[k] = temp;
			}
		}
	}
	return tableRows;
};

AssignmentReporter.getPpmArray =  function (table) {
	'use strict';
	var i,
		ppmArray = [];

	for (i = 1; i < table.length; i++) {
		ppmArray[i] = table[i][1];
	}
	return ppmArray;
};

AssignmentReporter.findInformation = function (decimals, multiplets, shifts, atomNH, addNumberOfNuclides, addMultiplicity, labeldj) {
	'use strict';

	function greaterThan(a, b) {
		return b - a;
	}

	var js, jArray, i, j, k, found, information = "", category, separ = ", ";

	for (j = 0; j < shifts.length; j++) {
		found = false;
		i = 0;
		information += shifts[j] + separ ;
		while (i < multiplets.count && !found) {

			if (multiplets.at(i).delta.toFixed(decimals) === shifts[j]) {
				found = true;

				if (addNumberOfNuclides || addMultiplicity) {
					//information += " (";

					if (addMultiplicity) {

						category = multiplets.at(i).category.toLowerCase();
						information += "S=" + category + separ;

					if (addNumberOfNuclides) {
						if (addMultiplicity) {
							//information += ", ";
						}

						//information += atomNH + "H";
						information += "N=" + atomNH + separ;

					}
					
						information += "L=" + labeldj ;
						if (category !== "m" && category !== "s") {
							js = multiplets.at(i).jList();
							information += separ + "J=";
							//information += ", ";
							jArray = [];
							for (k = 0; k < js.length; k++) {
								jArray[k] = js.at(k).toFixed(2);
							}

							jArray.sort(greaterThan); //sort descending

							for (k = 0; k < jArray.length; k++) {
								if (k === 0) {
									information += jArray[k];
								} else {
									information += "," + jArray[k];
								}
							}

							//information += " Hz";
						}
					}

					//information += "), ";
					information += ", ";
				} else {
					information += ", ";
				}
			} else {
				i++;
			}
		}

		if (!found) {
			information += ", ";
		}

	}

	if (information !== "") {
		information = information.substring(0, information.length - 2);
	}

	return information;
};

AssignmentReporter.atomIndexToString = function (aAtomLabel, aAtomNumber, aHIndex, aUseHIndex) {
	'use strict';
	var str = aAtomLabel;

	if (str === "") {
		str = aAtomNumber.toString();
	}
	if (aUseHIndex) {
		str += Str.hIndexToLetter(aHIndex);
	}
	return str;
};

AssignmentReporter.correlationToString = function (assignObject, aMolecule, protonDecimals, carbonDecimals, aAtom, aH, aCorrelation, aFormat) {
	'use strict';
	var i, cShift, cMeanShift, atomLabel,
		corrString = "",
		noEqAtoms,
		useHIndex,
		corrAtoms,
		precisionForShift;

	corrAtoms = assignObject.correlatedAtoms(aAtom, aH, aCorrelation);

	if (corrAtoms.length) {
		for (i = 0; i < corrAtoms.length; i++) {
			cShift = assignObject.chemShift(corrAtoms[i].atom, corrAtoms[i].indexes[0]);
			if (cShift !== undefined) {
				if (i > 0) {
					corrString  += ", ";
				}

				if (aFormat === 0) {
					atomLabel = aMolecule.atom(corrAtoms[i].atom).number;
					noEqAtoms = assignObject.notEqHs(corrAtoms[i].atom);

					if (corrAtoms[i].indexes.length === 1) {
						useHIndex = (noEqAtoms.length > 1);
						corrString += AssignmentReporter.atomIndexToString(atomLabel, corrAtoms[i].atom, corrAtoms[i].indexes[0], useHIndex);
					} else {
						corrString += atomLabel;
					}
				} else if (aFormat === 1) {

					atomLabel = aMolecule.atom(corrAtoms[i].atom).number;
					noEqAtoms = assignObject.notEqHs(corrAtoms[i].atom);
					cMeanShift = (cShift.max + cShift.min) / 2;

					if (corrAtoms[i].indexes[0] >= 1) {
						corrString += cMeanShift.toFixed(protonDecimals);

						if (corrAtoms[i].indexes[0] === 1 && noEqAtoms.length > 1) {
							atomLabel += "a";//"'"
						} else if (corrAtoms[i].indexes[0] === 2 && noEqAtoms.length > 1) {
							atomLabel += "b";//"''"
						}

						corrString += "(" + atomLabel + ")";
					} else if (corrAtoms[i].indexes.length === 1) {
						if (aMolecule.atom(corrAtoms[i].atom).elementSymbol === "C") {
							precisionForShift = carbonDecimals;
						} else {
							precisionForShift = protonDecimals;
						}
						corrString += cMeanShift.toFixed(precisionForShift);
						useHIndex = (noEqAtoms.length > 1);
						corrString += "(" + AssignmentReporter.atomIndexToString(atomLabel, corrAtoms[i].atom, corrAtoms[i].indexes[0], useHIndex) + ")";
					} else {
						corrString += cMeanShift.toFixed(protonDecimals);
						corrString += "(" + atomLabel + ")";
					}
				} else {

					atomLabel = aMolecule.atom(corrAtoms[i].atom).number;
					noEqAtoms = assignObject.notEqHs(corrAtoms[i].atom);



					if (corrAtoms[i].indexes[0] >= 1) {
						corrString += "H";

						corrString += atomLabel;

						if (corrAtoms[i].indexes[0] === 1 && noEqAtoms.length > 1) {
							corrString += "'";
						} else if (corrAtoms[i].indexes[0] === 2 && noEqAtoms.length > 1) {
							corrString += "''";
						}

						precisionForShift = protonDecimals;
					} else if (corrAtoms[i].indexes.length === 1) {

						corrString += aMolecule.atom(corrAtoms[i].atom).elementSymbol;
						useHIndex = (noEqAtoms.length > 1);
						corrString += AssignmentReporter.atomIndexToString(atomLabel, corrAtoms[i].atom, corrAtoms[i].indexes[0], useHIndex);
						if (aMolecule.atom(corrAtoms[i].atom).elementSymbol === "C") {
							precisionForShift = carbonDecimals;
						} else {
							precisionForShift = protonDecimals;
						}
					} else {
						corrString += "H";
						corrString += atomLabel;
						precisionForShift = protonDecimals;
					}
					cMeanShift = (cShift.max + cShift.min) / 2;
		//			corrString += "(" + cMeanShift.toFixed(precisionForShift) + ")";

				}
			}
		}
	}
	return corrString;
};

AssignmentReporter.findMultiplicity = function (decimals, multi, shift) {
	'use strict';

	function greaterThan(a, b) {
		return b - a;
	}

	var js, jArray, j,
		i = 0,
		found = false,
		multiplicity = "";

	while (i < multi.count && !found) {
		if (multi.at(i).delta.toFixed(decimals) === shift) {
			found = true;
			multiplicity = multi.at(i).category.toLowerCase();
			if (multiplicity !== "m" && multiplicity !== "s") {
				js = multi.at(i).jList();
				multiplicity += ",<i>J</i>=";
				jArray = [];
				for (j = 0; j < js.length; j++) {
					jArray[j] = js.at(j).toFixed(decimals);
				}

				jArray.sort(greaterThan); //sort descending

				for (j = 0; j < jArray.length; j++) {
					if (j === 0) {
						multiplicity += jArray[j];
					} else {
						multiplicity += "," + jArray[j];
					}
				}

				multiplicity += " Hz";
			}
		} else {
			i++;
		}
	}
	return multiplicity;
};


AssignmentReporter.assignmentReportCorrelations = function (decimals, aAssignmentObject, aMolecule, aAssignmentReporter, aMulti) {
	'use strict';

	//Deprecated

	var i, at, noEqHs, hIndex, atomLabel, atomRow, h, shift, hIsHeavyIndex, skip, shiftH, shiftH0, shiftH1, element,
		aCount = aMolecule.atomCount,
		tableRows = {},
		headerRow = [];

	if (aAssignmentReporter !== undefined) {
		for (i = 0; i < aAssignmentReporter.fCorrelations.length; i++) {
			headerRow.push(aAssignmentReporter.fCorrelationsDescription[i]);
		}
	}
	tableRows.header = headerRow;

	for (at = 1; at <= aCount; at++) {
		noEqHs = aAssignmentObject.notEqHs(at);
		skip = true;
		hIsHeavyIndex = false;
		atomLabel = aMolecule.atom(at).number;
		element = aMolecule.atom(at).elementSymbol;

		if (noEqHs.length === 0  && element === "C") {
			atomRow = [];
			atomRow[0] = AssignmentReporter.atomIndexToString(atomLabel, at);
			atomRow[1] = "";
			atomRow[2] = AssignmentReporter.findMultiplicity(decimals, aMulti, atomRow[1]);
			shift =  aAssignmentObject.chemShiftArr(at);
			if (shift) {
				if (shift[1]) {
					shiftH0 = Number((shift[0].max + shift[0].min) / 2).toFixed(decimals);
					shiftH1 = Number((shift[1].max + shift[1].min) / 2).toFixed(decimals);
					atomRow[3] = shiftH0 + "," + shiftH1;
				} else {
					atomRow[3] = Number((shift[0].max + shift[0].min) / 2).toFixed(decimals);
				}
			} else {
				atomRow[3] = "";
			}
			tableRows[atomRow[0]] = atomRow;
		} else {

			for (hIndex = 0; hIndex < noEqHs.length; hIndex++) {
				atomRow = [];
				atomRow[0] = AssignmentReporter.atomIndexToString(atomLabel, at);
				atomRow[1] = "";
				atomRow[2] = "";
				h = noEqHs[hIndex];
				if (h === 0) {
					hIsHeavyIndex = true;//H not attached to any C
				}
				shift =  aAssignmentObject.chemShiftArr(at, h);

				if (shift) {
					if (noEqHs.length > 1) {
						atomRow[0] = AssignmentReporter.atomIndexToString(atomLabel, at, h, true);
					} else if (noEqHs.length > 0) {
						atomRow[0] = AssignmentReporter.atomIndexToString(atomLabel, at, h, false);
					}
					skip = false;

					if (shift[1]) {
						shiftH0 = Number((shift[0].max + shift[0].min) / 2).toFixed(decimals);
						shiftH1 = Number((shift[1].max + shift[1].min) / 2).toFixed(decimals);
						atomRow[1] = shiftH0 + "," + shiftH1;

					} else {
						shiftH = Number((shift[0].max + shift[0].min) / 2).toFixed(decimals);
						if (atomRow[1] !== "") {
							shiftH = "," + shiftH;
							atomRow[1] += shiftH;
						} else {
							atomRow[1] = shiftH;
						}
					}
				}
				atomRow[2] = AssignmentReporter.findMultiplicity(decimals, aMulti, atomRow[1]);
				shift =  aAssignmentObject.chemShiftArr(at);
				if (!hIsHeavyIndex && shift && element === "C") {
					skip = false;

					if (shift[1]) {
						shiftH0 = Number((shift[0].max + shift[0].min) / 2).toFixed(decimals);
						shiftH1 = Number((shift[1].max + shift[1].min) / 2).toFixed(decimals);
						atomRow[3] = shiftH0 + "," + shiftH1;
					} else {
						atomRow[3] = Number((shift[0].max + shift[0].min) / 2).toFixed(decimals);
					}
				} else {
					atomRow[3] = "";
				}

				tableRows[atomRow[0]] = atomRow;
			}
		}
	}
	return tableRows;
};




if (this.MnUi && MnUi.scripts_nmr) {
	MnUi.scripts_nmr.scripts_nmr_ReportAssignments = assignmentReport;
}
