var canvas_width = 2500,
    controlpanel_width = 1300,
    geneGraph_width = canvas_width - controlpanel_width,
    transcriptRegion_height = 100,
    referenceLine_yPosition = transcriptRegion_height / 2,

    exon_height = 20,
    exon_yPosition = referenceLine_yPosition - exon_height / 2,
    exon_height_interval = 1,
    exon_editMode_height = exon_height + exon_height_interval,

    confirmButton_width = 100,
    confirmButton_height = 35,
    confirmButton_opacity = 0.1,
    cancelButton_width = confirmButton_width,
    cancelButton_height = confirmButton_height,
    cancelButton_opacity = confirmButton_opacity,

    editButton_width = 80,
    editButton_height = 40,
    editButton_xPosition = 1210,
    editButton_textX = editButton_xPosition + 13,
    editButton_yPosition = referenceLine_yPosition - editButton_height/2,
    editButton_textY = editButton_yPosition + editButton_height - 8,
    deleteButton_width = 80,
    deleteButton_height = 40,
    deleteButton_xPosition = 1295,
    deleteButton_textX = deleteButton_xPosition,
    deleteButton_yPosition = editButton_yPosition,
    deleteButton_textY = deleteButton_yPosition + deleteButton_height -8 ,
    addTranscript_textX = geneGraph_width/2 - 200,
addTranscript_textY = transcriptRegion_height/2 + 10,

    zeroOpacity = 0,
    fullOpacity = 1,
    greyout_opacity = 0.5,
    mouseover_opacity = 0.1,

    confirmbox_fontSize = 25;

var getExonOverlapDict = function(exons){
	    var exonOverlapDict = {};
	    for (var currentExon_index = 0; currentExon_index < exons.length; currentExon_index ++){
		for (var tmpExon_index = 0; tmpExon_index < exons.length; tmpExon_index ++){
		    if(currentExon_index != tmpExon_index){
			if( !( parseInt(exons[currentExon_index].start) >= parseInt(exons[tmpExon_index].end) || parseInt(exons[currentExon_index].end) <= parseInt(exons[tmpExon_index].start) ) ){
			    if(exonOverlapDict[currentExon_index]){
				exonOverlapDict[currentExon_index].push(tmpExon_index);
			    }
			    else{
				exonOverlapDict[currentExon_index] = new Array();
				exonOverlapDict[currentExon_index].push(tmpExon_index);
			    }
			}
		    }
		}
	    }
	    return exonOverlapDict;
};

var showTranscriptRegion = function(transcriptBackground, transcriptNr){
    transcriptBackground.mouseover = !transcriptBackground.mouseover;
    if(!transcriptBackground.clickState){
	d3.select("#transcriptBackground_" + transcriptNr)
	    .style("fill", "red")
	    .style("fill-opacity" , mouseover_opacity);
    }
};

var hideTranscriptRegion = function(transcriptBackground, transcriptNr){
    transcriptBackground.mouseover = !transcriptBackground.mouseover;
    if(!transcriptBackground.clickState){
	d3.select("#transcriptBackground_" + transcriptNr)
	    .style("fill", "red")
	    .style("fill-opacity", zeroOpacity);
    }
};

var chooseTranscriptRegion = function(transcriptBackground, transcriptNr, json, json_raw){
    transcriptBackground.clickState = !transcriptBackground.clickState;
    d3.select("#transcriptBackground_" + transcriptNr)
    .style("fill", "grey")
    .style("fill-opacity" , greyout_opacity);
    showTranscriptRegion(transcriptBackground, transcriptNr);

    if(transcriptBackground.clickState){
	// exclude this transcript and recalculate results
        // all we need to do is to set transcript.exons = []
	json.transcripts[transcriptNr].exons = [];
	if(debug_flag){
	    console.log(json);
	}
	getIsoformAbundance(json);
    }
    else{
	// include this transcript and recalculate results
        // add exons back to this transcript
	for (var i= 0; i < json_raw.transcripts[transcriptNr].exons.length; i++){
	    json.transcripts[transcriptNr].exons.push(json_raw.transcripts[transcriptNr].exons[i]);
	}
	getIsoformAbundance(json);
    }
};

var exonNotOverlap = function(exon_i, exon_j, exons){
    i_GreaterThan_j = parseInt(exons[exon_i].start) > parseInt(exons[exon_j].end);
    i_LessThan_j = parseInt(exons[exon_i].end) < parseInt(exons[exon_j].start);
    return (i_GreaterThan_j || i_LessThan_j);
};

var exonNotOverlapThisLevel = function(exon_i, exonLevel_j, exons){
    var notOverlapFlag = 1;
    for (var exon_index = 0; exon_index < exonLevel_j.length; exon_index ++){
	if(!exonNotOverlap(exon_i, exonLevel_j[exon_index], exons)){
	    notOverlapFlag = 0;
	}
    }
    return notOverlapFlag;
};

var getExonSizes = function(exons){
    var exonSizes = new Array();
    for(i = 0; i < exons.length; i++){
	exonSizes[i] = parseInt(exons[i].end) - parseInt(exons[i].start);
    }		
    return exonSizes;
};

var getExonOrders = function(exonSizes){
    var exonOrders = new Array();
    for(i = 0; i < exonSizes.length; i++){
	exonOrders[i] = 0;
	for (j =0 ; j<exonSizes.length;j++){
	    if(i != j){
		if(exonSizes[j] > exonSizes[i]){
		    exonOrders[i] += 1;
		}
		if(j < i){
		    if(exonSizes[j] == exonSizes[i]){
			exonOrders[i] += 1;
		    }
		}
	    }
	}
    }
    return exonOrders;
};

var getExonLevels = function(exons){
    var exonLevels = new Array();	    
    var exonSizes = getExonSizes(exons);
    var exonOrders = getExonOrders(exonSizes);

    var numberOflevels = 0;
    var needNewLevel = 1;
    for(i = 0; i< exonOrders.length; i++){ 
	needNewLevel = 1;
	if(i == 0){
	    exonLevels.push(new Array());
	    exonLevels[numberOflevels][i] = exonOrders.indexOf(i);
	}
	else{
	    for (j = 0; j<numberOflevels + 1;j++){
		if(exonNotOverlapThisLevel(exonOrders.indexOf(i), exonLevels[j], exons)){
		    exonLevels[j].push(exonOrders.indexOf(i));
		    needNewLevel = 0;
		    break;
		}
	    }
	    if(needNewLevel){
		numberOflevels += 1;
		exonLevels.push(new Array());
		exonLevels[numberOflevels][0] = exonOrders.indexOf(i);
	    }
	}
    }
    return exonLevels;
};

var getThisExonLevel = function (exonLevels, exon_index){
    for (i =0; i < exonLevels.length; i++){
	for (j=0; j<exonLevels[i].length; j++){
	    if(exonLevels[i][j] == exon_index){
		return i;
	    }
	}
    }
};

var colorOverlapExons = function(exon_index, exonOverlapDict){
    d3.select("#exon_edit_" + exon_index)
    .style("fill-opacity" , greyout_opacity);	
    if(exon_index in exonOverlapDict){
	for(j = 0; j < exonOverlapDict[exon_index].length; j++){
	    d3.select("#exon_edit_" + exonOverlapDict[exon_index][j])
		.style("fill-opacity" , fullOpacity);			    
	}
    }
};


var removeOverlapExons=	function (exons, exon_index, exonOverlapDict){
    if(exon_index in exonOverlapDict){
	for(j = 0; j < exonOverlapDict[exon_index].length; j++){
	    if(exons.indexOf(exonOverlapDict[exon_index][j]) != -1){
		exons.splice(exons.indexOf(exonOverlapDict[exon_index][j]) , 1);
	    }
	}
    }
};

var showEditRegion= function(d, transcriptNr, json){
    if(document.getElementById("editbox") != null){
	alert("Please confirm or cancel current editing task before proceed!");
    }
    else{
	var xScale = d3.scale.linear().domain([json.start, json.end]).range([0, geneGraph_width]);
	var xDiff = function(pt1, pt2) { return Math.abs(xScale(pt1) - xScale(pt2)); };

	var svg = d3.select("#geneGraph");
	var exonLevels = getExonLevels(json.exons);
	var exonOverlapDict = getExonOverlapDict(json.exons);
	var transcript_backup = jQuery.extend(true, {}, json.transcripts[transcriptNr]);

	var confirmButton_yPosition = exonLevels.length * exon_editMode_height + 10,
	confirmButton_textY = confirmButton_yPosition + confirmButton_height - 5,
	confirmButton_xPosition = geneGraph_width/2 - 20 - confirmButton_width,
	confirmButton_textX = confirmButton_xPosition + 8,
	cancelButton_yPosition = confirmButton_yPosition,
	cancelButton_textY = cancelButton_yPosition + cancelButton_height -5,
	cancelButton_xPosition = geneGraph_width/2 + 20,
	cancelButton_textX = cancelButton_xPosition + 15;

	var editbox_height = exonLevels.length * exon_editMode_height + 50;

	function redrawTranscript(transcriptNr){
	    var trans = d3.select("#transcript_" + transcriptNr);
	    trans.selectAll(".exon", function(exon){console.log(exon)}).remove();
	    trans
		.selectAll("rect.exon")
		.data(json.transcripts[transcriptNr].exons).enter()
		.append("rect")
		.attr("class", "exon")
		.attr("id", function(exonNr, arrIndex){return "exon_" + exonNr;})
		.attr("x", function(exon) {return xScale(parseInt(json.exons[exon].start));})
		.attr("y", exon_yPosition)
		.attr("width", function(exon) {
		    return xDiff( parseInt(json.exons[exon].end) , parseInt(json.exons[exon].start) )
		})
		.attr("height", exon_height);
	}

	function changeExon(exon, exon_index){
	    exon.clickState = !exon.clickState;
	    if(exon.clickState){
		if(json.transcripts[transcriptNr].exons.indexOf(exon_index) == -1){
		    removeOverlapExons(json.transcripts[transcriptNr].exons, exon_index, exonOverlapDict);
		    json.transcripts[transcriptNr].exons.push(exon_index);
		    redrawTranscript(transcriptNr);
		}
		colorOverlapExons(exon_index, exonOverlapDict);
		if(exon_index in exonOverlapDict){
		    for(i=0; i<exonOverlapDict[exon_index].length; i++){
			if(json.exons[exonOverlapDict[exon_index][i]].clickState){
			    json.exons[exonOverlapDict[exon_index][i]].clickState = !json.exons[exonOverlapDict[exon_index][i]].clickState;
			}		    
		    }
		}
	    }
	    else{
		json.transcripts[transcriptNr].exons.splice(json.transcripts[transcriptNr].exons.indexOf(exon_index), 1);
		redrawTranscript(transcriptNr);
		d3.select("#exon_edit_" + exon_index)
		    .style("fill-opacity" , fullOpacity);	
	    }
	}

	var allTranscripts = svg.selectAll("g");

	function restoreTranscript(){
	    json.transcripts[transcriptNr] = jQuery.extend(true, {}, transcript_backup);
	    redrawTranscript(transcriptNr);
	    for(exon_index =0; exon_index < json.exons.length; exon_index ++){
		if(json.exons[exon_index].clickState){
		    json.exons[exon_index].clickState = !json.exons[exon_index].clickState;
		}
	    }
	    d3.select("#editbox").remove();
	    if (transcriptNr != allTranscripts[0].length){
		for (j = transcriptNr + 1; j < allTranscripts[0].length; j++){
		    d3.select(allTranscripts[0][j]).attr("transform", "translate(0," + j*transcriptRegion_height +")");
		}
	    }
	}

	if (transcriptNr != allTranscripts[0].length){
	    for (j = transcriptNr + 1; j < allTranscripts[0].length; j++){
		d3.select(allTranscripts[0][j]).attr("transform", "translate(0,"+ (j*transcriptRegion_height + editbox_height) +")");
	    }
	}

	d3.select("#geneGraph")
	    .append("g")
	    .attr("id", "editbox")
	    .attr("transform", "translate(0," + (transcriptNr + 1)*transcriptRegion_height + ")")
	    .append("g")
	    .append("rect")
	    .attr("x", 0)
	    .attr("y", 0)
	    .attr("width", geneGraph_width)
	    .attr("height", editbox_height)
	    .attr("stroke", "red")
	    .attr("stroke-width", 2)
	    .style("fill-opacity" , zeroOpacity);

	d3.select("#editbox")
	    .append("g")
	    .attr("id", "confirmbox")
	    .append("text")
	    .attr("x", confirmButton_textX)
	    .attr("y", confirmButton_textY)
	    .text("Confirm")
	    .attr("font-size", confirmbox_fontSize);

	d3.select("#confirmbox")
	    .append("text")
	    .attr("x", cancelButton_textX)
	    .attr("y", cancelButton_textY)
	    .text("Cancel")
	    .attr("font-size", confirmbox_fontSize);

	d3.select("#confirmbox")
	    .append("rect")
	    .attr("id", "confirmButton")
	    .attr("x", confirmButton_xPosition)
	    .attr("y", confirmButton_yPosition)
	    .attr("width", confirmButton_width)
	    .attr("height", confirmButton_height)
	    .attr("stroke", "red")
	    .attr("stroke-width", 3)
	    .style("fill-opacity", zeroOpacity)
	    .on("click", function(confirm){startoverAlgo(json);});
	
	d3.select("#confirmbox")
	    .append("rect")
	    .attr("id", "cancelButton")
	    .attr("x", cancelButton_xPosition)
	    .attr("y", cancelButton_yPosition)
	    .attr("width", cancelButton_width)
	    .attr("height", cancelButton_height)
	    .attr("stroke", "red")
	    .attr("stroke-width", 3)
	    .style("fill-opacity", zeroOpacity)
	    .on("click", function(cancel){restoreTranscript();});
	
	d3.select("#editbox")
	    .append("g")
	    .selectAll("rect")
	    .data(json.exons).enter()
	    .append("rect")
	    .attr("id", function(d,i){return "exon_edit_" + i;})
	    .attr("x", function(exon) {
		return xScale(parseInt(exon.start));
	    })
	    .attr("y", function(exon, exon_index) {
		return getThisExonLevel(exonLevels, exon_index) * exon_editMode_height;
	    })
	    .attr("width", function(exon) {
		return xDiff( parseInt(exon.end) , parseInt(exon.start) );
	    })
	    .attr("height", exon_height)
	    .on("click", function(exon, exon_index){changeExon(exon, exon_index)});
    }
};

var deleteTranscript_startover = function(d, i, json){
    //console.log("this transcript is:" + i + " id:" + json.transcripts[i].id);
    json.transcripts.splice(i,1);
    startoverAlgo(json);
};

var addTranscript_startover = function(json){
    var newTranscript = {
	exons : [],
	id : "newTranscript"
    };
    json.transcripts.push(newTranscript);
    startoverAlgo(json);
};

var renderTranscripts = function(elem, json_raw){
    var json = jQuery.extend(true, {}, json_raw);

    var getEditRegion = function(d, i){
	showEditRegion(d, i, json);
    };

    var deleteTranscript = function(d, i){
	deleteTranscript_startover(d, i, json);
    };

    var addTranscript = function(){
	addTranscript_startover(json);
    };

    var switchTranscriptStatus = function(d, i){
	chooseTranscriptRegion(d, i, json, json_raw);
    };

    if(debug_flag){
	console.log(json);
    }

    var geneStart = parseInt(json.start);
    var geneEnd = parseInt(json.end);

    var svg = elem
    .attr("width", canvas_width)
    .attr("height",json.transcripts.length * 200 + 200)
    .append("svg")
    .attr("id", "geneGraph")
    .attr("width", elem.attr("width"))
    .attr("height", elem.attr("height"));

    var xScale = d3.scale.linear().domain([geneStart, geneEnd]).range([0, geneGraph_width]);
    var xDiff = function(pt1, pt2) { return Math.abs(xScale(pt1) - xScale(pt2)) }

    svg.selectAll("g")
    .data(json.transcripts).enter()
    .append("g")
    .attr("id", function(d, i){return "transcript_" + i;})
    .attr("transform", function(d,i) {
	    return "translate(0," + i*transcriptRegion_height + ")";
	})
    .selectAll("rect")
    .data(function(transcript){return transcript.exons;}).enter()
    .append("rect")
    .attr("class", "exon")
    .attr("id", function(exonNr, arrIndex){return "exon_" + exonNr;})
    .attr("x", function(exon) {return xScale(parseInt(json.exons[exon].start));})
    .attr("y", exon_yPosition)
    .attr("width", function(exon) {
	    return xDiff( parseInt(json.exons[exon].end) , parseInt(json.exons[exon].start) )
	})
    .attr("height", exon_height);

    svg.selectAll("g")
    .append("line")
    .attr("x1", xScale(geneStart))
    .attr("x2", xScale(geneEnd))
    .attr("y1", referenceLine_yPosition)
    .attr("y2", referenceLine_yPosition)
    .style("stroke", "#000");

    svg.selectAll("g")
    .append("rect")
    .attr("class", "area")
    .attr("id", function(d, i) {return "transcriptBackground_"+ i;})
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", geneGraph_width)
    .attr("height", transcriptRegion_height)
    .style("fill", "red")
    .style("fill-opacity", zeroOpacity)
    .on("mouseover", function(d,i){showTranscriptRegion(d, i)})
    .on("mouseout", function(d,i){hideTranscriptRegion(d, i)})
    .on("click", function(d, i){switchTranscriptStatus(d, i)});

    svg.selectAll("g")
    .append("text")
    .attr("x", editButton_textX)
    .attr("y", editButton_textY)
    .attr("font-size", 30)
    .text("Edit");

    svg.selectAll("g")
    .append("text")
    .attr("x", deleteButton_textX)
    .attr("y", deleteButton_textY)
    .attr("font-size", 30)
    .text("Delete");

    svg.selectAll("g")
    .append("rect")
    .attr("id", "editButton")
    .attr("x", editButton_xPosition)
    .attr("y", editButton_yPosition)
    .attr("width", editButton_width)
    .attr("height", editButton_height)
    .attr("stroke", "red")
    .attr("stroke-width", 3)
    .style("fill-opacity", zeroOpacity)
    .on("click", function(d, i){getEditRegion(d, i)});

    svg.selectAll("g")
    .append("rect")
    .attr("id", "deleteButton")
    .attr("x", deleteButton_xPosition)
    .attr("y", deleteButton_yPosition)
    .attr("width", deleteButton_width)
    .attr("height", deleteButton_height)
    .attr("stroke", "red")
    .attr("stroke-width", 3)
    .style("fill-opacity", zeroOpacity)
    .on("click", function(d, i){deleteTranscript(d, i)});

    
    svg.append("g")
    .attr("id", "addTranscript")
    .attr("transform", "translate(0," + json.transcripts.length*transcriptRegion_height + ")")
	.append("rect")
	.attr("id", "addButton")
	.attr("x", 0)
	.attr("y", 0)
	.attr("width", geneGraph_width)
	.attr("height", transcriptRegion_height)
	.attr("stroke", "red")
	.attr("stroke-width", 3)
	.style("fill-opacity", zeroOpacity)
	.on("click", addTranscript);

    svg.select("#addTranscript")
	.append("text")
	.attr("x", addTranscript_textX)
	.attr("y", addTranscript_textY)
	.attr("font-size", 55)
	.text("Add a new transcript");

};

var makeGeneGraph = function(json_raw){
	var elem = d3.select("#chart");
	elem.select("#geneGraph").remove();
	renderTranscripts(elem, json_raw);  
};

var startoverAlgo = function(json){
    for(exon_index =0; exon_index < json.exons.length; exon_index ++){
	if(json.exons[exon_index].clickState){
	    json.exons[exon_index].clickState = !json.exons[exon_index].clickState;
	}
    }  
    makeGeneGraph(json);
    getIsoformAbundance(json);
};