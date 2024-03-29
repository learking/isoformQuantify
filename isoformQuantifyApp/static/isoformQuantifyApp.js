var debug_flag = false;
var gene_mgi_id;
var arSelected = new Array();
var userQuery;

var jsonp = function(url) {
    var head = document.head;
    var script = document.createElement("script");
    script.setAttribute("src", url);
    head.appendChild(script);
    head.removeChild(script);
};

var showGeneID = function(geneInfo){
    var elem = d3.select("#chart");
    elem.select("#geneInfo").remove();
    
    var geneInfo_svg = elem
	.append("svg")
	.attr("id", "geneInfo")
	.attr("width", canvas_width)
	.attr("height", 60)
	.append("text")
	.attr("x", 0)
	.attr("y", 40)
	.attr("font-size", 40)
	.text(geneInfo);
};

var getGeneHits = function(geneHits){
    
    numHits = Object.keys(geneHits).length;
    // multiple hits
    if (numHits > 1){
	alert("Warning: search result is not unique");
    }

    for (var mgi_id in geneHits){
	gene_mgi_id = mgi_id;
	break;
    }

    // no hit
    if(gene_mgi_id == "Not found"){
	alert("Gene lookup service won't be able to locate this gene!");
    }
    else{
	var geneInfo = "Gene Symbol:" + geneHits[gene_mgi_id]["symbol"] + " || " + "MGI ID:" + gene_mgi_id + " || " + "ENSEMBL ID:" +  mgi2ensembl[gene_mgi_id];

	showGeneID(geneInfo);

	var gene_ensembl_id = mgi2ensembl[gene_mgi_id];

	d3.json("gene_" + gene_ensembl_id + ".json", function(json) {
	    try{
		makeGeneGraph(json);
		if(json.transcripts.length > 1){
		    getIsoformAbundance(json);
		}
	    }
	    catch(error){
		// one hit, but not within dataset
		alert("this gene does not exist");
	    }
	});  
    }

};

var selectDataset = function(select){
    var option = select.options[select.selectedIndex];
    if(!option){
	alert("must select at least one dataset!");
    }
    else{
	var ul = select.parentNode.getElementsByTagName('ul')[0];
	
	var choices = ul.getElementsByTagName('input');
	for (var i = 0; i < choices.length; i++)
	    if (choices[i].value == option.value)
		return;
	
	var li = document.createElement('li');
	var input = document.createElement('input');
	var text = document.createTextNode(option.firstChild.data);
	
	input.type = 'hidden';
	input.name = 'ingredients[]';
	input.value = option.value;

	li.appendChild(input);
	li.appendChild(text);
	li.setAttribute('onclick', 'removeOption(this);');     
	
	ul.appendChild(li);

	var options = document.all.tags("option");

	while (select.selectedIndex != -1)
	{
	    if (arSelected.indexOf(select.options[select.selectedIndex].value) == -1){
		arSelected.push(select.options[select.selectedIndex].value);
	    }
            select.options[select.selectedIndex].selected = false;
	}
	// You can use the arSelected array for further processing.                 
    }
};

var removeOption = function(option){
    arSelected.splice(arSelected.indexOf(option.childNodes[0].value), 1);
    option.parentNode.removeChild(option);
    jsonp("http://cbfg-dev/gene_name_lookup.php?genes=" + userQuery + "&callback=getGeneHits");
};

var calculateIsoforms = function (form) {
    if(form.geneSearchBox.value){
	userQuery = form.geneSearchBox.value.match(/\w+/)[0];
	selectDataset(form.fileSelection);
	jsonp("http://cbfg-dev/gene_name_lookup.php?genes=" + userQuery + "&callback=getGeneHits");
    }
    else{
	alert("must input a gene to start with!");
    }
};

var runScript = function(e, form){
    if (e.keyCode == 13) {
 	e.preventDefault();
	calculateIsoforms(form);
    }
};

var getIsoformAbundance = function(json){
    json['bamFiles'] = arSelected;
    $.ajax({  
	    url: "../getIsoformAbund",  
	    type: "POST",  
	    dataType: "json",
	    data: JSON.stringify(json),
	    success: function(data){
		showIsoformAbund(data);
		if(debug_flag){
		    console.log(data);              
		}
	    },  
	    error: function(){  
		console.log("post fail :-(");  
	    }  
	});
};

var showIsoformAbund = function(data){
    d3.selectAll("#transcript_abund").remove();
    d3.select("#transcript_abund_background").remove();

    d3.select("#geneGraph")
	.append("g")
	.attr("id", "transcript_abund_background");

    var trans_abund_background = d3.select("#transcript_abund_background");
    var trans_abund_background_height = Object.keys(data[0]).length * 100;
    for (var i = 0; i < arSelected.length; i++){
	trans_abund_background
	    .append("rect")
	    .attr("class", arSelected[i].split(".")[0])
	    .attr("x", 1385 + i*130)
	    .attr("y", 0)
	    .attr("width", 115)
	    .attr("height", trans_abund_background_height)
	    .style("fill-opacity", zeroOpacity);
    }

    for (var i = 0; i < data.length; i++){
	for (var key in data[i]){
	    var currentAbund = Math.round(parseFloat(data[i][key])*1000)/1000;
	    d3.select("#" + key).append("text")
		.attr("id", "transcript_abund")
		.attr("x", 1400 + i*130)
		.attr("y", 65)
		.attr("font-size", 40)
		.text(currentAbund.toString());
	}
    }
};