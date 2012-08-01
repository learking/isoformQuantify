var gene_mgi_id;
var fileSelected;

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

var calculateIsoforms = function (form) {
    var userQuery = form.geneSearchBox.value.match(/\w+/)[0];
    fileSelected = form.fileSelection.value;
    console.log(fileSelected);
    jsonp("http://cbfg-dev/gene_name_lookup.php?genes=" + userQuery + "&callback=getGeneHits");
};

var runScript = function(e, form){
    if (e.keyCode == 13) {
 	e.preventDefault();
	calculateIsoforms(form);
    }
};

var getIsoformAbundance = function(json){
    json['bamFile'] = fileSelected;
    $.ajax({  
	    url: "../getIsoformAbund",  
	    type: "POST",  
	    dataType: "json",
	    data: JSON.stringify(json),
	    success: function(data){
		showIsoformAbund(data);
		console.log(data);              
	    },  
	    error: function(){  
		console.log("post fail :-(");  
	    }  
	});
};

var showIsoformAbund = function(data){
    d3.selectAll("#transcript_abund").remove();
    for (var key in data){
	var currentAbund = Math.round(parseFloat(data[key])*1000)/1000;
	d3.select("#" + key).append("text")
	.attr("id", "transcript_abund")
	.attr("x", 1220)
	.attr("y", 60)
	.attr("font-size", 40)
	.text(currentAbund.toString());
    }
};