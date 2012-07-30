var calculateIsoforms = function (form) {
    var geneID = form.geneSearchBox.value.match(/\w+/)[0];
    var fileSelected = form.fileSelection.value;
    var geneJSON;
    var isoformAbundanceJSON;

    d3.json("gene_" + geneID + ".json", function(json) {
	try{
		geneJSON = json;
		makeGeneGraph(geneJSON);
		getIsoformAbundance(geneJSON);
	}
	catch(error){
	    alert("this gene does not exist");
	}
    });  

};

var runScript = function(e, form){
    if (e.keyCode == 13) {
 	e.preventDefault();
	calculateIsoforms(form);
    }
};

var getIsoformAbundance = function(geneJSON){
    $.ajax({  
	    url: "../getIsoformAbund",  
	    type: "POST",  
	    dataType: "json",
	    data: JSON.stringify(geneJSON),
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