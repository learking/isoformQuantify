var gene_mgi_id;

var jsonp = function(url) {
    var head = document.head;
    var script = document.createElement("script");
    script.setAttribute("src", url);
    head.appendChild(script);
    head.removeChild(script);
};

var getGeneHits = function(geneHits){
// no hit

// multiple hit

// one hit, but not within dataset

    numHits = Object.keys(geneHits).length;
    for (var mgi_id in geneHits){
	gene_mgi_id = mgi_id;
	break;
    }
    console.log(gene_mgi_id);
    console.log(mgi2ensembl[gene_mgi_id]);

    var gene_ensembl_id = mgi2ensembl[gene_mgi_id];

    d3.json("gene_" + gene_ensembl_id + ".json", function(json) {
	try{
		makeGeneGraph(json);
		getIsoformAbundance(json);
	}
	catch(error){
	    alert("this gene does not exist");
	}
    });  

};

var calculateIsoforms = function (form) {
    var userQuery = form.geneSearchBox.value.match(/\w+/)[0];
    var fileSelected = form.fileSelection.value;
    jsonp("http://cbfg-dev/gene_name_lookup.php?genes=" + userQuery + "&callback=getGeneHits");
};

var runScript = function(e, form){
    if (e.keyCode == 13) {
 	e.preventDefault();
	calculateIsoforms(form);
    }
};

var getIsoformAbundance = function(json){
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