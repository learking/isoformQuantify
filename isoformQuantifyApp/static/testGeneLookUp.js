var globalJSON;

var myJsonpCallback =function (json){
    globalJSON = json;
    console.log(globalJSON);
    writeJSON(globalJSON);
}

var jsonp = function(url) {
    var head = document.head;
    var script = document.createElement("script");

    script.setAttribute("src", url);
    head.appendChild(script);
    head.removeChild(script);
}

var calculateIsoforms = function(form){
    var geneID = form.geneSearchBox.value.match(/\w+/)[0];
    jsonp("http://cbfg-dev/gene_name_lookup.php?genes=" + geneID + "&callback=myJsonpCallback");
}

var writeJSON = function(globalJSON){
    console.log(globalJSON);
}