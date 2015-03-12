function getStimuli () {
    misc.checkParams( Array("path", "prefix"));
    files = api.files( params["path"] );
    
    r="";
    for(i=0; i<files.length; ++i) {
        print( files[i] );
        r=r+xml.stimulus( params["prefix"]+files[i], "datablock_" + files[i]); 
    }
    
    return r;
}


function getDatablocks() {
    misc.checkParams( Array("path", "device"));
    files = api.files( params["path"]);
    path=api.stripPath(params["path"]);
    
    r="";
    for(i=0; i<files.length; ++i) {
        print( files[i] );
        r=r+xml.datablock( "datablock_"+files[i], params["device"], path+files[i]); 
    }
    
    return r;
}


function getTrials() {
    misc.checkParams( Array("path", "screen", "answer"));
    files = api.files( params["path"]);
    path=api.stripPath(params["path"]);
    
    r="";
    for(i=0; i<files.length; ++i) {
        print( files[i] );
        r=r+xml.tag("trial", {id: "trial_"+files[i]}, 
            xml.tag("answer", params["answer"]),
            xml.tag("screen", {id: params["screen"]} ) ,
            xml.tag("stimulus", {id: "stimulus_"+files[i]})); 
    }
    
    return r;
}

