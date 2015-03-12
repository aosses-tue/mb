function vocoder_makeexperiment
% Copyright Tom Francart, 2009

sourcepath='/mnt/l/speechmaterials/english/bkb/';
targetpath_base='/home/tom/temp/audis-cisim/';


for list=1:2
    targetpath=[targetpath_base 'list' num2str(list) '/'];
    if (~exist(targetpath))
       mkdir(targetpath, '.'); 
    end
    copyfile('sentencetest.apx', targetpath);
    copyfile('sentencetest.js', targetpath);
    for sentence=1:16
        nr=num2str(100*list+sentence);
        if (nr<999)
            nr=['0' nr];
        end
        sourcefile=sprintf('%s/bkbe%s.wav', sourcepath, nr);
        [d,fs]=wavread(sourcefile);
        r=vocoder(d(:,1),fs);
        wavwrite(r/max(abs(r)), fs, [targetpath 'sentence' nr '.wav']);
    end
end


