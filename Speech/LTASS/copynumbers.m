function copynumbers(a, speaker)
% a is matrix van nummers (lijst, item)
% speaker is speaker code

sourcedir='l:\speechmaterials\numbers\';
targetdir='c:\tom\woordenlijst\numbers\';


for lijst=1:size(a,1)
    disp(['mkdir ' targetdir speaker '\' num2str(lijst)]);

    for file=1:size(a,2)
        disp(['copy ' sourcedir speaker num2str( a(lijst,file)) '.wav ' targetdir speaker '\' num2str(lijst)]);
    end

end