for i = 34:35
    audio = wavread(['C:\Temp\Track-' num2str(i,'%2.2d') '.wav']);
    audio2 = wavread(['C:\Temp\Tracks\Track' num2str(i,'%2.2d') '.wav']);
    dl(i) = length(audio) - length(audio2);
    e2(i) = any(audio2(:,2)-audio(7:(6+length(audio2)),2));
    e(i) = any(audio2(:,1)-audio(7:(6+length(audio2)),1));
end