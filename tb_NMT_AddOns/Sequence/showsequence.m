function showsequence(qic,total_rate,handle)

if (nargin<2)
    total_rate=NaN;
end
if (nargin<3)
    handle=gca;
end

if (isfield(qic,'periods'))
    if (length(qic.periods)>1)
        time=cumsum(qic.periods)/1e6;
        time=[0; time(1:end-1)];
    else
        time=[0:length(qic.channels)-1]*qic.periods/1e6;
        if (~isnan(total_rate) && 1e6/total_rate ~= qic.periods)
            warning('Rate mismatch');
        end
    end
else
    if (nargin<2 || isnan(total_rate))
        error('Total rate not specified');
    end
    time=[0:length(qic.channels)-1]/total_rate;
end
hold(handle, 'on');

if (~isfield(qic,'channels') && isfield(qic,'electrodes'))
    qic.channels=-qic.electrodes+23
    qic.channels(qic.channels==23)=0;
end
if (~isfield(qic,'magnitudes') && isfield(qic,'current_levels'))
    qic.magnitudes=qic.current_levels/255;
end
    

% plot(time(odd), qic.channels(odd)+qic.magnitudes(odd), '.');
% plot(time(even), qic.channels(even)+qic.magnitudes(even), '.r');
% odd=find( mod(qic.channels,2));
% even=find( mod(qic.channels,2)==0);

for c=1:22
    I=find(qic.channels==c & qic.magnitudes>=0);
    if (~isempty(I))
        if (mod(c,2))
            color='r';
        else
            color='b';
        end
        h=stem(handle, time(I), qic.magnitudes(I)+c, color, 'Basevalue', c, ...
            'Marker', 'none', 'LineWidth', 1);
        delete(get(h,'baseline'));
    end
    
    % emphasize T level stimulation
    I=find(qic.channels==c & qic.magnitudes==0);
    if (~isempty(I))
        if (mod(c,2))
            color='r';
        else
            color='b';
        end
        h=plot(handle,time(I), qic.magnitudes(I)+c, color, 'Marker', '.', 'LineStyle', 'none');
    end
end

xlabel(handle, 'Time (s)');
ylabel(handle, {'Magnitude' '+ #electrode'});

% draw lines between electrodes
%maxx=max(length(wav)/fs, length(qic.channels)/total_rate);
maxx=time(end);
for i=1:22
   plot(handle,[0 maxx], [i i], ':', 'Color', [1 1 1]*0.5, 'LineWidth', 0.1);
end
ylim(handle, [1 23]);