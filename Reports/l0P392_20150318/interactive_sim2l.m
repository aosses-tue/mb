% Interactive Dynamic Field Simulator
% Sebastian Schneegans
% Ruhr-Universitaet Bochum, 2008


function interactive_sim2l

%%%%%%%%%%%%%%%%%%%%%%
% default parameters %
%%%%%%%%%%%%%%%%%%%%%%

fieldSize=181; % should be odd

tau_u = 20; tau_v = 5; % time constants
beta_u = 5; beta_v = 5; % steepness of sigmoid output functions
h_u = -5; h_v = -5; % resting levels of fields

c_uu = 20; sigma_uu = 5;
c_uv = 15; sigma_uv = 10; g_uv = 0;
c_vu = 15; sigma_vu = 5;

q_u = 0.0; q_v = 0.0;
sigma_q = 0;

iw1 = 5; ip1 = 60; ih1 = 0;
iw2 = 5; ip2 = 120; ih2 = 0;


%%%%%%%%%%%%%%%%%%%%
% program settings %
%%%%%%%%%%%%%%%%%%%%

% define how broad the interaction kernels should be
kernelWidthMultiplier = 3;

% define parameters that will be saved to or loaded from .mat file
saveParamNames = {'fieldSize','tau_u','tau_v','beta_u','beta_v','h_u','h_v',...
  'c_uu','c_uv','c_vu','g_uv','sigma_uu','sigma_uv','sigma_vu','q_u', 'q_v','sigma_q',...
  'iw1','ip1','ih1','iw2','ip2','ih2'};
saveParamStr = cellToString(saveParamNames);

paramFile = 0;
paramPath = 0;
advParamsChanged = false;
initialize = true;
loop = false;
sliderChanged = [];
quit = false;


%%%%%%%%%%%%%%
% create GUI %
%%%%%%%%%%%%%%

% note: some changes to the sliders must be done during the initialization
% (setting slider range for fieldSize-dependent values ip1 and ip2)

fig = figure('Position',[50,50,700,700],'Name','Neural Field',...
  'Color','w','NumberTitle','off','MenuBar','none');

% create axes for field plots
topAxes = axes('Position',[0.1 0.65 0.85 0.3]);
bottomAxes = axes('Position',[0.1 0.3 0.85 0.3]);


% create sliders for model parameters
controlFieldHeight = 0.045;
controlFieldWidth = 0.3;
sliderWidth = 0.2;
gapWidth = 0.01;
textWidth = controlFieldWidth - sliderWidth - gapWidth;


controlParamNames = {'h_u','h_v','c_uu','c_uv','c_vu','g_uv','q_u','q_v',...
  'iw1','ip1','ih1','iw2','ip2','ih2'};
controlPosX = [0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2] * controlFieldWidth;
controlPosY = [4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0] * controlFieldHeight;
controlMin = [-10, -10, 0, 0, 0, 0, 0, 0, 0.01, 0, 0, 0.01, 0, 0];
controlMax = [0, 0, 100, 100, 100, 5, 0.5, 0.5, 20, fieldSize-1, 20, 20, fieldSize-1, 20];
textFormat = {'%0.1f', '%0.1f', '%0.1f', '%0.1f', '%0.1f', '%0.1f', ...
  '%0.2f', '%0.2f', '%0.1f', '%0.0f', '%0.1f', '%0.1f', '%0.0f', '%0.1f'};

nControlParams = length(controlParamNames);
sliders = zeros(nControlParams, 1);
textFields = zeros(nControlParams, 1);
for i = 1 : nControlParams
  eval(['tmp = ' controlParamNames{i} ';']);
  sliders(i) = uicontrol(fig, 'Style', 'Slider', 'Units', 'Norm', 'Position', ...
    [controlPosX(i)+textWidth+gapWidth, controlPosY(i), sliderWidth, controlFieldHeight], ...
    'Value', tmp, 'Min', controlMin(i), 'Max', controlMax(i), 'Callback', @sliderCallback);
  textFields(i) = uicontrol(fig,'Style','Text','Units','Norm','HorizontalAlignment', 'left', ...
    'String',[controlParamNames{i} '=' num2str(tmp, textFormat{i})], ...
    'Position',[controlPosX(i)+gapWidth, controlPosY(i) textWidth controlFieldHeight]);
end


% create buttons for general functions
buttonLabels = {'Advanced', 'Reset', 'Load', 'Save', 'Quit'};
buttonToolTips = {'Advanced model parameters', 'Reset field activities', 'Load parameters', 'Save parameters', 'Quit'};
nButtons = length(buttonLabels);
buttons = zeros(nButtons, 1);
for i = 1 : nButtons
  buttons(i) = uicontrol(fig,'Style','Pushbutton','String',buttonLabels{i},...
    'Units','Norm','Position',...
    [3*controlFieldWidth, 5/5*(5-i)*controlFieldHeight, 1-3*controlFieldWidth, 5/5*controlFieldHeight],...
    'Callback',@buttonCallback,'ToolTip',buttonToolTips{i});
end


% prepare GUI for advanced model parameters
advButtonHeight = 0.25;
advControlFieldHeight = (1 - advButtonHeight) / 5;
advControlFieldWidth = 1/3;
advEditWidth = advControlFieldWidth / 2;
advTextWidth = advControlFieldWidth  - advEditWidth;

advControlParams = {'fieldSize', 'sigma_q', 'tau_u', 'tau_v', 'beta_u', 'beta_v', ...
  'sigma_uu', 'sigma_uv', 'sigma_vu'};
advControlPosX = [0 1 0 1 0 1 0 1 2] * advControlFieldWidth;
advControlPosY = [4 4 3 3 2 2 1 1 1] * advControlFieldHeight;
minSigma = 10^(-5);
advParamsMin = [1, 0, 1, 1, 0, 0, minSigma, minSigma, minSigma]; % minimum values for parameters

nAdvControlParams = length(advControlParams);
advEditFields = zeros(nAdvControlParams);
advTextFields = zeros(nAdvControlParams);
advFig = 0;
advCancelButton = 0;
advOkButton = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop with inialization and simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while ~quit
  
  % load parameters from file
  if paramFile
    fileParamNames = who('-file', [paramPath paramFile]);
    missingParamNames = setdiff(saveParamNames, fileParamNames);
    if ~isempty(missingParamNames)
      msgbox(['The follwing model parameters were not found in the file: ' ...
        cellToString(missingParamNames) ' Previously assigned values will be used for these parameters.'],...
        'Warning', 'warn');
    end
    loadParamNames = intersect(saveParamNames, fileParamNames);
    if ~isempty(loadParamNames)
      eval(['load ''' paramPath paramFile ''' ' cellToString(loadParamNames)]);
    end
    paramFile = 0;
    initialize = 1;
  end
  
  
  % change advanced parameters
  if advParamsChanged
    paramsOk = true;
    for i = 1 : nAdvControlParams
      tmp = str2double(get(advEditFields(i), 'String'));
      if isnan(tmp) || tmp < advParamsMin(i)
        msgbox(['Value for parameter ' advControlParams{i} ' is outside the valid range.'],...
          'Error', 'error');
        paramsOk = false;
      end
    end
    if paramsOk
      for i = 1 : nAdvControlParams
        tmp = str2double(get(advEditFields(i), 'String'));
        eval([advControlParams{i} ' = tmp;']);
      end
      close(advFig);
      initialize = true;
    end
    advParamsChanged = false;
  end
 
  
  % model initialization
  if initialize
    
    % make sure that input positions are within valid range
    ip1 = min(ip1, fieldSize-1);
    ip2 = min(ip2, fieldSize-1);
    
    % set values of sliders and text fields
    for i = 1 : nControlParams
      if strcmp(controlParamNames{i}, 'ip1') || strcmp(controlParamNames{i}, 'ip2')
        set(sliders(i), 'Max', fieldSize-1);
      end
      
      eval(['tmp = ' controlParamNames{i} ';']);
      if get(sliders(i), 'Min') > tmp
        set(sliders(i), 'Min', floor(tmp));
      elseif get(sliders(i), 'Max') < tmp
        set(sliders(i), 'Max', ceil(tmp));
      end
      set(sliders(i), 'Value', tmp);
      set(textFields(i), 'String', [controlParamNames{i} '=' num2str(tmp, textFormat{i})]);
    end


    % initialization of dynamic fields
    field_u = zeros(1, fieldSize) + h_u;
    field_v = zeros(1, fieldSize) + h_v;
    output_u = sigmoid(field_u, beta_u, 0.0);
    output_v = sigmoid(field_v, beta_v, 0.0);

    input = circshift(ih1 * exp(-0.5*((1:fieldSize) - floor(fieldSize/2)).^2 / iw1.^2), ...
      [1; round(ip1) - floor(fieldSize/2) + 1]) + ...
      circshift(ih2 * exp(-0.5*((1:fieldSize) - floor(fieldSize/2)).^2 / iw2.^2), ...
      [1; round(ip2) - floor(fieldSize/2) + 1]);

    
    % create interaction kernels
    kSize_uu = min(round(kernelWidthMultiplier * sigma_uu), floor((fieldSize-1)/2));
    kSize_uv = min(round(kernelWidthMultiplier * sigma_uv), floor((fieldSize-1)/2));
    kSize_vu = min(round(kernelWidthMultiplier * sigma_vu), floor((fieldSize-1)/2));

    extIndex_uu = [fieldSize-kSize_uu+1:fieldSize, 1:fieldSize, 1:kSize_uu];
    extIndex_uv = [fieldSize-kSize_uv+1:fieldSize, 1:fieldSize, 1:kSize_uv];
    extIndex_vu = [fieldSize-kSize_vu+1:fieldSize, 1:fieldSize, 1:kSize_vu];

    kernel_uu = gaussNorm(-kSize_uu : kSize_uu, 0, sigma_uu);
    kernel_uv = gaussNorm(-kSize_uv : kSize_uv, 0, sigma_uv);
    kernel_vu = gaussNorm(-kSize_vu : kSize_vu, 0, sigma_vu);
    
    kSize_q = min(round(kernelWidthMultiplier * sigma_q), floor((fieldSize-1)/2));
    extIndex_q = [fieldSize-kSize_q+1:fieldSize, 1:fieldSize, 1:kSize_q];
    kernel_q = gaussNorm(-kSize_q : kSize_q, 0, sigma_q);
    
    % plot graphs
    axes(topAxes);
    cla;
    hold on;
    plot([0,fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    actPlot_u = plot(0:fieldSize-1,field_u,'color','b','Linewidth',3);
    outPlot_u = plot(0:fieldSize-1,10*output_u,'color','r','Linewidth',1);
    inPlot = plot(0:fieldSize-1,input+h_u,'color','g','Linewidth',1);
    set(gca,'ylim',[-15,15],'xlim',[0,fieldSize-1],'Ytick',[-10,0,10]);
    ylabel('u field','Fontsize',12);
    hold off;

    axes(bottomAxes);
    cla;
    hold on;
    plot([0,fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    actPlot_v = plot(0:fieldSize-1,field_v,'color','b','Linewidth',3);
    outPlot_v = plot(0:fieldSize-1,10*output_v,'color','r','Linewidth',1);
    set(gca,'ylim',[-15,15],'xlim',[0,fieldSize-1],'Ytick',[-10,0,10]);
    ylabel('v field','Fontsize',12);
    hold off;

    initialize = false;
    loop = true;
  end
  
  pause(0.1);
  
  % loop over time steps
  while loop
    output_u = sigmoid(field_u, beta_u, 0.0);
    output_v = sigmoid(field_v, beta_v, 0.0);
    totalOutput_v = sum(output_v);

    fieldConv_uu = conv2(1, kernel_uu, output_u(extIndex_uu), 'valid');
    fieldConv_uv = conv2(1, kernel_uv, output_v(extIndex_uv), 'valid');
    fieldConv_vu = conv2(1, kernel_vu, output_u(extIndex_vu), 'valid');

    % prepapre field noise (with spatial correlation if required)
    noise_u = q_u * randn(1, fieldSize);
    noise_v = q_v * randn(1, fieldSize);
    if sigma_q > 0
      noise_u = conv2(1, kernel_q, noise_u(extIndex_q), 'valid');
      noise_v = conv2(1, kernel_q, noise_v(extIndex_q), 'valid');
    end
    
    field_u = field_u + 1/tau_u * ( -field_u + h_u + input ...
      + c_uu*fieldConv_uu - c_uv*fieldConv_uv - g_uv*totalOutput_v) ...
      + noise_u;
    field_v = field_v + 1/tau_v * ( -field_v + h_v ...
      + c_vu*fieldConv_vu) + noise_v;
    
    set(actPlot_u,'Ydata',field_u);
    set(outPlot_u,'Ydata',10*output_u);
    set(actPlot_v,'Ydata',field_v);
    set(outPlot_v,'Ydata',10*output_v);


    drawnow;
    pause(0.01);
  end
end
close(fig);


%%%%%%%%%%%%%%%%%%%%%%
% callback functions %
%%%%%%%%%%%%%%%%%%%%%%

function buttonCallback(hObject, eventdata)
  buttonPressed = find(hObject == buttons);

  switch buttonPressed
    case {1} % Advanced
      loop = false;
      for j = 1 : nControlParams
        set(sliders(j), 'Enable', 'off');
      end
      for j = 1 : nButtons
        set(buttons(j), 'Enable', 'off');
      end
     
      advFig = figure('Position',[50,50,400,200],'Name','Advanced model parameters',...
        'Color','w','NumberTitle','off','MenuBar','none','DeleteFcn',@advCloseCallback);
      for j = 1 : nAdvControlParams
        eval(['tmp = ' advControlParams{j} ';']);
        advEditFields(j) = uicontrol(advFig, 'Style', 'Edit', 'Units', 'Norm', 'Position', ...
          [advControlPosX(j)+advTextWidth, advControlPosY(j)+advButtonHeight, advEditWidth, advControlFieldHeight], ...
          'String',  num2str(tmp), 'Min', 0, 'Max', 0);
        advTextFields(j) = uicontrol(advFig,'Style','Text','Units','Norm','HorizontalAlignment', 'left', ...
          'String',[advControlParams{j} ':'], ...
          'Position',[advControlPosX(j), advControlPosY(j)+advButtonHeight, advTextWidth, advControlFieldHeight]);
      end
      advCancelButton = uicontrol(advFig, 'Style', 'Pushbutton', 'String', 'Cancel', 'Units', 'Norm', ...
        'Position', [0.25, 0, 0.25, advButtonHeight], 'Callback', @advCancelCallback);
      advOkButton = uicontrol(advFig, 'Style', 'Pushbutton', 'String', 'Ok', 'Units', 'Norm', ...
        'Position', [0.5, 0, 0.25, advButtonHeight], 'Callback', @advOkCallback);  
    case {2} % Reset
      field_u = zeros(1, fieldSize) + h_u;
      field_v = zeros(1, fieldSize) + h_v;
    case {3} % Load
      [paramFile, paramPath] = uigetfile('*.mat', 'Load parameter file');
      loop = (paramFile == 0);
    case {4} % Save
      [paramFile, paramPath] = uiputfile('*.mat', 'Save parameter file');
      if paramFile
        eval(['save ''' paramPath paramFile ''' ' saveParamStr ' -v6']);
        paramFile = 0;
      end
    case {5} % Quit
      quit = true;
      loop = false;
  end
end


% update paramter values after slider changed
function sliderCallback(hObject, eventdata)
  sliderChanged = find(hObject == sliders);

  paramName = controlParamNames{sliderChanged};
  tmp = get(sliders(sliderChanged), 'Value');
  set(textFields(sliderChanged), 'String', [paramName '=' num2str(tmp, textFormat{sliderChanged})]);
  eval([paramName '= tmp;']);
  if paramName(1) == 'i' || strcmp(paramName, 'h_u')
    input = circshift(ih1 * exp(-0.5*((1:fieldSize) - floor(fieldSize/2)).^2 / iw1.^2), ...
      [1; round(ip1) - floor(fieldSize/2) + 1]) + ...
      circshift(ih2 * exp(-0.5*((1:fieldSize) - floor(fieldSize/2)).^2 / iw2.^2), ...
      [1; round(ip2) - floor(fieldSize/2) + 1]);
    set(inPlot, 'Ydata',input+h_u);
  end
end


% close advanced parameters windows after cancel
function advCancelCallback(hObject, eventdata)
  close(advFig);
end


% continue simulation after advanced parameter setting has been cancelled
function advCloseCallback(hObject, eventdata)
  for j = 1 : nControlParams
    set(sliders(j), 'Enable', 'on');
  end
  for j = 1 : nButtons
    set(buttons(j), 'Enable', 'on');
  end

  loop = true;
end



% update parameters after advanced parameter setting has been confirmed
function advOkCallback(hObject, eventdata)
  advParamsChanged = true;
end

end


%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary functions %
%%%%%%%%%%%%%%%%%%%%%%%

% create single string from cell array of strings
% elements of cell array are seperated by blanks
function str = cellToString(c)
  str = '';
  for i = 1 : length(c)
    str = [str c{i} ' '];
  end
  deblank(str);
end


