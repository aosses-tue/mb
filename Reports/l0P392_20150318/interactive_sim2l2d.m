% Interactive Dynamic Field Simulator
% Sebastian Schneegans
% Ruhr-Universitaet Bochum, 2008


function interactive_sim2l2d

%%%%%%%%%%%%%%%%%%%%%%
% default parameters %
%%%%%%%%%%%%%%%%%%%%%%

fieldSize=101; % should be odd

tau_u = 20; tau_v = 5; % time constants
h_u = -5; h_v = -5; % resting levels of fields
beta_u = 5; beta_v = 5; % steepness of sigmoid output functions

c_uu = 25; sigma_uu = 5;
c_uv = 25; sigma_uv = 10; g_uv = 0;
c_vu = 20; sigma_vu = 5;

q_u = 0.0; q_v = 0.0;
sigma_q = 0;

iw1 = 4; ip1x = 25; ip1y = 25; ih1 = 0;
iw2 = 4; ip2x = 75; ip2y = 75; ih2 = 0;
iw3x = 4; ip3x = 25; ih3 = 0;
iw4y = 4; ip4y = 25; ih4 = 0;


%%%%%%%%%%%%%%%%%%%%
% program settings %
%%%%%%%%%%%%%%%%%%%%

% define how broad the interaction kernels should be
kernelWidthMultiplier = 3;

% define parameters that will be saved to or loaded from .mat file
saveParamNames = {'fieldSize','tau_u','tau_v','beta_u','beta_v','h_u','h_v',...
  'c_uu','c_uv','c_vu','g_uv','sigma_uu','sigma_uv','sigma_vu','q_u','q_v','sigma_q',...
  'iw1','ip1x','ip1y','ih1','iw2','ip2x','ip2y','ih2','iw3x','ip3x','ih3','iw4y','ip4y','ih4'};
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
% (setting slider range for fieldSize-dependent values ip*)

fig = figure('Position',[50,50,900,600],'Name','Neural Field',...
  'Color','w','NumberTitle','off','MenuBar','none');

% create axes for field plots
leftAxes = axes('Position',[0.08 0.4 0.36 0.54]);
rightAxes = axes('Position',[0.56 0.4 0.36 0.54]);
colorAxes = axes('Position',[0.3 0.32 0.4 0.02]);
imagesc([-15, 15], [0, 0], linspace(-15, 15, 256), [-15 15]);
set(colorAxes, 'YTick', [], 'XTick', [-15 -10 -5 0 5 10 15]);
colormap(jet(256));


% create sliders for model parameters
controlFieldHeight = 0.04;
controlFieldWidth = 0.22;
sliderWidth = 0.15;
gapWidth = 0.005;
textWidth = controlFieldWidth - sliderWidth - gapWidth;

controlParamNames = {'h_u','h_v','c_uu','c_uv','c_vu','g_uv','q_u','q_v',...
  'iw1','ip1x','ip1y','ih1','iw2','ip2x','ip2y','ih2',...
  'iw3x','ip3x','ih3','iw4y','ip4y','ih4'};
controlPosX = [0, 1, 0, 1, 2, 3, 0, 1, ...
  0, 1, 2, 3, 0, 1, 2, 3, ...
  0, 1, 3, 0, 1, 3] * controlFieldWidth;
controlPosY = [6, 6, 5, 5, 5, 5, 4, 4, ...
  3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0] * controlFieldHeight;
controlMin = [-10, -10, 0, 0, 0, 0, 0, 0, ...
  0.01, 0, 0, 0, 0.01, 0, 0, 0, ...
  0.01, 0, 0, 0.01, 0, 0];
controlMax = [0, 0, 100, 100, 100, 1, 0.5, 0.5, ...
  20, fieldSize-1, fieldSize-1, 20, 20, fieldSize-1, fieldSize-1, 20 ...  
  20, fieldSize-1, 20, 20, fieldSize-1, 20];
textFormat = {'%0.1f', '%0.1f', '%0.1f', '%0.1f', '%0.1f', '%0.2f', '%0.2f', '%0.2f',...
  '%0.1f', '%0.0f', '%0.0f', '%0.1f', '%0.1f', '%0.0f', '%0.0f', '%0.1f', ...
  '%0.1f', '%0.0f', '%0.1f', '%0.1f', '%0.0f', '%0.1f'};

nControlParams = length(controlParamNames);
sliders = zeros(nControlParams, 1);
textFields = zeros(nControlParams, 1);
for i = 1 : nControlParams
  eval(['tmp = ' controlParamNames{i} ';']);
  sliders(i) = uicontrol(fig, 'Style', 'Slider', 'Units', 'Norm', 'Position', ...
    [controlPosX(i)+textWidth+gapWidth, controlPosY(i), sliderWidth, controlFieldHeight], ...
    'Value', tmp, 'Min', controlMin(i), 'Max', controlMax(i), 'Callback', @sliderCallback);
  textFields(i) = uicontrol(fig,'Style','Text','Units','Norm','HorizontalAlignment', 'left', ...
    'String',[controlParamNames{i} '=' num2str(tmp, textFormat{i})], 'BackgroundColor', 'w', ...
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
    [4*controlFieldWidth, 7/5*(5-i)*controlFieldHeight, 1-4*controlFieldWidth, 7/5*controlFieldHeight],...
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
    halfSize = floor(fieldSize/2);
    
    % make sure that input positions are within valid range
    ip1x = min(ip1x, fieldSize-1);
    ip1y = min(ip1y, fieldSize-1);
    ip2x = min(ip2x, fieldSize-1);
    ip2y = min(ip2y, fieldSize-1);
    ip3x = min(ip3x, fieldSize-1);
    ip4y = min(ip4y, fieldSize-1);
    
    % set values of sliders and text fields
    for i = 1 : nControlParams
      if strcmp(controlParamNames{i}, 'ip1x') || strcmp(controlParamNames{i}, 'ip1y') ...
          || strcmp(controlParamNames{i}, 'ip2x') || strcmp(controlParamNames{i}, 'ip2y')
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
    field_u = zeros(fieldSize, fieldSize) + h_u;
    field_v = zeros(fieldSize, fieldSize) + h_v;
    output_u = sigmoid(field_u, beta_u, 0.0);
    output_v = sigmoid(field_v, beta_v, 0.0);

    
    input = circshift(ih1 * gauss2d(1:fieldSize, 1:fieldSize, halfSize+1, halfSize+1, iw1, iw1), ...
      [round(ip1x) - halfSize; round(ip1y) - halfSize]) + ...
      circshift(ih2 * gauss2d(1:fieldSize, 1:fieldSize, halfSize+1, halfSize+1, iw2, iw2), ...
      [round(ip2x) - halfSize; round(ip2y) - halfSize]) + ...
      circshift(ih3 * gauss2d(1:fieldSize, 1:fieldSize, halfSize+1, halfSize+1, iw3x, inf), ...
      [round(ip3x) - halfSize; 0]) + ...
      circshift(ih4 * gauss2d(1:fieldSize, 1:fieldSize, halfSize+1, halfSize+1, inf, iw4y), ...
      [0; round(ip4y) - halfSize]);

    
    % create interaction kernels
    kernelRange_uu = min(round(kernelWidthMultiplier * sigma_uu), floor((fieldSize-1)/2));
    kernelRange_uv = min(round(kernelWidthMultiplier * sigma_uv), floor((fieldSize-1)/2));
    kernelRange_vu = min(round(kernelWidthMultiplier * sigma_vu), floor((fieldSize-1)/2));

    extIndex_uu = [fieldSize-kernelRange_uu+1:fieldSize, 1:fieldSize, 1:kernelRange_uu];
    extIndex_uv = [fieldSize-kernelRange_uv+1:fieldSize, 1:fieldSize, 1:kernelRange_uv];
    extIndex_vu = [fieldSize-kernelRange_vu+1:fieldSize, 1:fieldSize, 1:kernelRange_vu];

    kernel_uu = gaussNorm(-kernelRange_uu : kernelRange_uu, 0, sigma_uu);
    kernel_uv = gaussNorm(-kernelRange_uv : kernelRange_uv, 0, sigma_uv);
    kernel_vu = gaussNorm(-kernelRange_vu : kernelRange_vu, 0, sigma_vu);
    
    kernelRange_q = min(round(kernelWidthMultiplier * sigma_q), floor((fieldSize-1)/2));
    extIndex_q = [fieldSize-kernelRange_q+1:fieldSize, 1:fieldSize, 1:kernelRange_q];
    kernel_q = gaussNorm(-kernelRange_q : kernelRange_q, 0, sigma_q);
    
    % plot graphs
    axes(leftAxes);
    cla;
    actPlot_u = imagesc(0:fieldSize-1, 0:fieldSize-1, field_u, [-15, 15]);
    title('u field','Fontsize',12);

    axes(rightAxes);
    cla;
    actPlot_v = imagesc(0:fieldSize-1, 0:fieldSize-1, field_v, [-15, 15]);
    title('v field','Fontsize',12);

    initialize = false;
    loop = true;
  end
  
  pause(0.1);
  
  % loop over time steps
  while loop
    output_u = sigmoid(field_u, beta_u, 0.0);
    output_v = sigmoid(field_v, beta_v, 0.0);
    totalOutput_v = sum(sum(output_v));

    % perform convolutions of padded fields in both dimensions
    conv_uu = conv2(kernel_uu', 1, output_u(extIndex_uu, :), 'valid');
    conv_uu = conv2(1, kernel_uu, conv_uu(:, extIndex_uu), 'valid');
    
    conv_vu = conv2(kernel_vu', 1, output_u(extIndex_vu, :), 'valid');
    conv_vu = conv2(1, kernel_vu, conv_vu(:, extIndex_vu), 'valid');
    
    conv_uv = conv2(kernel_uv', 1, output_v(extIndex_uv, :), 'valid');
    conv_uv = conv2(1, kernel_uv, conv_uv(:, extIndex_uv), 'valid');

    % prepapre field noise (with spatial correlation if required)
    noise_u = q_u * randn(fieldSize, fieldSize);
    noise_v = q_v * randn(fieldSize, fieldSize);
    if sigma_q > 0
      noise_u = conv2(kernel_q', 1, noise_u(extIndex_q, :), 'valid');
      noise_u = conv2(1, kernel_q, noise_u(:, extIndex_q), 'valid');
      noise_v = conv2(kernel_q', 1, noise_v(extIndex_q, :), 'valid');
      noise_v = conv2(1, kernel_q, noise_v(:, extIndex_q), 'valid');
    end
        
    field_u = field_u + 1/tau_u * ( -field_u + h_u + input ...
      + c_uu*conv_uu - c_uv*conv_uv - g_uv*totalOutput_v) ...
      + noise_u;
    field_v = field_v + 1/tau_v * ( -field_v + h_v ...
      + c_vu*conv_vu) + noise_v;
    
    set(actPlot_u,'CData',field_u);
    set(actPlot_v,'CData',field_v);

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
      field_u = zeros(fieldSize, fieldSize) + h_u;
      field_v = zeros(fieldSize, fieldSize) + h_v;
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
    input = circshift(ih1 * gauss2d(1:fieldSize, 1:fieldSize, halfSize+1, halfSize+1, iw1, iw1), ...
      [round(ip1x) - halfSize; round(ip1y) - halfSize]) + ...
      circshift(ih2 * gauss2d(1:fieldSize, 1:fieldSize, halfSize+1, halfSize+1, iw2, iw2), ...
      [round(ip2x) - halfSize; round(ip2y) - halfSize]) + ...
      circshift(ih3 * gauss2d(1:fieldSize, 1:fieldSize, halfSize+1, halfSize+1, iw3x, inf), ...
      [round(ip3x) - halfSize; 0]) + ...
      circshift(ih4 * gauss2d(1:fieldSize, 1:fieldSize, halfSize+1, halfSize+1, inf, iw4y), ...
      [0; round(ip4y) - halfSize]);
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


