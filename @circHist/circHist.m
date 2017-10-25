classdef CircHist < handle
    %CircHist   Class representing a figure with a circular histogram. Constructing an
    %object creates a polar-coordinates axes containing a histogram. Circular statistics
    %(average angle, resultant vector length, Rayleigh test of uniformity and
    %circular-linear correlation) are automatically calculated using the CircStat toolbox
    %and saved as object properties. Note that this is a handle class, but that properties
    %of the plot can be accessed via properties and methods of the created object.
    %
    %   Requirements: CircStat toolbox (mathworks.com/matlabcentral/fileexchange/10676)
    %                 Ephys class
    %
    %   Usage:  CircHist(data,edges);
    %           CircHist(data,edges,Name,Value);
    %           obj = CircHist(___);
    %
    %   Notes and instructions:
    %   * To change the axis limits, use obj.setRLim([lower,upper]).
    %   * To change the scale-label, use obj.scaleBar.Label.String = 'my label'.
    %   * To add a tilted label to the degree-axis outside of the plot area, use
    %     obj.thetaLabel('my label',location).
    %   * To change visual properties, either use the name-value pairs for the constructor
    %     as specified below, or access the graphics-objects, namely obj.polarAxs for the
    %     coordinate system (font size, line width, etc.) and obj.scaleBar for the scale
    %     bar (line thickness, label, etc.). The scale bar is drawn anew each time the
    %     figure-window's size is changed; note that this may not work flawlessly and
    %     check that the scale bar matches the coordinate-grid after changing the figure
    %     size.
    %   * To adjust the bars, standard-deviation whiskers, phimax line and r line after
    %     plotting, use their handles which are saved as properties. Access them using dot
    %     notation, e.g., h = obj.phimaxH.lineWidth, and delete them using delete(h).
    %   * Note that the constructed CircHist-object handle is stored in AX.UserData, where
    %     AX is the axes the histogram has been plotted into. This may be useful if you
    %     plot a series of histograms and forget to store the CircHist-objects.
    %   * Each standard-deviation whisker consists of a long line representing the
    %     magnitude of the deviation and a very short, thick line that marks the tip of
    %     the deviation-line. Both line types are comprised in the handle-array obj.stdH;
    %     however, they can be separately accessed by using
    %     findobj(obj.polarAxs,'Tag',TYPETAG), where TYPETAG is either 'stdWhisk' for the
    %     "main" line or 'stdWhiskEnd' for the tips.
    %   * Consider creating a new figure for each histogram because there may be
    %     side-effects regarding the axis and the scale bar if the same axes-object or
    %     figure-window are used.
    %
    %
    %   Methods:
    %       setRLim([lower,upper])      Change axis limits (usage:
    %                                   obj.setRLim([lower,upper])) (get current limits by
    %                                   calling rlim)
    %       setThetaLabel(txt,location) Adds (or updates) a label saying TXT outside of
    %                                   the plot at the location specified by LOCATION,
    %                                   which may be one of the following characters:
    %                                   'topleft', 'topright', 'bottomleft'(default if
    %                                   omitted), 'bottomright'. (usage:
    %                                   obj.setThetaLabel('Direction','bottomright') ).
    %                                   Specify TXT as a cell array of characters to add
    %                                   line breaks. Access the created text-object via
    %                                   obj.thetaLabel.
    %       colorBar                    Change bar color (usage: obj.colorBar = newColor)
    %       barWidth                    Change bar width (usage: obj.barWidth = newWidth)
    %       colorStd                    Change standard-deviation-whisker color (usage as
    %                                   above)
    %       stdWidth                    Change standard-deviation-whisker width (usage as
    %                                   above)
    %       toPdf(fileName)             Save as pdf-file with the specified file name
    %       toPng(fileName,resol)       Save as png-file with the specified file name at
    %                                   the optionally specified resolution (e.g. '-r300')
    %       drawCirc(rho,lineSpec)      Draws a circle in the plot with radius RHO, line
    %                                   appearance optionally specified by LINESPEC (see
    %                                   doc LineSpec). Optionally returns the handle to
    %                                   the graphics object. Example usage:
    %                                   h = obj.drawCirc(15,'-g','LineWidth',4)
    %
    %   ---Required input:
    %   data            Either one vector of angle samples (distribution), a cell array of
    %                   such samples, or a two-column matrix of N length where N is the
    %                   number of bins, the first column contains the bar heights and the
    %                   second column the standard deviations. If the second column is
    %                   omitted, the standard deviation is considered zero. According to
    %                   the input data, specify the dataType property.
    %   edges           Edges (°) of histogram, e.g. [0:20:360] for 20° bins.
    %
    %   ---Optional Name-Value pair input:
    %   dataType        'distribution'(default)/'psth'. Type of input data: 'distribution'
    %                   treats the input data as distributions of degree values; the data
    %                   of each input vector are binned inside the specified edges, each
    %                   bin is averaged and the standard deviation is calculated. 'psth'
    %                   treats the data as already binned data in Hz, e.g. from a
    %                   peri-stimulus-time-histogram; if the input matrix has a second
    %                   column, it is taken as the standard-deviation values for each bar.
    %
    %   psthType        'count'(default)/'frequency'. Based on this, the plotted values
    %                   are either counts per bin or counts per second (Hz). Note that
    %                   BINSIZESEC may need to be specified, based on DATATYPE and
    %                   PSTHTYPE.
    %
    %   binSizeSec      Width of bins in seconds. Note that this needs to be specified if
    %                   DATATYPE is 'distribution' and PSTHTYPE is 'frequency', as well as
    %                   when DATATYPE is 'psth' and PSTHTYPE is 'count'.
    %
    %   areAxialData    True(default)/false, specifies whether or not input data are
    %                   axial; else they are considered circular. This is taken into
    %                   account for statistical computations (axial data are multiplied by
    %                   2 before calculation).
    %
    %   includePhimax   'on'(default)/'off', plots average angle.
    %
    %   phimax          Numeric value of the average angle. Should be specified if
    %                   DATATYPE is 'psth' and INCLUDEPHIMAX is 'on'.
    %
    %   includeR        'on'(default)/'off', plots resultant vector length (r) as a black
    %                   overlay bar on top of the average angle. The black bar's length
    %                   equals the r-value in percent of the coordinate-system diameter.
    %
    %   r               Numeric value of the resultant vector length. Should be specified
    %                   if DATATYPE is 'psth' and INCLUDER is 'on'.
    %
    %   baseLineShift   Numeric value (default = 2) specifying the factor that scales the
    %                   size of the offset in the center of the plot. The default value
    %                   produces nice results; lower the value (negative values allowed)
    %                   to increase the size, increase the value to decrease the size.
    %                   Specify it as NaN to produce an offset of 0. (This functionality
    %                   might need improvement, maybe just implement it as a percent value
    %                   of the plot diameter)
    %
    %   adjustSlope     Slope that defines how strong optical properties such as bar width
    %                   scale with bin-size; default = 0.3.
    %
    %   ax              Axes handle to plot diagram into; becomes POLARAXS property. Note
    %                   that the referenced axes must be a 'polaraxes'. (experimental
    %                   feature, working in principle, but the scale sometimes misbehaves).
    %
    %   colorBar        Color of bars (default = [0 .45 .74]; (Matlab blue)).
    %   colorStd        Color of standard-deviation lines (default = 'k').
    %   colorPhimax     Color of phimax line (default = [.85 .33 .1]; (orange)).
    %   colorR          Color of r line (default = 'k').
    %   fontSize        Font size of axis labels (default = 13).
    %
    %
    %
    %  ---Author: Frederick Zittrell
    %
    % See also polaraxes polarplot
    properties (SetAccess = immutable)
        data            % Required input: Data.
        edges           % Required input: Histogram edges;
                
        dataType        % Optional input; 'distribution'(default)/'psth'
        psthType        % Optional input; 'frequency'(default)/'count'
        binSizeSec      % Optional input; Width of bins (s)
        includePhimax   % Optional input; 'on'(default)/'off'
        phimax          % Optional input; Numeric value of the average angle
        includeR        % Optional input; 'on'(default)/'off'
        r               % Optional input; Numeric value of the resultant vector length
        baseLineShift   % Optional input; scaling factor for the size of the offset in the middle of the diagram
        adjustSlope     % Optional input; Slope for scaling of visual properties with bin size
        areAxialData    % Optional input; True(default)/false
        
        polarAxs        % Polaraxes handle. Change visual properties such as line width of the axes here.
        figH            % Handle to figure where diagram is plotted
        
        psthData        % PSTH data as plotted; 1st column average counts, 2nd column standard deviations
        rayleighP       % P-value of Rayleigh test of uniformity
        rayleighZ       % Z-value of Rayleigh test of uniformity
        corrAnP         % P-value of correlation analysis
        corrAnR         % R-value of correlation analysis (square this to get the coefficient of determination)
    end
    
    properties
        scaleBar        % Handle of scale bar. Use to access visual properties.
        axisLabel       % Label of scale bar as originally set (change the scaleBar.Label.String property to adjust)
        thetaLabel      % Label of the degree-axis (text-object, constructed via TEXT)
        
        phimaxH         % Handle to the phimax line
        rH              % Handle to the r line
        barH            % Array of handles to the bars (which are line objects)
        stdH            % Array of handles to the standard-deviation whiskers (which are line objects)
        
        colorBar        % Optional input; Color of bars (default = [0 .45 .74]; (Matlab blue))
        colorStd        % Optional input; Color of standard-deviation lines (default = 'k')
        colorPhimax     % Optional input; Color of phimax line (default = [.85 .33 .1]; (orange))
        colorR          % Optional input; Color of r line (default = 'k')
        fontSize        % Optional input; Font size of axis labels (default = 13)
        
        barWidth        % Width of bars. Change this property to adjust the bar width after plotting.
        stdWidth        % Width of standard-deviation whiskers. For adjustment after plotting.
        
        whiteCircH      % Handles to white bars that obscure the center of the axis
    end
    
    methods
        %% constructor
        function self = CircHist(data,edges,varargin)   
            %% CircHist obj = CircHist(data,edges,Name,Value)
            
            %% validate and parse input
            validateattributes(data,{'numeric','cell'},{'nonempty'});
            validateattributes(edges,{'numeric'},{'vector'});
            self.data  = data;
            self.edges = edges;
            
            if exist('isColorSpec','file')
                validColor = @isColorSpec; % custom function, may not be present
            else
                validColor = @(x)validateattributes(x,{'numeric','char'},{'vector'});
            end
            
            validScalarNum = @(N) isnumeric(N) & isscalar(N);
            
            %default values
            def.dataType      = 'distribution';
            def.psthType      = 'count';
            def.binSizeSec    = [];
            def.includePhimax = 'on';
            def.phimax        = [];
            def.includeR      = 'on';
            def.r             = [];
            def.baseLineShift = 2;
            def.adjustSlope   = 0.3;
            def.axialData     = true;
            def.axes          = [];
            def.colorBar      = [0 .45 .74]; %matlab blue
            def.colorStd      = 'k';
            def.colorPhimax   = [.85 .33 .1];
            def.colorR        = 'k';
            def.fontSize      = 13;
            
            pr = inputParser;
            addOptional(pr,'dataType',def.dataType,...
                @(str) any(strcmpi(str,{'distribution','psth'})));
            addOptional(pr,'psthType',  def.psthType,...
                @(str) any(strcmpi(str,{'frequency','count'})));
            addOptional(pr,'binSizeSec',def.binSizeSec,...
                @(x) validateattributes(x,{'numeric'},{'scalar'}));
            
            addParameter(pr,'includePhimax',def.includePhimax);
            addParameter(pr,'phimax'       ,def.phimax,validScalarNum);
            addParameter(pr,'includeR'     ,def.includeR);
            addParameter(pr,'r'            ,def.r,validScalarNum);
            addParameter(pr,'baseLineShift',def.baseLineShift,validScalarNum);
            addParameter(pr,'adjustSlope'  ,def.adjustSlope,validScalarNum);
            addParameter(pr,'areAxialData' ,def.axialData,...
                @(x) validateattributes(x,{'logical'},{'scalar'}));
            addParameter(pr,'ax'           ,def.axes...
                ,@(x) validateattributes(x,{'matlab.graphics.axis.PolarAxes'},{'scalar'}));
            
            addParameter(pr,'colorBar',   def.colorBar,   validColor);
            addParameter(pr,'colorStd',   def.colorStd,   validColor);
            addParameter(pr,'colorPhimax',def.colorPhimax,validColor);
            addParameter(pr,'colorR',     def.colorR,     validColor);
            addParameter(pr,'fontSize',   def.fontSize,   validScalarNum);
            
            parse(pr,varargin{:});
            
            self.dataType   = pr.Results.dataType;
            areDistribData  = strcmpi(self.dataType,'distribution');
            binSizeSec      = pr.Results.binSizeSec;
            self.binSizeSec = binSizeSec;
            self.psthType   = pr.Results.psthType;
            isFrequency     = strcmpi(self.psthType,'frequency');
            
            self.includePhimax = strcmpi(pr.Results.includePhimax,'on');
            includePhimax      = self.includePhimax;
            self.phimax        = pr.Results.phimax;
            phimax             = self.phimax;
            self.includeR      = strcmpi(pr.Results.includeR,'on');
            includeR           = self.includeR;
            self.r             = pr.Results.r;
            r                  = self.r;
            self.baseLineShift = pr.Results.baseLineShift;
            baseLineShift      = self.baseLineShift;
            self.adjustSlope   = pr.Results.adjustSlope;
            adjustSlope        = self.adjustSlope;
            self.areAxialData  = pr.Results.areAxialData;
            areAxialData       = self.areAxialData;
            ax                 = pr.Results.ax;
            
            self.colorBar    = pr.Results.colorBar;
            self.colorStd    = pr.Results.colorStd;
            self.colorPhimax = pr.Results.colorPhimax;
            colorPhimax      = self.colorPhimax;
            self.colorR      = pr.Results.colorR;
            colorR           = self.colorR;
            self.fontSize    = pr.Results.fontSize;
            fontSize         = self.fontSize;
            
            %% validate input permutations
            if isempty(binSizeSec) && isFrequency && areDistribData
                error(['For frequency data, specify bin width in seconds with ',...
                    'the binSizeSec parameter.']);
            end
            if isempty(binSizeSec) && ~isFrequency && ~areDistribData
                error(['PSTH-type frequency-data input needs binSizeSec to be ',...
                    'specified in order to display data as counts per bin.']);
            end
            
            %%
            if isFrequency %set axis label
                self.axisLabel = 'Frequency /Hz';
            else
                self.axisLabel = 'Counts per bin';
            end
            % initialize theta-label
            self.thetaLabel = text; % empty
            
            %%
            % deduce histogram data from edges
            binSizeDeg = abs(edges(1) - edges(2));
            binCentersDeg = edges(1:end-1) + binSizeDeg/2;
            
            %% operations on input data based on dataType
            if areDistribData
                if isnumeric(data) %if it is only a vector, pack it into a cell
                    validateattributes(data,{'numeric'},{'vector'});
                    data = {data(:)};
                end
                
                nSamples = numel(data);
                psth = nan(numel(edges)-1,nSamples);
                for s = 1:nSamples %calculate PSTHs
                    psthStruct = Ephys.psth(data{s},binSizeDeg);
                    psth(:,s) = psthStruct.count;
                end
                
                if isFrequency %convert to frequency
                    psth = psth / binSizeSec; end
                
                %calculate means and standard deviations
                psthData(:,1) = mean(psth,2);
                psthData(:,2) = std(psth,0,2);
                %phimax, r, rayleigh
                degPool = vertcat(data{:}); %column-vector of all data points
                self.phimax  = Ephys.phimax(degPool,areAxialData);
                phimax       = self.phimax;
                self.r       = Ephys.r(degPool,areAxialData);
                r            = self.r;
                [self.rayleighP,self.rayleighZ] = Ephys.rayl(degPool,areAxialData);
                pRayl = self.rayleighP;
                zRayl = self.rayleighZ;
            else %PSTH-data
                nSamples = NaN; %not feasible
                if isvector(data) %if it is a vector, use zeros for standard deviation
                    psthData = data(:); %columnize
                    psthData(:,2) = 0;
                else
                    psthData = data;
                end
                pRayl = NaN; %not feasible only from PSTH-data
                zRayl = NaN;
                
                if ~isFrequency %convert back to counts per bin
                    psthData = psthData * binSizeSec; end
            end
            self.psthData = psthData;
            % correlation analysis
            [self.corrAnR,self.corrAnP] = ...
                Ephys.corrAn(binCentersDeg(:),psthData(:,1),areAxialData);
            corrAnR = self.corrAnR;
            corrAnP = self.corrAnP;
            
            %% initialize figure, set visual properties
            if ~isempty(ax) && isvalid(ax)
                currFig = ax.Parent;
            else
                currFig = gcf;
            end
            currFig.Visible = 'off';
            self.figH = currFig;
            set(0,'currentfigure',currFig);
            polarplot(0);hold on
            set(currFig,'color',[1,1,1]) % white background
            self.polarAxs = gca;
            
            % self-reference in property for hyper-redundancy (this is actually quite
            % handy if you want to retrieve the circHist-object from a figure)
            self.polarAxs.UserData.circHistObj = self;
            
            polarAxs = self.polarAxs;
            polarAxs.ThetaZeroLocation = 'top';
            polarAxs.Tag = 'Polar';
            
            lineSp = '-'; %continuous lines
            self.barWidth = adjustSlope * binSizeDeg + 4; %bar width
            self.stdWidth = self.barWidth / 3;
            lineWPhimax = adjustSlope * binSizeDeg/10 + 2;
            lineWR = lineWPhimax * 1.5;
            
            %% draw bars and whiskers
            self.drawBars
            %%
            circR = max(rlim); %radius of plot in data units
            
            %calculate baseline shift from center; depends on bin size for aesthetic reasons
            if ~isnan(baseLineShift)
                %shift baseline by the fraction of the full radius and this value
                baseLineShiftDiv = adjustSlope * binSizeDeg + baseLineShift;
                baseLineOffset = -round(circR/baseLineShiftDiv); %baseline shift in data units
            else
                baseLineOffset = 0;
            end
            
            self.drawBaseLine
            
            %shift radius baseline
            polarAxs.RLim = [baseLineOffset,circR];
            
            %% phimax and r
            % based on AREAXIALDATA, the phimax- and r-lines are axes in the histogram, or
            % lines with the phimax as direction
            phimaxRad = deg2rad(phimax);
            if ~isempty(phimax) && includePhimax %plot phimax
                if areAxialData
                    thetaPhimax = [phimaxRad,phimaxRad+pi];
                    rhoPhimax = [circR,circR];
                else
                    thetaPhimax = [phimaxRad,phimaxRad];
                    rhoPhimax = rlim;
                end
                self.phimaxH = polarplot(self.polarAxs,thetaPhimax,rhoPhimax,lineSp...
                    ,'lineWidth',lineWPhimax,'color',colorPhimax,'Tag','phimax');
            end
            if ~isempty(r) && ~isempty(phimax) && includeR
                % make vector length relative to plot radius (after shift)
                rNorm = r * range(rlim) + baseLineOffset;
                if areAxialData
                    thetaR = [phimaxRad,phimaxRad+pi];
                    rohR = [rNorm,rNorm];
                else
                    thetaR = [phimaxRad,phimaxRad];
                    rohR = [min(rlim),rNorm];
                end
                self.rH = polarplot(self.polarAxs,thetaR,rohR,lineSp...
                    ,'lineWidth',lineWR,'color',colorR,'Tag','r');
            end
            
            %% edit axes
            polarAxs.Color = 'none'; %no background
            labels = polarAxs.ThetaTickLabel;
            % add °-sign to labels
            if ~strcmp(labels{1}(end),'°') % in the case of subsequent plotting into the same axis
                for n = 1:numel(labels)
                    labels{n} = [labels{n},'°'];end
                polarAxs.ThetaTickLabel=labels;
            end
            polarAxs.ThetaAxis.FontSize = fontSize;
            
            polarAxs.LineWidth = 1;
            polarAxs.GridColor = 'k';
            polarAxs.GridAlpha = 0.5;
            %nifty way to dynamically create minor grid lines in 10° spacing, skipping the
            %major grid lines. No one will ever want to understand this.
            minorGrid = 10:10:350;
            minorGrid = minorGrid(logical(mod(minorGrid,polarAxs.ThetaTick(2))));
            polarAxs.ThetaAxis.MinorTickValues = minorGrid;
            polarAxs.ThetaMinorGrid='on';
            polarAxs.MinorGridColor = 'k';
            polarAxs.MinorGridAlpha = 0.5;
            
            % draw white circle as background in the center
            self.drawWhiteCirc;
            
            %% title
            details = sprintf(['N = %u , phimax = %.2f°, r = %.4f\np_{Rayl} = %.3f, '...
                ,'Z_{Rayl} = %.4f, p_{C. an.} = %.3f, R²_{C. an.} = %.3f']...
                ,nSamples,phimax,r,pRayl,zRayl,corrAnP,corrAnR*corrAnR);
            title(details,'FontSize',9)
            
            %%
            colormap white %for axis appearance
            polarAxs.RTickLabel = [];
            
            if isempty(ax) % adjust figure-window size
                currFig.OuterPosition(1) = currFig.OuterPosition(1)/2; 
                currFig.OuterPosition(2) = currFig.OuterPosition(2)/2;
                currFig.OuterPosition(3) = currFig.OuterPosition(3)*1.2;
                currFig.OuterPosition(4) = currFig.OuterPosition(4)*1.5;
            end
            
            self.drawScale;
            currFig.SizeChangedFcn = {@self.redrawScale,self};

            currFig.Visible = 'on';
        end
        %%
        %%
        %% drawBars
        function drawBars(self)
            % Draws the bars. Needs to be done each time the r-axis limits change because
            % it's stupid.
            
            delete(self.barH);
            
            binSizeDeg = abs(self.edges(1) - self.edges(2));
            binCentersDeg = self.edges(1:end-1) + binSizeDeg/2;
            
            lowerLim = self.polarAxs.RLim(1); % lower limit of radius-axis
            
            % properties
            lineSp = '-';
            lineWBar = self.barWidth;
            clrBar = self.colorBar;
            
            % standard dev. bars and whiskers
            % whisker-endings consist of a short line starting at (whisker-end - whiskLen)
            % and ending at the whisker-end
            lineWStd = self.stdWidth;
            lineWStdWhisk = lineWBar * 0.7; %width
            clrStd = self.colorStd;
            whiskLen = .15;
            
            % draw
            nBars = numel(binCentersDeg);
            for n = 1:nBars
                currAng = deg2rad(binCentersDeg(n)); %angle of bar in rad
                currVal = self.psthData(n,1);
                if currVal == 0 || currVal <= lowerLim % skip bar-drawing
                    continue, end
                currStd = self.psthData(n,2);
                
                % specify the bar: a line between the base line of the plot (at theta =
                % angle and rho = 0) and the point specified by the angle and the
                % radius-value (at theta = angle and rho = bar height)
                thetaVal = [currAng,currAng];
                rhoVal = [0,currVal];
                % adjust bar-offset depending on axis-offset
                if lowerLim > 0, rhoVal(1) = lowerLim; end
                polarplot(self.polarAxs,thetaVal,rhoVal...
                    ,lineSp,'lineWidth',lineWBar,'color',clrBar,'Tag','histBar');
                
                if currStd > 0
                    % specify the standard deviation whisker: a line between the outer end
                    % of the bar (at theta = angle and rho = bar height) and the same end
                    % extended by the standard deviation (at theta = angle and rho = bar
                    % height + standard deviation)
                    thetaStd = [currAng,currAng];
                    rhoStd = [currVal,currVal+currStd];
                    polarplot(self.polarAxs,thetaStd,rhoStd,lineSp,'lineWidth',lineWStd...
                        ,'color',clrStd,'Tag','stdWhisk');
                    
                    % whisker-endings: short lines thicker than the whisker, starting at
                    % (whisker-end - whiskLen) and ending at the whisker-end
                    thetaStdWhisk = thetaStd;
                    rhoStdWhisk = [currVal+currStd-whiskLen,currVal+currStd];
                    polarplot(self.polarAxs,thetaStdWhisk,rhoStdWhisk,lineSp,'lineWidth',...
                        lineWStdWhisk,'color',clrStd,'Tag','stdWhiskEnd');
                end
            end
            self.barH = findobj(self.polarAxs,'Tag','histBar');
            self.stdH = findobj(self.polarAxs,'Tag','stdWhisk','-or','Tag','stdWhiskEnd');
            if ~isempty(self.phimaxH) && isvalid(self.phimaxH)
                uistack(self.phimaxH,'top');
            end
            if ~isempty(self.rH) && isvalid(self.rH)
                uistack(self.rH,'top');
            end
        end
        %% drawBaseLine
        function drawBaseLine(self)
            % Draws a base line at 0. This function should be called before drawWhiteCirc
            % is called, else the base line will be obscured by the white area
            
            delete(findobj(self.polarAxs,'Tag','baseLine')); % delete old if present
            % skip drawing if lower axis-limit > 0
            if self.polarAxs.RLim(1) > 0, return, end
            
            % properties
            tag = 'baseLine';
            lineSp = '-';
            lineW = 1;
            clr = self.colorBar; %baseline color, same as bar-color
            
            binSizeDeg = abs(self.edges(1) - self.edges(2));
            thetaBase = 0:deg2rad(binSizeDeg):2*pi; %degree-values between bars
            rhoBase = zeros(numel(thetaBase),1);
            polarplot(self.polarAxs,thetaBase,rhoBase,lineSp,'lineWidth',lineW...
                ,'color',clr,'Tag',tag);
            uistack(findobj(self.polarAxs,'Tag',tag),'bottom');
            uistack(findobj(self.polarAxs,'Tag',tag),'up');
        end
        %% drawScale; draw scale bar
        function drawScale(self)
            %drawScale Draws the scale bar. Used as SizeChangedFcn for figure so it is
            %drawn each time the figure size is changed. The scale bar is actually a
            %colorbar-object, thus it does not behave as neat as a conventional axis.
            
            self.figH.Visible = 'off'; % as recommended for SizeChangedFcn operations
            pAx = self.polarAxs;
            scl = self.scaleBar;
            
            %default font size. This line should be completely useless, but for some
            %reason, there are strange effects without it. The fact that the bar-drawing
            %is executed via SIZECHANGEDFCN seems to interrupt the link between the font
            %properties of the bar and the font properties of the parent axes, leading to
            %strang effects.
            fontsz = 13; %#ok

            initialDraw = isempty(scl);
            if initialDraw % if the bar is drawn for the first time, use predefined label and font size
                label = self.axisLabel;
                fontsz = self.fontSize;
            else % read out (assumingly predefined) properties and use these for the next drawing
                label = scl.Label.String;
                fontsz = scl.Label.FontSize;
                oldLineWidth = scl.LineWidth;
            end
            
            delete(scl); %deletes the previous bar
            
            scl = colorbar('Location','manual');
            self.scaleBar = scl; % re-assign
            
            % with this link, the label font-name is changed with the corresponding
            % axes-property using set(gca,'FontName',fontName), which is the default
            % behavior of colorbar labels.
            pAx.UserData.fontNameLink =  linkprop([pAx,scl],'FontName');
            
            sclUnitsOld = scl.Units;
            scl.Units = 'pixels';
            polarAxsUnitsOld = pAx.Units;
            pAx.Units = 'pixels';
            
            polarPos = pAx.Position; %position property = [left,bottom,width,height]
            scl.Position(1) = polarPos(1);
            polarBot = polarPos(2);
            polarHeight = polarPos(4);
            
            sclWidth = 0.015;
            sclHeight = polarHeight/2; %halve the height -> bar spans from center to top
            lowerLim = pAx.RLim(1);
            if lowerLim < 0
                offset = sclHeight * abs(lowerLim)/range(pAx.RLim);
            else, offset = 0;
            end
            sclHeight = sclHeight - offset; %adjust baseline shift
            sclBot = polarBot + polarHeight - sclHeight; %calculate baseline position
            scl.Position = [pAx.Position(1),sclBot,sclWidth,sclHeight];
            %move scale left if it overlays the polar-axis label (hard-coded, empirical values, sadly)
            if offset < 15, moveLeft = 15;
            else,           moveLeft =  0;
            end
            scl.Position(1) = pAx.Position(1) - moveLeft; %move
            
            scl.Label.String = label;
            lowerLimScl = lowerLim;
            if lowerLimScl < 0, lowerLimScl = 0; end
            scl.Limits = [lowerLimScl,max(pAx.RLim)];
            sclTicks = pAx.RTick; %use ticks as produced by the polarplot function
            sclTicks = sclTicks(sclTicks >= lowerLimScl);
            scl.Ticks = sclTicks;
            scl.Box = 'off';
            scl.TickLength = 0.04;
            scl.FontSize = fontsz;
            
            scl.Units = sclUnitsOld;
            pAx.Units = polarAxsUnitsOld;
            
            if ~initialDraw
                scl.LineWidth = oldLineWidth;
            end
            self.figH.Visible = 'on';
        end
        %% draw white circle in the middle to obscure the axis-lines
        function drawWhiteCirc(self)
            if ~isempty(self.whiteCircH), delete(self.whiteCircH); end
            
            baseLineOffset = rlim;
            baseLineOffset = baseLineOffset(1);
            if baseLineOffset >= 0, return, end
            
            rhoWhite = [baseLineOffset,0];
            thetaWhite = deg2rad(0:2:358);
            for t = 1:numel(thetaWhite)
                currAng = thetaWhite(t);
                polarplot(self.polarAxs,[currAng,currAng],rhoWhite...
                    ,'linestyle','-','color','w','linewidth',3,'Tag','whiteCirc');
            end
            self.whiteCircH = findobj(self.polarAxs,'Tag','whiteCirc');
            uistack(self.whiteCircH,'bottom'); % move lines to lowest graphical layer
        end
        %% change scale limits
        function setRLim(self,limits)
            %setRLim Change scale limits specified by the two-element vector LIMITS ==
            %[lower,upper]. Get the current limits by calling rlim.
            %
            % circHistObj.setRLim(limits); where LIMITS == [lower,upper]
            %
            rlim(limits); %change limits
            
%             if limits(1) > 0
%                 rDataOffset = limits(1);
%             else
%                 rDataOffset = 0;
%             end
%             
%             nBars = numel(self.barH);
%             for n = 1:nBars
%                 currBar = self.barH(n);
%                 barHeight = currBar.RData(2);
%                 currBar.RData = [rDataOffset,barHeight];                
%             end

            % update line data
            self.drawBars
            if isvalid(self.phimaxH)
                if self.areAxialData
                    self.phimaxH.RData(:) = limits(2);
                else
                    self.phimaxH.RData = limits;
                end
            end
            if isvalid(self.rH)
                rNorm = self.r * range(limits) + limits(1);
                if self.areAxialData
                    self.rH.RData(:) = rNorm;
                else
                    self.rH.RData = [limits(1),rNorm];
                end                
            end
            
            % Sometimes, a warning with this identifier is issued for obscure reasons. To
            % suppress this warning, it is converted to an error by calling
            % warning('error',_). This error can then be caught and ignored.
            wrn = warning('error','MATLAB:callback:error'); %#ok<CTPCT>
            try   self.drawScale; % update scale
            catch ME
                if ~strcmp(ME.identifier,wrn.identifier)
                    rethrow(ME); end % re-throw error if it is not this error
            end
            warning(wrn); % set error back to being a warning
            
            self.drawBaseLine
            self.drawWhiteCirc; % adjust background
        end
        %% set theta-axis label
        function setThetaLabel(self,txt,location)
            %setThetaLabel Add (or update) a label to the theta-axis, label-text specified by TXT, location specified by LOCATION, which may be 'bottomleft'(default if omitted), 'bottomright', 'topleft' or 'topright'.
            %
            % circHistObj.setThetaLabel('My label','topright');
            %
            
            assert(ischar(txt) || isstring(txt) || iscellstr(txt) ...
                || all(cellfun(@isstring,txt))...
                ,'TXT input must be a CHAR, STRING or cell array of CHARS or STRINGS.');
            
            locations = {'bottomleft','bottomright','topleft','topright'};
            if nargin < 3, location = locations{1};
            else,          location = lower(location);
                assert(any(strcmp(location,locations)) ...
                    ,['LOCATION input must be one of the following: ' ...
                    ,strjoin(locations,',')]);
            end
            
            % based on THETADIR and THETAZEROLOCATION, the label-theta angle needs to be adjusted so the label is in die specified corner
            if strcmp(self.polarAxs.ThetaDir,'counterclockwise')
                  thetaDirSign = +1;
            else, thetaDirSign = -1;
            end
            
            switch self.polarAxs.ThetaZeroLocation
                case 'top',     thOffsetZero =   0;
                case 'right',   thOffsetZero =  90;
                case 'bottom',  thOffsetZero = 180;
                case 'left',    thOffsetZero = 270;
            end
            thOffsetZero = thOffsetZero * thetaDirSign;
            
            switch location
                case 'topleft',     thOffsetLoc =   0; txtRot =  45;
                case 'bottomleft',  thOffsetLoc =  90; txtRot = -45;
                case 'bottomright', thOffsetLoc = 180; txtRot =  45;
                case 'topright',    thOffsetLoc = 270; txtRot = -45;
            end
            thOffsetLoc = thOffsetLoc * thetaDirSign;
            
            th = 45 * thetaDirSign + thOffsetLoc + thOffsetZero;
            
            rlims = self.polarAxs.RLim;
            delete(self.thetaLabel);
            self.thetaLabel = text(deg2rad(th),rlims(2) + range(rlims)*0.2,txt ...
                ,'HorizontalAlignment','center','VerticalAlignment','cap' ...
                ,'Rotation',txtRot,'FontSize',self.fontSize);
            % make label-color be linked to theta-axis color (default behavior for labels)
            % for some reason, the link only works when THETALABEL.COLOR is changed, not
            % the other way around ...
            self.polarAxs.UserData.thetaLabelColorLink = ...
                linkprop([self.polarAxs.ThetaAxis,self.thetaLabel],'Color');
        end
        %% change bar color
        function set.colorBar(self,color)
            self.colorBar = color;
            lineObjArr = findobj(self.polarAxs,'Tag','histBar','-or','Tag','baseLine');
            set(lineObjArr,'color',color);
        end
        %% change bar width
        function set.barWidth(self,width)
            self.barWidth = width;
            set(self.barH,'lineWidth',width);
        end
        %% change whisker color
        function set.colorStd(self,color)
            % Note to self: The color- and width-changing functions could also be
            % implemented by linking the respective property of all objects using linkprop
            % and then only changing the property of one.
            %
            %   circHistObj.colorStd([.5,.5,.5]);
            self.colorStd = color;
            set(self.stdH,'color',color);
        end
        %% change whisker width
        function set.stdWidth(self,width)
            % Since each whisker consists of two line-objects with different widths (the
            % "main" line and the ending), the width of the ending is scaled
            % proportionally to the width-change of the "main" line.
            %
            %   circHistObj.stdWidth(2);
            
            oldWidth = self.stdWidth;
            self.stdWidth = width;
            
            %happens at object construction; this is NOT elegant and I am sorry
            if isempty(oldWidth), return; end
            
            scalingFactor = width / oldWidth;
            
            %get handles of different line objects
            maskLines = strcmp({self.stdH.Tag},'stdWhisk');
            stdMain = self.stdH(maskLines);
            stdEnding = self.stdH(~maskLines);
            
            oldEndingWidth = stdEnding(1).LineWidth; %scale proportionally
            newEndingWidth = oldEndingWidth * scalingFactor;
            
            set(stdMain,'lineWidth',width);
            set(stdEnding,'lineWidth',newEndingWidth);
        end
        %% save to pdf
        function toPdf(self,fileName)
            %toPdf  Save histogram as (FILENAME).pdf.
            %
            %   circHistObj.toPdf(fileName);
            
            if exist('toPdf','file') % call custom function if available
                toPdf(self.figH,fileName);
            else
                print(self.figH,fileName,'-dpdf','-fillpage','-painters');
            end
        end
        %% save to png
        function toPng(self,fileName,resol)
            %toPng Save histogram as (FILENAME).png at the optionally specified resolution
            %(default = 90 dpi). Specify RESOL as a string of the pattern '-r90'.
            %
            %   obj.toPng(filename);
            %   obj.toPng(filename,resol);
            
            if nargin < 3
                resol = '-r90';
            end
            
            if exist('toPng','file') % call custom function if available
                toPng(self.figH,fileName,resol);
            else
                print(self.figH,fileName,'-dpng','-opengl',resol);
            end
        end
        %% drawCirc
        function hOut = drawCirc(self,rho,varargin)
            %drawCirc Draws a circle in the plot, radius specified by RHO, appearance
            %optionally specified by additional parameters which must be Name-Value pairs
            %as accepted by POLARPLOT. Optionally returns the graphics-object handle.
            %
            %   obj.drawCirc(rho);
            %   obj.drawCirc(rho,Name,Value);
            %   h = obj.drawCirc(___);
            %
            % See also polarplot
            
            assert(isnumeric(rho) && isscalar(rho) ...
                ,'First input RHO must be a numerical scalar.');
                        
            theta = linspace(0,2*pi,100);
            h = polarplot(theta,repmat(rho,size(theta)),varargin{:});
            
            if isvalid(self.phimaxH)
                uistack(self.phimaxH,'top'); end
            if isvalid(self.rH)
                uistack(self.rH,'top'); end
            
            if nargout > 0 
                hOut = h; end
        end
    end
    methods (Static)
        %%
        function redrawScale(~,~,circHistObj)
           set(0,'currentfigure',gcbf);
           circHistObj.drawScale;
       end
    end
end