classdef CircHist < handle
    %CircHist   Class representing a figure with a circular histogram. Constructing an
    % object creates a polar-coordinates axes containing a histogram. Circular statistics
    % (average angle, 95 % confidence interval, resultant vector length, Rayleigh test of
    % uniformity and circular-linear correlation) are automatically calculated using the
    % CircStat toolbox and stored as object properties. Note that this is a handle class,
    % but that properties of the plot can be accessed via properties and methods of the
    % created object. Input data may be either angular distributions or already-binned
    % data.
    %
    %   Requirements: CircStat toolbox (mathworks.com/matlabcentral/fileexchange/10676)
    %
    %   Usage:  CircHist(data,edges);
    %           CircHist(data,nBins);
    %           CircHist(data,edges,Name,Value);
    %           CircHist(data,nBins,Name,Value);
    %           obj = CircHist(___);
    %
    %   Notes and instructions:
    %   * To change the radius-axis limits, use obj.setRLim([lower,upper]); to change the
    %     angle-axis limits, use obj.thetaLim = [lower,upper].
    %   * To change the scale-label, use obj.axisLabel = 'my label'.
    %   * To add a tilted label to the degree-axis outside of the plot area, use
    %     obj.setThetaLabel('my label',location).
    %   * To change visual properties, either use the name-value pairs for the constructor
    %     as specified below, or access the graphics-objects, namely obj.polarAxs for the
    %     coordinate system (font size, line width, etc.) and obj.scaleBar for the scale
    %     bar (line thickness etc.). The scale bar is drawn anew each time the
    %     figure-window's size is changed; note that this may not work flawlessly and
    %     assure that the scale bar matches the coordinate-grid after changing the figure
    %     size.
    %   * To adjust the bars, standard-deviation whiskers, average-angle line and r line
    %     after plotting, use their handles which are saved as properties. Access them
    %     using dot notation, e.g., h = obj.avgAngH.lineWidth, and delete them using
    %     delete(h) (not: obj.avgAngH.lineWidth = []).
    %   * Note that the constructed CircHist-object handle is stored in
    %     obj.polarAxs.UserData.circHistObj. This may be useful if you plot a series of
    %     histograms and forget to store the CircHist-objects.
    %   * Each standard-deviation whisker consists of a long line representing the
    %     magnitude of the deviation and a very short, thick line that marks the tip of
    %     the deviation-line. Both line types are comprised in the handle-array obj.stdH;
    %     however, they can be separately accessed by using
    %     findobj(obj.polarAxs,'Tag',TYPETAG), where TYPETAG is either 'stdWhisk' for the
    %     "main" line or 'stdWhiskEnd' for the tips.
    %   * Analogously, access the 95 % confidence-interval line(s) by accessing the
    %     handles in obj.avgAngCiH or by using findobj(obj.polarAxs,'Tag',TYPETAG), where
    %     TYPETAG is either 'avgAngCiWhisk' or 'avgAngCiWhiskEnd'.
    %   * Consider creating a new figure for each histogram because there may be
    %     side-effects regarding the axis and the scale bar if the same axes-object or
    %     figure-window are used.
    %   * If you want the angle-axis units to be in radians, first change the
    %     POLARAXES.THETAAXIS.TICKLABELSMODE property to 'auto', then change
    %     THETAAXISUNITS to 'radians'.
    %   * The axis-grid-lines in the center of the plot are obscured by white bars for
    %     cosmetic reasons. To access the bars, e.g. to change their color, refer to
    %     obj.whiteDiskH.
    %
    %
    %   Methods and noteworthy properties:
    %       setRLim([lower,upper])      Change axis limits (usage:
    %                                   obj.setRLim([lower,upper])) (get current limits by
    %                                   calling rlim)
    %       setThetaLabel(txt,location) Adds (or updates) a label saying TXT outside of
    %                                   the plot at the location specified by LOCATION,
    %                                   which may be one of the following: 'topleft',
    %                                   'topright', 'bottomleft'(default if omitted),
    %                                   'bottomright'. (usage:
    %                                   obj.setThetaLabel('Direction','bottomright') ).
    %                                   Specify TXT as a cell array of characters to add
    %                                   line breaks. Access the created text-object via
    %                                   obj.thetaLabel.
    %       scaleBarSide                Change the side ('left'/'right') of the
    %                                   theta-axis scale bar (usage: obj.scaleBarSide =
    %                                   'right').
    %       colorBar                    Change bar color (usage: obj.colorBar = newColor).
    %       barWidth                    Change bar width (usage: as above).
    %       colorStd                    Change standard-deviation-whisker color (usage: as
    %                                   above).
    %       stdWidth                    Change standard-deviation-whisker width (usage: as
    %                                   above).
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
    %   data            Either one vector of angle samples (distribution) in degree (!), a
    %                   cell array of such samples, or a two-column matrix of N length
    %                   where N is the number of bins, the first column contains the
    %                   counts per bin and the second column the standard deviations. If
    %                   the second column is omitted, the standard deviations are
    %                   considered zero. According to the input data, specify the DATATYPE
    %                   property accordingly.
    %   edges           Edges of histogram, either specified by a vector of degree-values,
    %                   e.g. [0:20:360] for 20� bins, or by a scalar integer specifying
    %                   the number of bins in a full circle.
    %
    %   ---Optional Name-Value pair input:
    %   dataType        'distribution'(default)/'histogram'. Type of input data:
    %                   'distribution' treats the input data as distributions of angles
    %                   (in degree); the data of each input vector are binned inside the
    %                   specified edges, each bin is averaged and the standard deviation
    %                   is calculated. 'histogram' treats the data as already-binned data
    %                   (counts per bin); if the input matrix has a second column, it is
    %                   taken as the standard-deviation values for each bar. If you want
    %                   to plot PSTH-data (peri-stimulus-time-histogram) where your
    %                   bin-heights should be converted to frequencies, specify HISTTYPE
    %                   and BINSIZESEC accordingly.
    %
    %   histType        'count'(default)/'frequency'. Based on this, the plotted values
    %                   are either counts per bin or converted to counts per second
    %                   (feasible for PSTH-data). Note that BINSIZESEC needs to be
    %                   specified if HISTTYPE is 'frequency'.
    %
    %   binSizeSec      Width of bins in seconds in case a conversion from angles to time
    %                   is feasible. Note that this needs to be specified if HISTTYPE is
    %                   'frequency'.
    %
    %   areAxialData    True/false(default), specifies whether or not input data are
    %                   axial; else they are considered circular. This is taken into
    %                   account for statistical computations (axial data are multiplied by
    %                   2 before calculation) and determines whether the average-angle
    %                   line is plotted as an axis (axial data) or as a direction from the
    %                   plot center (circular data).
    %
    %   scaleBarSide    'left'(default)/'right', specifies the side of the rho-axis scale
    %                   bar. Access the scalebar (actually a COLORBAR object) via
    %                   OBJ.SCALEBAR. Change the side after plotting by changing
    %                   OBJ.SCALEBARSIDE to 'left' or 'right'.
    %
    %   thetaLim        [lowerDeg,upperDeg] (default == [0,360]); specifies the angle-axis
    %                   limits of the polarplot. May be changed after plotting.
    %
    %   drawAvgAng      'on'(default)/'off', plots average angle as a line.
    %
    %   avgAng          Numeric value (�) of the average angle. If specified, AVGANG is
    %                   not calculated from the data and this value is taken and plotted
    %                   instead.
    %
    %   drawAvgAngCi    'on'(default)/'off', plots 95 % confidence interval as a circle
    %                   segment outside of the plot.
    %
    %   avgAngCi        Numeric value (�) of 95 % confidence interval (AVGANG +/-
    %                   AVGANGCI) to plot; same behavior as for AVGANG.
    %
    %   drawR           'on'(default)/'off', plots resultant vector length (r) as a black
    %                   overlay bar on top of the average angle. The black bar's length
    %                   equals the r-value in percent of the average-angle-line length.
    %
    %   r               Numeric value of the resultant vector length; same behavior as for
    %                   AVGANG.
    %
    %   baseLineOffset  Numeric value (default = 2) specifying the factor that scales the
    %                   size of the offset in the center of the plot. The default value
    %                   produces nice results; lower the value (negative values allowed)
    %                   to increase the size, increase the value to decrease the size.
    %                   Specify it as NaN to produce an offset of 0. (This functionality
    %                   might need improvement, maybe just implement it as a percent value
    %                   of the plot diameter; as of now, it depends on bin-size and
    %                   ADJUSTSLOPE)
    %
    %   adjustSlope     Slope that defines how strongly optical properties such as bar
    %                   width scale with bin-size; default = 0.3.
    %
    %   ax              Axes handle to plot diagram into; becomes POLARAXS property. Note
    %                   that the referenced axes must be a POLARAXES. (experimental
    %                   feature, working in principle, but the scale sometimes
    %                   misbehaves).
    %
    %   colorBar        Color of bars (default = [0 .45 .74]; (Matlab blue)).
    %   colorStd        Color of standard-deviation lines (default = 'k').
    %   colorAvgAng     Color of average-angle line (default = [.85 .33 .1]; (orange)).
    %   colorR          Color of r line (default = 'k').
    %   fontSize        Font size of axis labels (default = 13).
    %
    %
    %
    %  ---Author: Frederick Zittrell
    %
    % See also polaraxes polarplot CircStat
    properties (SetAccess = immutable)
        data            % Required input: Data.
        edges           % Required input: Histogram edges;
                
        dataType        % Optional input; 'distribution'(default)/'histogram'
        histType        % Optional input; 'frequency'(default)/'count'
        binSizeSec      % Optional input; Width of bins (s)
        drawAvgAng      % Optional input; 'on'(default)/'off', converted to boolean
        avgAng          % Optional input; Numeric value of the average angle
        drawAvgAngCi    % Optional input; 'on'(default)/'off', converted to boolean
        avgAngCi        % Optional input; Numeric value of the 95 % confidence interval bounds (average +/ bounds)
        drawR           % Optional input; 'on'(default)/'off', converted to boolean
        r               % Optional input; Numeric value of the resultant vector length
        baseLineOffset  % Optional input; scaling factor for the size of the offset in the middle of the diagram
        adjustSlope     % Optional input; Slope for scaling of visual properties with bin size
        areAxialData    % Optional input; True/false(default)
        
        polarAxs        % Polaraxes handle. Change visual properties such as line width of the axes here.
        figH            % Handle to figure where diagram is plotted
        
        histData        % Histogram data as plotted; 1st column average counts, 2nd column standard deviations
        rayleighP       % P-value of Rayleigh test of uniformity
        rayleighZ       % Z-value of Rayleigh test of uniformity
        corrAnP         % P-value of correlation analysis
        corrAnR         % R-value of correlation analysis (square this to get the coefficient of determination)
    end
    
    properties
        scaleBar        % Handle of scale bar. Use to access visual properties.
        axisLabel       % Label of scale bar as originally set (change the scaleBar.Label.String property to adjust)
        scaleBarSide    % Side of the scale bar, either 'left'(default) or 'right'
        thetaLabel      % Label of the degree-axis (text-object, constructed via TEXT)
        thetaLim        % Angle-axis limits [lower,upper]
        
        avgAngH         % Handle to the average-angle line
        avgAngCiH       % Handle to confidence interval lines
        rH              % Handle to the r line
        barH            % Array of handles to the bars (which are line objects)
        stdH            % Array of handles to the standard-deviation whiskers (which are line objects)
        
        colorBar        % Optional input; Color of bars (default = [0 .45 .74]; (Matlab blue))
        colorStd        % Optional input; Color of standard-deviation lines (default = 'k')
        colorAvgAng     % Optional input; Color of average-angle line (default = [.85 .33 .1]; (orange))
        colorR          % Optional input; Color of r line (default = 'k')
        colorAvgAngCi   % Optional input; Color of confidence interval line (default = 'k')
        fontSize        % Optional input; Font size of axis labels (default = 13)
        
        barWidth        % Width of bars. Change this property to adjust the bar width after plotting.
        stdWidth        % Width of standard-deviation whiskers. For adjustment after plotting.
        
        whiteDiskH      % Handles to white bars that obscure the center of the axis
    end
    
    methods
        %% constructor
        function self = CircHist(data,edges,varargin)   
            %CircHist constructor
            %       obj = CircHist(data,edges,Name,Value)
            
            %% validate and parse input
            validateattributes(data,{'numeric','cell'},{'nonempty'});
            validateattributes(edges,{'numeric'},{'vector'});
            self.data  = data;
            if isscalar(edges)
                edges = 0 : (360/edges) : 360; end
            self.edges = edges;
            
            if exist('isColorSpec','file')
                  validColor = @isColorSpec; % custom function, may not be present
            else, validColor = @(x)validateattributes(x,{'numeric','char'},{'vector'});
            end
            
            % validates string as 'on' or 'off' and returns TRUE for 'on', FALSE for 'off'
            isOnOff = @(str) strcmp('on',validatestring(str,{'on','off'}));
            
            validScalarNum = @(N) isnumeric(N) & isscalar(N);
            
            %default values
            def.dataType       = 'distribution';
            def.histType       = 'count';
            def.binSizeSec     = [];
            def.thetaLim       = [0,360];
            def.drawAvgAng     = 'on';
            def.avgAng         = [];
            def.drawAvgAngCi   = 'on';
            def.avgAngCi       = [];
            def.drawR          = 'on';
            def.r              = [];
            def.baseLineOffset = 2;
            def.adjustSlope    = 0.3;
            def.axialData      = false;
            def.axes           = [];
            def.colorBar       = [0 .45 .74]; %matlab blue
            def.colorStd       = 'k';
            def.colorAvgAng    = [.85 .33 .1];
            def.colorR         = 'k';
            def.colorAvgAngCi  = 'k';
            def.fontSize       = 13;
            def.scaleBarSide   = 'left';
            
            pr = inputParser;
            addOptional(pr,'dataType',def.dataType,...
                @(str) any(strcmpi(str,{'distribution','histogram'})));
            addOptional(pr,'histType',  def.histType,...
                @(str) any(strcmpi(str,{'frequency','count'})));
            addOptional(pr,'binSizeSec',def.binSizeSec,...
                @(x) validateattributes(x,{'numeric'},{'scalar'}));
            
            addParameter(pr,'drawAvgAng'    ,def.drawAvgAng);
            addParameter(pr,'avgAng'        ,def.avgAng,validScalarNum);
            addParameter(pr,'drawAvgAngCi'  ,def.drawAvgAngCi);
            addParameter(pr,'avgAngCi'      ,def.avgAngCi);
            addParameter(pr,'drawR'         ,def.drawR);
            addParameter(pr,'r'             ,def.r,validScalarNum);
            addParameter(pr,'baseLineOffset',def.baseLineOffset,validScalarNum);
            addParameter(pr,'adjustSlope'   ,def.adjustSlope,validScalarNum);
            addParameter(pr,'areAxialData'  ,def.axialData,...
                @(x) validateattributes(x,{'logical'},{'scalar'}));
            addParameter(pr,'ax'            ,def.axes...
                ,@(x) validateattributes(x,{'matlab.graphics.axis.PolarAxes'},{'scalar'}));
            
            addParameter(pr,'scaleBarSide', ...
                def.scaleBarSide,@(str) any(strcmpi(str,{'left','right'})));
            addParameter(pr,'thetaLim',def.thetaLim, ...
                @(x) validateattributes(x,{'numeric'},{'numel',2}));
            addParameter(pr,'colorBar',   def.colorBar,   validColor);
            addParameter(pr,'colorStd',   def.colorStd,   validColor);
            addParameter(pr,'colorAvgAng',def.colorAvgAng,validColor);
            addParameter(pr,'colorR',     def.colorR,     validColor);
            addParameter(pr,'colorAvgAngCi',def.colorAvgAngCi,validColor);
            addParameter(pr,'fontSize',   def.fontSize,   validScalarNum);
            
            parse(pr,varargin{:});
            
            self.dataType   = pr.Results.dataType;
            areDistribData  = strcmpi(self.dataType,'distribution');
            binSizeSec      = pr.Results.binSizeSec;
            self.binSizeSec = binSizeSec;
            self.histType   = pr.Results.histType;
            isFrequency     = strcmpi(self.histType,'frequency');
            
            thetaLim             = pr.Results.thetaLim;
            self.drawAvgAng      = isOnOff(pr.Results.drawAvgAng);
            drawAvgAng           = self.drawAvgAng;
            self.avgAng          = pr.Results.avgAng;
            avgAng               = self.avgAng;
            self.drawAvgAngCi    = isOnOff(pr.Results.drawAvgAngCi);
            drawAvgAngCi         = self.drawAvgAngCi;
            self.avgAngCi        = pr.Results.avgAngCi;
            avgAngCi             = self.avgAngCi;
            
            self.drawR          = isOnOff(pr.Results.drawR);
            drawR               = self.drawR;
            self.r              = pr.Results.r;
            r                   = self.r;
            self.baseLineOffset = pr.Results.baseLineOffset;
            baseLineOffset      = self.baseLineOffset;
            self.adjustSlope    = pr.Results.adjustSlope;
            adjustSlope         = self.adjustSlope;
            self.areAxialData   = pr.Results.areAxialData;
            areAxialData        = self.areAxialData;
            ax                  = pr.Results.ax;
            
            scaleBarSide        = pr.Results.scaleBarSide;
            self.colorBar       = pr.Results.colorBar;
            self.colorStd       = pr.Results.colorStd;
            self.colorAvgAng    = pr.Results.colorAvgAng;
            colorAvgAng         = self.colorAvgAng;
            self.colorR         = pr.Results.colorR;
            colorR              = self.colorR;
            self.colorAvgAngCi  = pr.Results.colorAvgAngCi;
            colorAvgAngCi       = self.colorAvgAngCi;
            self.fontSize       = pr.Results.fontSize;
            fontSize            = self.fontSize;
            
            %% validate HISTTYPE-BINSIZESEC combination
            if isempty(binSizeSec) && isFrequency
                error('To obtain frequency-data (counts per second), specify BINSIZESEC.');
            end
            %%
            % deduce bin data from edges
            binSizeDeg = abs(edges(1) - edges(2));
            binSizeRad = deg2rad(binSizeDeg);
            binCentersDeg = edges(1:end-1) + binSizeDeg/2;
            binCentersRad = deg2rad(binCentersDeg');
            
            %% use this function handle to transform axial data on demand
            % axial dimension, input to CIRCSTAT functions
            if areAxialData, axialDim = 2; else, axialDim = 1; end
            axTrans = @(x)circ_axial(x,axialDim);
            
            %% operations on input data based on dataType
            if areDistribData
                if isnumeric(data) %if it is only a vector, pack it into a cell
                    validateattributes(data,{'numeric'},{'vector'});
                    data = {data(:)};
                end
                
                nSamples = numel(data);
                binnedData = nan(numel(edges)-1,nSamples);
                for s = 1:nSamples % calculate histogram counts
                    if exist('histcounts','file') || exist('histcounts','builtin')
                        [counts,~] = histcounts(data{s},edges);
                    else % for Matlab versions <2014
                        [counts,~] = histc(data{s},edges);
                        counts(end) = [];
                    end
                    binnedData(:,s) = counts;
                end
                
                if isFrequency %convert to counts per second
                    binnedData = binnedData / binSizeSec; end
                
                % calculate means and standard deviations
                histData(:,1) = mean(binnedData,2);
                histData(:,2) = std(binnedData,0,2);
                % avgAng, r, rayleigh
                degPool = vertcat(data{:}); % column-vector of all data points
                radPool = deg2rad(degPool);
                if isempty(avgAng)
                    self.avgAng = mod(rad2deg(circ_mean(axTrans(radPool))),360) /axialDim;
                    avgAng      = self.avgAng;
                end
                if isempty(r)
                    self.r = circ_r(axTrans(radPool));
                    r      = self.r;
                end
                if isempty(avgAngCi)
                    self.avgAngCi = ...
                        rad2deg(circ_confmean(axTrans(radPool),0.05)) / axialDim;
                    avgAngCi      = self.avgAngCi;
                end
                
                [self.rayleighP,self.rayleighZ] = circ_rtest(axTrans(radPool));
                pRayl = self.rayleighP;
                zRayl = self.rayleighZ;
            else % histogram data
                nSamples = NaN; %not feasible
                if isvector(data) %if it is a vector, use zeros for standard deviation
                    histData = data(:); %columnize
                    histData(:,2) = 0;
                else
                    histData = data;
                end
                histCounts = histData(:,1);
                
                if isempty(avgAng)
                    self.avgAng  = ...
                        rad2deg(circ_mean(axTrans(binCentersRad),histCounts)) / axialDim;
                    avgAng       = self.avgAng;
                end
                if isempty(r)
                    self.r = circ_r(axTrans(binCentersRad),histCounts,binSizeRad);
                    r      = self.r;
                end
                if isempty(avgAngCi)
                    self.avgAngCi = rad2deg(circ_confmean( ...
                        axTrans(binCentersRad),0.05,histCounts,binSizeRad)) / axialDim;
                    avgAngCi  = self.avgAngCi;
                end
                
                [self.rayleighP,self.rayleighZ] = circ_rtest( ...
                    axTrans(binCentersRad),histCounts,binSizeRad);
                pRayl = self.rayleighP;
                zRayl = self.rayleighZ;
                
                if isFrequency % convert to counts per second
                    histData = histData / binSizeSec; end
            end
            self.histData = histData;
            % correlation analysis
            [self.corrAnR,self.corrAnP] = ...
                circ_corrcl(axTrans(binCentersRad),histData(:,1));
            corrAnR = self.corrAnR;
            corrAnP = self.corrAnP;
            
            %% initialize figure, set visual properties
            % axis labels
            if isFrequency, self.axisLabel = 'Counts per second';
            else,           self.axisLabel = 'Counts per bin';   end
            % initialize theta-label
            self.thetaLabel = text; % empty
            
            if ~isempty(ax) && isvalid(ax), figH = ax.Parent;
            else,                           figH = gcf; end
            
            figH.Visible = 'off';
            self.figH = figH;
            set(0,'currentfigure',figH);
            polarplot(0);hold on
            set(figH,'color',[1,1,1]) % white background
            self.polarAxs = gca;
            polarAxs = self.polarAxs;
            
            % self-reference in property for hyper-redundancy (this is actually quite
            % handy if you want to retrieve the CircHist-object from a figure)
            polarAxs.UserData.circHistObj = self;
            
            self.thetaLim = thetaLim;
            polarAxs.ThetaAxisUnits = 'degree';
            polarAxs.ThetaZeroLocation = 'top';
            polarAxs.Tag = 'Polar';
            
            polarAxs.Units = 'normalized';
            % decrease axes-dimensions to provide room for the scale bar
            polarAxs.Position(3) = polarAxs.Position(3) * 0.95;
            polarAxs.Position(4) = polarAxs.Position(4) * 0.95;
            self.scaleBarSide = scaleBarSide;
                        
            lineSp = '-'; %continuous lines
            self.barWidth = adjustSlope * binSizeDeg + 4; %bar width
            self.stdWidth = self.barWidth / 3;
            lineWAvgAng = adjustSlope * binSizeDeg/10 + 2;
            lineWR = lineWAvgAng * 1.5;
            
            %% draw bars and whiskers
            self.drawBars
            %%
            circR = max(rlim); %radius of plot in data units
            
            %calculate baseline shift from center; depends on bin size for aesthetic reasons
            if ~isnan(baseLineOffset)
                %shift baseline by the fraction of the full radius and this value
                baseLineOffsetDiv = adjustSlope * binSizeDeg + baseLineOffset;
                baseLineOffset = -round(circR/baseLineOffsetDiv); %baseline shift in data units
            else, baseLineOffset = 0;
            end
            
            self.drawBaseLine
            
            %shift radius baseline
            polarAxs.RLim = [baseLineOffset,circR];
            
            %% average angle and r
            % based on AREAXIALDATA, the average-angle and r lines are axes in the
            % histogram, or lines with the average angle as direction
            avgAngRad = deg2rad(avgAng);
            if ~isempty(avgAng) && drawAvgAng %plot average angle
                if areAxialData
                    thetaAvgAng = [avgAngRad,avgAngRad+pi];
                    rhoAvgAng = [circR,circR];
                else
                    thetaAvgAng = [avgAngRad,avgAngRad];
                    rhoAvgAng = rlim;
                end
                self.avgAngH = polarplot(self.polarAxs,thetaAvgAng,rhoAvgAng,lineSp...
                    ,'lineWidth',lineWAvgAng,'color',colorAvgAng,'Tag','avgAng');
            end
            if ~isempty(r) && ~isempty(avgAng) && drawR
                % make vector length relative to plot radius (after shift)
                rNorm = r * range(rlim) + baseLineOffset;
                if areAxialData
                    thetaR = [avgAngRad,avgAngRad+pi];
                    rohR = [rNorm,rNorm];
                else
                    thetaR = [avgAngRad,avgAngRad];
                    rohR = [min(rlim),rNorm];
                end
                self.rH = polarplot(self.polarAxs,thetaR,rohR,lineSp...
                    ,'lineWidth',lineWR,'color',colorR,'Tag','r');
            end
            
            % distance between confidence-interval line and plot border in percent
            avgAngCiWhiskOffset = 0.02;
            self.polarAxs.UserData.avgAngCiWhiskOffset = avgAngCiWhiskOffset;

            avgAngCiRad = deg2rad(avgAngCi);
            if ~isempty(avgAngCi) && drawAvgAngCi
                self.avgAngCiH = gobjects(0);
                avgAngCiPlotArgs = {'lineStyle',lineSp,'color',colorAvgAngCi ...
                    ,'Clipping','off'};
                
                whiskWidthEnd = self.stdWidth * 0.7; % width of whisker-endings
                thetaStepN = ceil(avgAngCi / 0.2); % plot line in 0.2�-steps
                thetaAvgAngCi = linspace( ...
                    avgAngRad - avgAngCiRad,avgAngRad + avgAngCiRad,thetaStepN);
                
                rhoAvgAngCi = max(rlim) + avgAngCiWhiskOffset * range(rlim);
                avgAngCiWhiskLen = avgAngCiWhiskOffset*2/3 * range(rlim);
                rhoAvgAngCiWhiskEnd = ...
                    [rhoAvgAngCi + avgAngCiWhiskLen, rhoAvgAngCi - avgAngCiWhiskLen];
                rhoAvgAngCi = repmat(rhoAvgAngCi,thetaStepN,1);
                
                self.avgAngCiH(1,1) = polarplot(thetaAvgAngCi,rhoAvgAngCi,avgAngCiPlotArgs{:} ...
                    ,'lineWidth',self.stdWidth,'Tag','avgAngCiWhisk');
                self.avgAngCiH(1,2) = polarplot([thetaAvgAngCi(1),thetaAvgAngCi(1)] ...
                    ,rhoAvgAngCiWhiskEnd,avgAngCiPlotArgs{:} ...
                    ,'lineWidth',whiskWidthEnd,'Tag','avgAngCiWhiskEnd');
                self.avgAngCiH(1,3) = polarplot([thetaAvgAngCi(end),thetaAvgAngCi(end)] ...
                    ,rhoAvgAngCiWhiskEnd,avgAngCiPlotArgs{:} ...
                    ,'lineWidth',whiskWidthEnd,'Tag','avgAngCiWhiskEnd');
                
                if areAxialData % plot again with mirrored angles
                    thetaAvgAngCi = thetaAvgAngCi + pi;
                    self.avgAngCiH(2,1) = polarplot(thetaAvgAngCi,rhoAvgAngCi ...
                        ,'lineWidth',self.stdWidth,avgAngCiPlotArgs{:},'Tag','avgAngCiWhisk');
                    self.avgAngCiH(2,2) = polarplot([thetaAvgAngCi(1),thetaAvgAngCi(1)] ...
                        ,rhoAvgAngCiWhiskEnd,avgAngCiPlotArgs{:} ...
                        ,'lineWidth',whiskWidthEnd,'Tag','avgAngCiWhiskEnd');
                    self.avgAngCiH(2,3) = polarplot([thetaAvgAngCi(end),thetaAvgAngCi(end)] ...
                        ,rhoAvgAngCiWhiskEnd,avgAngCiPlotArgs{:} ...
                        ,'lineWidth',whiskWidthEnd,'Tag','avgAngCiWhiskEnd');
                end
            end
            
            %% edit axes
            polarAxs.Color = 'none'; %no background
            labels = polarAxs.ThetaTickLabel;
            % add �-sign to labels
            if ~strcmp(labels{1}(end),'�') % in the case of subsequent plotting into the same axis
                for n = 1:numel(labels)
                    labels{n} = [labels{n},'�'];end
                polarAxs.ThetaTickLabel=labels;
            end
            polarAxs.ThetaAxis.FontSize = fontSize;
            
            polarAxs.LineWidth = 1;
            polarAxs.GridColor = 'k';
            polarAxs.GridAlpha = 0.5;
            % nifty way to dynamically create minor grid lines in 10� spacing, skipping
            % the major grid lines. No one will ever want to understand this.
            minorGrid = 10:10:350;
            minorGrid = minorGrid(logical(mod(minorGrid,polarAxs.ThetaTick(2))));
            polarAxs.ThetaAxis.MinorTickValues = minorGrid;
            polarAxs.ThetaMinorGrid = 'on';
            polarAxs.MinorGridColor = 'k';
            polarAxs.MinorGridAlpha = 0.5;
            
            % draw white disk as background in the center
            self.drawWhiteDisk;
            
            %% title
            details = sprintf(['N = %u , avgAng = %.2f�\\pm%.2f�, r = %.4f\np_{Rayl} = %.3f, '...
                ,'Z_{Rayl} = %.4f, p_{C. an.} = %.3f, R�_{C. an.} = %.3f']...
                ,nSamples,avgAng,avgAngCi,r,pRayl,zRayl,corrAnP,corrAnR*corrAnR);
            title(details,'FontSize',9);
            %%
            colormap white %for axis appearance
            polarAxs.RTickLabel = [];
            
            if isempty(ax) % adjust figure-window size
                figH.OuterPosition(1) = figH.OuterPosition(1)/2; 
                figH.OuterPosition(2) = figH.OuterPosition(2)/2;
                figH.OuterPosition(3) = figH.OuterPosition(3)*1.2;
                figH.OuterPosition(4) = figH.OuterPosition(4)*1.5;
            end
            
            self.drawScale;
            figH.SizeChangedFcn = {@self.redrawScale,self};
            figH.Visible = 'on';
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
                currVal = self.histData(n,1);
                if currVal == 0 || currVal <= lowerLim % skip bar-drawing
                    continue, end
                currStd = self.histData(n,2);
                
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
            if ~isempty(self.avgAngH) && isvalid(self.avgAngH)
                uistack(self.avgAngH,'top');
            end
            if ~isempty(self.rH) && isvalid(self.rH)
                uistack(self.rH,'top');
            end
        end
        %% drawBaseLine
        function drawBaseLine(self)
            % Draws a base line at 0. This function should be called before drawWhiteDisk
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
            
            sclLeft = strcmp(self.scaleBarSide,'left');

            %default font size. This line should be completely useless, but for some
            %reason, there are strange effects without it. The fact that the bar-drawing
            %is executed via SIZECHANGEDFCN seems to interrupt the link between the font
            %properties of the bar and the font properties of the parent axes, messing
            %something up.
            fontsz = 13; %#ok

            initialDraw = isempty(scl);
            if initialDraw % if the bar is drawn for the first time, use predefined label and font size
                hasLabel = true;
                label = self.axisLabel;
                fontsz = self.fontSize;
            else % read out (assumingly predefined) properties and use these for the next drawing
                hasLabel = ~isempty(scl.Label);
                if hasLabel
                    label = scl.Label.String;
                    fontsz = scl.Label.FontSize;
                else
                    label = '';
                    fontsz = self.fontSize;
                end
                oldLineWidth = scl.LineWidth;
            end
            
            delete(scl); %deletes the previous bar
            
            scl = colorbar('Location','east');
            self.scaleBar = scl; % re-assign
            
            % with this link, the label font-name is changed with the corresponding
            % axes-property using set(gca,'FontName',fontName), which is the default
            % behavior of colorbar labels.
            pAx.UserData.fontNameLink = linkprop([pAx,scl],'FontName');
            
            sclUnitsOld = scl.Units;
            scl.Units = 'pixels';
            polarAxsUnitsOld = pAx.Units;
            pAx.Units = 'pixels';
            
            sclWidth = 0.015;
            
            polarPos = pAx.Position; %position property = [left,bottom,width,height]
            polarLeft   = polarPos(1);
            polarBot    = polarPos(2);
            polarWidth  = polarPos(3);
            polarHeight = polarPos(4);
            polarRight  = polarLeft + polarWidth;                        
            
            % this works because the plot-area is always a square, thus the lower value of
            % the dimensions equals the plot-circle diameter
            polarPlotDiam = min([polarWidth,polarHeight]); % diameter of plot
            if polarHeight > polarPlotDiam
                % distance between upper plot-edge and upper figure-edge
                  polarHeightOffset = (polarHeight - polarPlotDiam)/2;
            else, polarHeightOffset = 0;
            end
            
            sclHeight = polarPlotDiam/2; % diameter/2 -> bar spans from center to top
            % adjust scale-height to theta-offset (scale starts at 0)
            lowerLim = pAx.RLim(1);
            if lowerLim < 0, offset = sclHeight * abs(lowerLim)/range(pAx.RLim);
            else,            offset = 0;
            end
            sclHeight = sclHeight - offset;
            % calculate bottom-position of scale bar
            sclBot = polarBot + polarHeight - polarHeightOffset - sclHeight;
            
            % use the margins between plot-area and axis-labels to adjust the
            % scale-bar left-position
            aBit = 10; % just a tiny little bit of space
            if sclLeft
                polarEdge = polarLeft;  sideMarg = pAx.TightInset(1) - sclWidth; sclSideSign = -1;
            else
                polarEdge = polarRight; sideMarg = pAx.TightInset(3) + sclWidth; sclSideSign = +1;
            end
            botMarg = pAx.TightInset(2); % margin at bottom
            % a margin > 0 means that there is no space beyond the axis labels, so the
            % scale-bar must be located at the border of the labels. A side margin = 0
            % means that there is (excess) space, but because the plot-area automatically
            % resizes to be a square, the bottom margin gives information about the actual
            % side margin, so the scale-bar is placed at the side border of the plot plus
            % the bottom margin.
            if sideMarg > 0, sclLeft = polarEdge + sclSideSign * (sideMarg + aBit);
            else
                sclLeft = pAx.OuterPosition(3)/2 ... % mid-line of axes-area
                    + (sclSideSign * polarPlotDiam)/2 ... % +/- radius of plot
                    + sclSideSign * (botMarg + aBit); % +/- (margin + a bit more)
            end
            scl.Position = [sclLeft,sclBot,sclWidth,sclHeight];

            if hasLabel, scl.Label.String = label; end
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
        %% draw white disk in the middle to obscure the axis-lines
        function drawWhiteDisk(self)
            if ~isempty(self.whiteDiskH), delete(self.whiteDiskH); end
            
            offset = self.polarAxs.RLim(1);
            if offset >= 0, return, end % nothing to obscure
            
            pAx = self.polarAxs;
            rhoWhite = [offset,0];
            thetaLims = pAx.ThetaAxis.Limits;
            if strcmp(pAx.ThetaAxisUnits,'radians'), thetaLims = rad2deg(thetaLims); end
            % 2� steps between white lines
            thetaWhite = deg2rad(thetaLims(1):2:thetaLims(2));
            for t = 1:numel(thetaWhite)
                currAng = thetaWhite(t);
                polarplot(pAx,[currAng,currAng],rhoWhite...
                    ,'linestyle','-','color','w','linewidth',3,'Tag','whiteDisk');
            end
            self.whiteDiskH = findobj(pAx,'Tag','whiteDisk');
            uistack(self.whiteDiskH,'bottom'); % move lines to lowest graphical layer
        end
        %% change scale limits
        function setRLim(self,limits)
            %setRLim Change scale limits specified by the two-element vector LIMITS ==
            %[lower,upper]. Get the current limits by calling rlim.
            %
            % circHistObj.setRLim(limits); where LIMITS == [lower,upper]
            %
            rlim(limits); %change limits
            
            % update line data
            self.drawBars
            if isvalid(self.avgAngH)
                if self.areAxialData, self.avgAngH.RData(:) = limits(2);
                else,                 self.avgAngH.RData = limits;       end
            end
            if isvalid(self.rH)
                rNorm = self.r * range(limits) + limits(1);
                if self.areAxialData, self.rH.RData(:) = rNorm;
                else,                 self.rH.RData = [limits(1),rNorm]; end
            end
            if all(isvalid(self.avgAngCiH))
                oldRho = self.avgAngCiH(1,1).RData(1);
                newRho = limits(2) + self.polarAxs.UserData.avgAngCiWhiskOffset * range(limits);
                avgAngCiRhoOffset = newRho - oldRho;
                arrayfun(@(h)set(h,'RData',h.RData + avgAngCiRhoOffset),self.avgAngCiH);
            end
            
            if isvalid(self.thetaLabel)
                self.setThetaLabel(self.thetaLabel.String ...
                    ,self.thetaLabel.UserData.location);
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
            self.drawWhiteDisk; % adjust background
        end
        %% change angle-axis limits
        function set.thetaLim(self,limits)
            %ThetaLim Changes angle-axis limits to the specified two-element vector LIMITS
            %== [lower,upper].
            %
            % circHistObj.ThetaLim = limits; where LIMITS == [lower,upper]
            %
            self.polarAxs.ThetaAxis.Limits = limits;
            self.thetaLim = limits;
        end
        %% set theta-axis label
        function setThetaLabel(self,txt,location)
            %setThetaLabel Add (or update) a label to the theta-axis, label-text specified
            % by TXT, location specified by LOCATION, which may be 'bottomleft'(default if
            % omitted), 'bottomright', 'topleft' or 'topright'.
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
            
            % based on THETADIR and THETAZEROLOCATION, the label-theta angle needs to be
            % adjusted so the label is in the specified corner
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
            self.thetaLabel.UserData.location = location;
            % make label-color be linked to theta-axis color (default behavior for labels)
            % for some reason, the link only works when THETALABEL.COLOR is changed, not
            % the other way around ...
            self.polarAxs.UserData.thetaLabelColorLink = ...
                linkprop([self.polarAxs.ThetaAxis,self.thetaLabel],'Color');
        end
        %% change scale bar side
        function set.scaleBarSide(self,side)
            %   obj.scaleBarSide(side); where SIDE is either 'left' or 'right'
            assert(any(strcmpi(side,{'left','right'})) ...
                ,'Property SCALEBARSIDE must be ''left'' or ''right''.');
            if strcmpi(side,self.scaleBarSide), return; end % nothing to do
            
            pAx = self.polarAxs;
            isLeft = strcmpi(side,'left');
            % initially, the plot is not positioned in the center of the axes, but
            % slightly right, and there seems to be no way to center it, so the first
            % adjustment needs to be different than subsequent adjustments
            if isempty(self.scaleBarSide) % first call during construction
                if isLeft, offset = 0.05; else, offset = -0.05; end
            else
                if isLeft, offset = 0.1;  else, offset = -0.1; end
            end
            
            self.scaleBarSide = side;
            unitsOld = pAx.Units;
            pAx.Units = 'normalized';            
            pAx.Position(1) = pAx.Position(1) + offset;
            pAx.Units = unitsOld;
            if ~isempty(self.scaleBar) && isvalid(self.scaleBar), self.drawScale; end
        end
        %% change scale bar label
        function set.axisLabel(self,newLabel)
            self.axisLabel = newLabel;
            if ~isempty(self.scaleBar) && isvalid(self.scaleBar)
                self.scaleBar.Label.String = newLabel; end
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
            % Change the standard-deviation main-line width. Since each whisker consists
            % of two line-objects with different widths (the "main" line and the ending),
            % the width of the ending is scaled proportionally to the width-change of the
            % "main" line.
            %
            %   circHistObj.stdWidth(2);
            
            oldWidth = self.stdWidth;
            self.stdWidth = width;
            
            %happens at object construction
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
            else, print(self.figH,fileName,'-dpdf','-fillpage','-painters');
            end
        end
        %% save to png
        function toPng(self,fileName,resol)
            %toPng Save histogram as (FILENAME).png at the optionally specified resolution
            %(default = 90 dpi). Specify RESOL as a string of the pattern '-r90'.
            %
            %   obj.toPng(filename);
            %   obj.toPng(filename,resol);
            
            if nargin < 3, resol = '-r90'; end
            
            if exist('toPng','file') % call custom function if available
                toPng(self.figH,fileName,resol);
            else, print(self.figH,fileName,'-dpng','-opengl',resol);
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
            
            if isvalid(self.avgAngH), uistack(self.avgAngH,'top'); end
            if isvalid(self.rH),      uistack(self.rH,'top'); end
            
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