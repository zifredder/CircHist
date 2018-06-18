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
    %   Requirements: CircStat toolbox by Philipp Berens & Marc J. Velasco
    %                 (http://www.jstatsoft.org/v31/i10 or
    %                 http://www.mathworks.com/matlabcentral/fileexchange/10676 or
    %                 https://github.com/circstat/circstat-matlab)
    %
    %   Usage:  CircHist(data);
    %           CircHist(data,edges);
    %           CircHist(data,nBins);
    %           CircHist(___,Name,Value);
    %           obj = CircHist(___);
    %
    %
    %   ---Required Input Argument:
    %       data        Either one vector of angle samples (distribution) in degree (!), a
    %                   cell array of such samples, or a two-column matrix of N length
    %                   where N is the number of bins, the first column contains the
    %                   counts per bin and the second column the standard deviations (or,
    %                   generally, error-whisker lenghts). If the second column is
    %                   omitted, the standard deviations are considered zero. According to
    %                   the input data, specify the DATATYPE property accordingly.
    %
    %
    %  ---Optional Input Argument:
    %       edges       Edges of histogram bins, either specified by a vector of bin-edges
    %                   in degree, e. g. [0:20:360] for 20-degree bins, or by a scalar
    %                   positive integer specifying the number of bins between 0 and 360
    %                   degree, e. g. 18 for 20-degree bins. If omitted: For
    %                   already-binned data (DATATYPE == 'histogram'), the bins are
    %                   deduced from the number of data points. For distribution-data, the
    %                   HISTCOUNTS default binning-algorithm is used to calculate
    %                   bin-widths (only available in MATLAB R2014b and newer; see
    %                   HISTCOUNTS). Note, however, that with this procedure, the bins
    %                   might not span the whole circle and consider possible implications
    %                   for the statistics. For distribution-data that contain >1 sample,
    %                   the bin width is determined based on the sample with the most data
    %                   points.
    %
    %
    %   ---Optional Name-Value Pair Input:
    %   dataType        'distribution'(default)/'histogram'. Type of input data:
    %                   'distribution' treats the input data as distributions of angles
    %                   (degree); the data of each input vector are binned inside the
    %                   specified edges, each bin is averaged and the standard deviation
    %                   is calculated. 'histogram' treats the data as already-binned data
    %                   (counts per bin); if the input matrix has a second column, it is
    %                   taken as the standard-deviation values for each bar. If you want
    %                   to plot PSTH-data (peri-stimulus-time-histogram) where your
    %                   bin-heights should be converted to frequencies, specify HISTTYPE
    %                   and BINSIZESEC accordingly. In the 'histogram'-case, EDGES does
    %                   not need to be specified, as it results from the number of data
    %                   points.
    %
    %   histType        'count'(default)/'frequency'. Based on this, the plotted values
    %                   are either treated as counts per bin or converted to counts per
    %                   second (feasible for PSTH-data). Note that BINSIZESEC needs to be
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
    %                   obj.scaleBar. Change the side after plotting by changing
    %                   OBJ.SCALEBARSIDE to 'left' or 'right'.
    %
    %   drawAvgAng      True(default)/false, plots average angle as a line.
    %
    %   avgAng          Numeric value (degree) of the average angle. If specified, AVGANG
    %                   is not calculated from the data and this value is taken and
    %                   plotted instead.
    %
    %   drawAvgAngCi    True(default)/false, plots 95 % confidence interval as a circle
    %                   segment outside of the plot.
    %
    %   avgAngCi        Numeric value (degree) of 95 % confidence interval (AVGANG +/-
    %                   AVGANGCI) to plot; same behavior as for AVGANG.
    %
    %   drawR           True(default)/false, plots resultant vector length (r) as a black
    %                   overlay bar on top of the average angle line. The black bar's
    %                   length equals the r-value in percent of the average-angle-line
    %                   length.
    %
    %   r               Numeric value of the resultant vector length (between 0 and 1);
    %                   same behavior as for AVGANG.
    %
    %   baseLineOffset  Numeric value (default == 2) specifying the factor that scales the
    %                   size of the offset in the center of the plot. The default value
    %                   produces nice results; lower the value (negative values allowed)
    %                   to increase the size, increase the value to decrease the size.
    %                   Specify it as NaN to produce an offset of 0. (This functionality
    %                   might need improvement, maybe just implement it as a percent value
    %                   of the plot diameter; as of now, it depends on bin-size and
    %                   ADJUSTSLOPE)
    %
    %   adjustSlope     Slope that defines how strongly optical properties such as bar
    %                   width scale with bin-size; default == 0.3.
    %
    %   ax              Axes handle to plot diagram into; becomes POLARAXS property. Note
    %                   that the referenced axes must be a POLARAXES. (Experimental
    %                   feature, working in principle, but the scale sometimes
    %                   misbehaves).
    %
    %   colorBar        Color of bars (default == [0 .45 .74]; (standard MATLAB blue)).
    %                   Can be adjusted after plotting by changing the property of the
    %                   same name.
    %   colorStd        Color of standard-deviation lines (default == 'k'). Can be
    %                   adjusted after plotting.
    %   colorAvgAng     Color of average-angle line (default == [.85 .33 .1]; (orange)).
    %                   Can be adjusted after plotting.
    %   colorR          Color of r line (default == 'k'). Can be adjusted after plotting.
    %   fontSize        Font size of axis labels (default == 13). Corresponds to the
    %                   FONTSIZE property of the POLARAXES. Can be adjusted after
    %                   plotting.
    %
    %
    %   ---Output:
    %       obj         CircHist object containing the reference to the created diagram.
    %                   Use this to adjust the diagram after plotting (see methods and
    %                   properties below).
    %
    %
    %   ---Methods and Noteworthy Properties:
    %       setRLim([lower,upper])      Change rho-axis limits (usage:
    %                                   obj.setRLim([lower,upper])) (get current limits by
    %                                   calling RLIM). Do not use RLIM to change the
    %                                   limits because this does not update diagram
    %                                   elements that depend on the diagram boundaries.
    %       setThetaLabel(txt,location) Adds (or updates) a label saying TXT outside of
    %                                   the plot at the location specified by LOCATION,
    %                                   which may be one of the following: 'topleft',
    %                                   'topright', 'bottomleft'(default if omitted),
    %                                   'bottomright'. (usage:
    %                                   obj.setThetaLabel('Direction','bottomright') ).
    %                                   Specify TXT as a cell array of characters to add
    %                                   line breaks. Access the created text-object via
    %                                   obj.thetaLabel.
    %       drawArrow(theta,rho)        Draws arrows from the plot center at specified
    %                                   angles with specified lengths; usage:
    %                                   obj.drawArrow(theta,rho,Name,Value), where THETA
    %                                   is a scalar or vector of angles in degree and RHO
    %                                   is a vector of the same size specifying the arrow
    %                                   length in data units. RHO may be specified as []
    %                                   to have the arrows end at the plot edge.
    %                                   Name-Value pairs as accepted by ANNOTATION (see
    %                                   DOC ANNOTATION) can be passed to adjust arrow
    %                                   properties.
    %       updateArrows                Updates arrow-positions. Call this in case the
    %                                   arrows misbehave: obj.updateArrows;
    %       scaleBarSide                Change the side ('left'/'right') of the
    %                                   rho-axis scale bar (usage: obj.scaleBarSide =
    %                                   'right').
    %       colorBar                    Change histogram-bar color (usage: obj.colorBar =
    %                                   newColor).
    %       barWidth                    Change histogram-bar width (usage: as above).
    %       colorStd                    Change standard-deviation-whisker color (usage: as
    %                                   above).
    %       stdWidth                    Change standard-deviation-whisker width (usage: as
    %                                   above).
    %       toPdf(fileName)             Save as pdf-file with the specified file name.
    %       toPng(fileName,resol)       Save as png-file with the specified file name at
    %                                   the optionally specified resolution (e.g.
    %                                   '-r300').
    %       drawCirc(rho,lineSpec)      Draws a circle in the plot with radius RHO, line
    %                                   appearance optionally specified by LINESPEC (see
    %                                   DOC LINESPEC). Optionally returns the handle to
    %                                   the graphics object. Example usage:
    %                                   h = obj.drawCirc(15,'-g','LineWidth',4);
    %       drawScale                   Updates the rho-axis scale; call this, e. g.,
    %                                   after the theta-axis ticks have changed (usage:
    %                                   obj.drawScale).
    %       exampleCircHist             Creates example diagrams.
    %
    %
    %   ---Notes and Instructions:
    %   * For more details on a method, read its doc page, e. g.: doc CircHist.drawArrow
    %   * To change the rho-axis limits, use obj.setRLim([lower,upper]); to change the
    %     angle-axis limits, use obj.thetalim([lower,upper]) or just thetalim(_).
    %   * To change the scale-label, use obj.axisLabel = 'my label'.
    %   * To change visual properties, either use the Name-Value pairs for the constructor
    %     as specified below, or access the graphics objects via their handles that are
    %     stored as properties of the CircHist object, such as obj.polarAxs for the
    %     coordinate system (font size, line width, etc.) and obj.scaleBar for the scale
    %     bar (line thickness etc.). The scale bar is drawn anew each time the
    %     figure-window's size is changed; note that this may not work flawlessly and
    %     assure that the scale bar matches the coordinate-grid after changing the figure
    %     size. Consider calling obj.drawScale after changing properties to update the
    %     scale bar. Remove the scale bar with delete(obj.scaleBar).
    %   * To adjust the bars, standard-deviation whiskers, average-angle line and r line
    %     after plotting, use their handles which are saved as properties. Access them
    %     using dot notation, e. g., h = obj.avgAngH.lineWidth, and delete them using
    %     delete(h) (not: obj.avgAngH.lineWidth = []).
    %   * Note that the constructed CircHist-object handle is stored in
    %     obj.polarAxs.UserData.circHistObj which can be accessed via the POLARAXES it has
    %     been plotted into. This may be useful if you plot a series of histograms and
    %     forget to store the CircHist-objects.
    %   * Each standard-deviation whisker consists of a long line representing the
    %     magnitude of the deviation and a very short, thick line that marks the tip of
    %     the deviation-line. Both line types are comprised in the handle-array obj.stdH;
    %     however, they can be separately accessed by calling
    %     findobj(obj.polarAxs,'Tag',TYPETAG), where TYPETAG is either 'stdWhisk' for the
    %     "main" line or 'stdWhiskEnd' for the tips.
    %   * Analogously, access the 95 % confidence-interval line(s) by accessing the
    %     handles in obj.avgAngCiH or by using findobj(obj.polarAxs,'Tag',TYPETAG), where
    %     TYPETAG is either 'avgAngCiWhisk' or 'avgAngCiWhiskEnd'.
    %   * Consider creating a new figure for each histogram because there may be
    %     side-effects regarding the axis and the scale bar if the same axes-object or
    %     figure-window are used. When plotting multiple histograms into the same figure
    %     by using the 'ax' Name-Value pair, expect issues such as overlapping graphical
    %     elements.
    %   * To change the angle-axis units into radians, call
    %     obj.polarAxs.ThetaAxisUnits = 'radians';
    %   * The axis-grid-lines in the center of the plot are obscured by white bars for
    %     cosmetic reasons. To access these bars, e. g. to change their color, refer to
    %     obj.whiteDiskH.
    %   * To change the layer order of the graphical objects, use UISTACK, e. g.:
    %     uistack(obj.whiteDiskH,'top'). Note, however, that arrows as created from
    %     DRAWARROW are always on top of the histogram objects because they are anchored
    %     to the figure and thus are independent from objects in the histogram axes.
    %   * To get neat auto-completion for Name-Value pairs to work, add
    %     functionSignatures.json to the same directory as the @CircHist folder, or add
    %     its contents to your existing file; see also:
    %     https://mathworks.com/help/matlab/matlab_prog/customize-code-suggestions-and-completions.html
    %
    %
    %  ---Author: Frederick Zittrell
    %
    % See also CircStat polaraxes polarplot thetatickformat thetaticks rticks thetalim
    % uistack

    properties (SetAccess = immutable)
        data            % Required input: Data in degree.
        edges           % Required input: Histogram-bin edges or number of histogram-bins.
                
        dataType        % Optional input; 'distribution'(default)/'histogram'.
        histType        % Optional input; 'frequency'(default)/'count'.
        binSizeSec      % Optional input; Width of bins (s).
        drawAvgAng      % Optional input; True(default)/false.
        avgAng          % Optional input; Numeric value of the average angle.
        drawAvgAngCi    % Optional input; True(default)/false.
        avgAngCi        % Optional input; Numeric value of the 95 % confidence interval bounds (average +/ bounds).
        drawR           % Optional input; True(default)/false.
        r               % Optional input; Numeric value of the resultant vector length.
        baseLineOffset  % Optional input; scaling factor for the size of the offset in the middle of the diagram.
        adjustSlope     % Optional input; Slope for scaling of visual properties with bin size.
        areAxialData    % Optional input; True/false(default).
        
        polarAxs        % Handle to POLARAXES that containes the histogram. Change visual properties such as line width of the axes here.
        figH            % Handle to figure in which the histogram is plotted.
        
        histData        % Histogram data as plotted; 1st column average counts, 2nd column standard deviations.
        rayleighP       % P-value of Rayleigh test of uniformity.
        rayleighZ       % Z-value of Rayleigh test of uniformity.
        corrAnP         % P-value of correlation analysis.
        corrAnR         % R-value of correlation analysis (square this to get the coefficient of determination).
    end
    
    properties
        scaleBar        % Handle of scale bar. Use to access visual properties.
        axisLabel       % Label of scale bar. Can be set using obj.axisLabel = 'new label'.
        scaleBarSide    % Side of the scale bar, either 'left'(default) or 'right'.
        thetaLabel      % Label of the degree-axis (TEXT object).
        
        avgAngH         % Handle to the average-angle line.
        avgAngCiH       % Handle to confidence interval lines.
        rH              % Handle to the r line.
        barH            % Array of handles to the bars (LINE objects).
        stdH            % Array of handles to the standard-deviation whiskers (LINE objects).
        arrowH          % Array of handles to drawn arrows (ANNOTATION objects).
        
        colorBar        % Optional input; Color of bars (default = [0 .45 .74]; (standard MATLAB blue)).
        colorStd        % Optional input; Color of standard-deviation lines (default = 'k').
        colorAvgAng     % Optional input; Color of average-angle line (default = [.85 .33 .1]; (orange)).
        colorR          % Optional input; Color of r line (default = 'k').
        colorAvgAngCi   % Optional input; Color of confidence interval line (default = 'k').
        fontSize        % Optional input; Font size of axis labels (default = 13).
        
        barWidth        % Width of bars. Change this property to adjust the bar width after plotting.
        stdWidth        % Width of standard-deviation whiskers. For adjustment after plotting.
        
        whiteDiskH      % Handles to white bars that obscure the center of the axis.
    end
    
    methods
        %% constructor
        function self = CircHist(data,varargin)   
            %CircHist constructor
            %       obj = CircHist(data[,edges][,Name,Value])
            
            %% empty-object case
            if nargin == 0, return; end
            
            %% validate and parse input
            if exist('circ_axial','file') ~= 2
                error('CircStat toolbox not found, please install (see doc CircHist).');
            end
            
            if exist('isColorSpec','file')
                  validColor = @isColorSpec; % custom function, may not be present
            else, validColor = @(x)validateattributes(x,{'numeric','char'},{'vector'});
            end
            
            % HISTCOUNTS was introduced in MATLAB R2014b; use HISTC for older versions
            hasHistcounts = exist('histcounts','file') || exist('histcounts','builtin');
            
            % validates input to be scalar and either logical or a numeric 0 or 1
            validLogical01 = @(tf) isscalar(tf) ...
                && (islogical(tf) || isnumeric(tf) && ismember(tf,[0,1]));
            validScalarNum = @(N) isscalar(N) && isnumeric(N);
            
            % default values
            def.edges          = [];
            def.dataType       = 'distribution';
            def.histType       = 'count';
            def.binSizeSec     = [];
            def.drawAvgAng     = true;
            def.avgAng         = [];
            def.drawAvgAngCi   = true;
            def.avgAngCi       = [];
            def.drawR          = true;
            def.r              = [];
            def.baseLineOffset = 2;
            def.adjustSlope    = 0.3;
            def.axialData      = false;
            def.axes           = [];
            def.colorBar       = [0 .45 .74]; % matlab blue
            def.colorStd       = 'k';
            def.colorAvgAng    = [.85 .33 .1];
            def.colorR         = 'k';
            def.colorAvgAngCi  = 'k';
            def.fontSize       = 13;
            def.scaleBarSide   = 'left';
            
            pr = inputParser;
            pr.FunctionName = 'CircHist';
            
            addRequired(pr,'data' ...
                ,@(x) validateattributes(x,{'numeric','cell'},{'nonempty'}));
            addOptional(pr,'edges',def.edges ...
                ,@(x) validateattributes(x,{'numeric'},{'vector','increasing'}));
            
            addParameter(pr,'dataType',def.dataType ...
                ,@(str) any(strcmpi(str,{'distribution','histogram'})));
            addParameter(pr,'histType',  def.histType ...
                ,@(str) any(strcmpi(str,{'count','frequency'})));
            addParameter(pr,'binSizeSec',def.binSizeSec ...
                ,@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
            
            addParameter(pr,'drawAvgAng'    ,def.drawAvgAng,validLogical01);
            addParameter(pr,'avgAng'        ,def.avgAng,validScalarNum);
            addParameter(pr,'drawAvgAngCi'  ,def.drawAvgAngCi,validLogical01);
            addParameter(pr,'avgAngCi'      ,def.avgAngCi,validScalarNum);
            addParameter(pr,'drawR'         ,def.drawR,validLogical01);
            addParameter(pr,'r'             ,def.r,validScalarNum);
            addParameter(pr,'baseLineOffset',def.baseLineOffset,validScalarNum);
            addParameter(pr,'adjustSlope'   ,def.adjustSlope,validScalarNum);
            addParameter(pr,'areAxialData'  ,def.axialData,validLogical01);
            addParameter(pr,'ax'            ,def.axes ...
                ,@(x) validateattributes(x,{'matlab.graphics.axis.PolarAxes'},{'scalar'}));
            
            addParameter(pr,'scaleBarSide' ...
                ,def.scaleBarSide,@(str) any(strcmpi(str,{'left','right'})));
            addParameter(pr,'colorBar',   def.colorBar,   validColor);
            addParameter(pr,'colorStd',   def.colorStd,   validColor);
            addParameter(pr,'colorAvgAng',def.colorAvgAng,validColor);
            addParameter(pr,'colorR',     def.colorR,     validColor);
            addParameter(pr,'colorAvgAngCi',def.colorAvgAngCi,validColor);
            addParameter(pr,'fontSize',   def.fontSize,   validScalarNum);
            
            parse(pr,data,varargin{:});
            
            data            = pr.Results.data;
            edges           = pr.Results.edges;
            if isscalar(edges), edges = 0 : (360/edges) : 360; end
            
            self.dataType   = pr.Results.dataType;
            areDistribData  = strcmpi(self.dataType,'distribution');
            binSizeSec      = pr.Results.binSizeSec;
            self.binSizeSec = binSizeSec;
            self.histType   = pr.Results.histType;
            isFrequency     = strcmpi(self.histType,'frequency');
            
            self.drawAvgAng      = pr.Results.drawAvgAng;
            drawAvgAng           = self.drawAvgAng;
            self.avgAng          = pr.Results.avgAng;
            avgAng               = self.avgAng;
            self.drawAvgAngCi    = pr.Results.drawAvgAngCi;
            drawAvgAngCi         = self.drawAvgAngCi;
            self.avgAngCi        = pr.Results.avgAngCi;
            avgAngCi             = self.avgAngCi;
            
            self.drawR          = pr.Results.drawR;
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
            
            self.arrowH         = gobjects(0); % empty handle array
            %% validate HISTTYPE-BINSIZESEC combination
            if isempty(binSizeSec) && isFrequency
                error('To obtain frequency-data (counts per second), specify BINSIZESEC.');
            end
            %% validate that input data match DATATYPE
            if areDistribData
                assert(isvector(data) || all(cellfun(@(c)isvector(c) && isnumeric(c),data))...
                    ,['For distribution-data, input variable must be either a cell' ...
                    ,' array of samples or a single vector of samples.']);
            else, assert(isvector(data) || size(data,2) <= 2 ...
                    ,['For histogram-data, input variable must be either a vector of' ...
                    ,' bin-values (counts or frequencies) or a N-by-2 matrix  where N' ...
                    ,' is the number of bins, the first column contains the bin-values' ...
                    ,' and the second column contains the standard deviation of the' ...
                    ,' respective bin.']);
            end
            
            %% use this function handle to transform axial data on demand
            % axial dimension, input to CIRCSTAT functions
            if areAxialData, axialDim = 2; else, axialDim = 1; end
            axTrans = @(x)circ_axial(x,axialDim);
            
            %% operations on input data based on dataType
            if areDistribData
                if isnumeric(data) % if it is only a vector, pack it into a cell
                    validateattributes(data,{'numeric'},{'vector'});
                    data = {data(:)};
                end
                % wrap angles into [0,360[, necessary for correct binning
                data = cellfun(@(c)mod(c,360),data,'UniformOutput',false);
                self.data = data;
                nSamples = numel(data);
                
                % edge calculation based on input
                if isempty(edges)
                    assert(hasHistcounts,['HISTCOUNTS not available, probably because ' ...
                        ,'you are running a MATLAB version below R2014b. In this case, the ' ...
                        ,'second input argument EDGES must be specified.']);
                    % use sample with highest number of data points for auto-binning
                    idxMostData = find(max(cellfun(@numel,data)));
                    [~,edges] = histcounts(data{idxMostData}); %#ok<FNDSB>
                end
                
                % deduce bin data from edges
                binSizeDeg = abs(edges(1) - edges(2));
                binCentersDeg = edges(1:end-1) + binSizeDeg/2;
                binCentersRad = deg2rad(binCentersDeg');
                
                binnedData = nan(numel(edges)-1,nSamples);
                for iS = 1:nSamples % calculate histogram counts
                    if hasHistcounts, [counts,~]  = histcounts(data{iS},edges);
                    else,             [counts,~]  = histc(data{iS},edges);
                                      counts(end) = [];
                    end
                    binnedData(:,iS) = counts;
                end
                
                if isFrequency % convert to counts per second
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
                
            else % already-binned data
                nSamples = NaN; % not feasible
                self.data = data;
                if isvector(data) % if it is a vector, use zeros for standard deviation
                    histData = data(:); % columnize
                    histData(:,2) = 0;
                else
                    histData = data;
                end
                histCnts = histData(:,1);
                
                % edge calculations based on input
                % for already-binned data: number of bins equals number of data points
                if isempty(edges),  edges = 0 : (360/numel(histCnts)) : 360;
                end
                
                % deduce bin data from edges
                binSizeDeg = abs(edges(1) - edges(2));
                binSizeRad = deg2rad(binSizeDeg);
                binCentersDeg = edges(1:end-1) + binSizeDeg/2;
                binCentersRad = deg2rad(binCentersDeg');
                
                if isempty(avgAng)
                    self.avgAng  = ...
                        rad2deg(circ_mean(axTrans(binCentersRad),histCnts)) / axialDim;
                    avgAng       = self.avgAng;
                end
                if isempty(r)
                    self.r = circ_r(axTrans(binCentersRad),histCnts,binSizeRad);
                    r      = self.r;
                end
                if isempty(avgAngCi)
                    self.avgAngCi = rad2deg(circ_confmean( ...
                        axTrans(binCentersRad),0.05,histCnts,binSizeRad)) / axialDim;
                    avgAngCi  = self.avgAngCi;
                end
                
                [self.rayleighP,self.rayleighZ] = circ_rtest( ...
                    axTrans(binCentersRad),histCnts,binSizeRad);
                pRayl = self.rayleighP;
                zRayl = self.rayleighZ;
                
                if isFrequency % convert to counts per second
                    histData = histData / binSizeSec; end
            end
            self.edges = edges;
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
            
            if ~isempty(ax) && isvalid(ax)
                figH = ax.Parent;
                polarAxs = ax;
            else
                figH = gcf;
                polarplot(0);
                polarAxs = gca;
            end
            hold(polarAxs,'on');
            
            figVisibleState = figH.Visible;
            figH.Visible = 'off';
            self.figH = figH;
%             set(0,'currentfigure',figH);
%             polarplot(0);hold on
            set(figH,'color',[1,1,1]) % white background
%             self.polarAxs = gca;
            self.polarAxs = polarAxs;
            
            % self-reference in property for hyper-redundancy (this is actually quite
            % handy if you want to retrieve the CircHist-object from a figure)
            polarAxs.UserData.circHistObj = self;
            
            polarAxs.ThetaAxisUnits = 'degrees';
            
            if exist('thetatickformat','file')
                  thetatickformat(polarAxs,'degrees'); % not supported below MATLAB R2016b
            else, polarAxs.ThetaAxis.TickLabelFormat = '%g\x00B0';
            end
            
            polarAxs.ThetaZeroLocation = 'top';
            polarAxs.Tag = 'Polar';
            
            polarAxs.Units = 'normalized';
            % decrease axes-dimensions to provide room for the scale bar
            polarAxs.Position(3) = polarAxs.Position(3) * 0.95;
            polarAxs.Position(4) = polarAxs.Position(4) * 0.95;
            self.scaleBarSide = scaleBarSide;
                        
            lineSp = '-'; % continuous lines
            self.barWidth = adjustSlope * binSizeDeg + 4; % bar width
            self.stdWidth = self.barWidth / 3;
            lineWAvgAng = adjustSlope * binSizeDeg/10 + 2;
            lineWR = lineWAvgAng * 1.5;
            
            %% draw bars and whiskers
            self.drawBars
            %%
            circR = max(rlim(polarAxs)); %radius of plot in data units
            
            % calculate baseline shift from center; depends on bin size for aesthetic reasons
            if ~isnan(baseLineOffset)
                % shift baseline by the fraction of the full radius and this value
                baseLineOffsetDiv = adjustSlope * binSizeDeg + baseLineOffset;
                baseLineOffset = -round(circR/baseLineOffsetDiv); % baseline shift in data units
            else, baseLineOffset = 0;
            end
            
            self.drawBaseLine;
            
            % shift radius baseline
            polarAxs.RLim = [baseLineOffset,circR];
            
            %% average angle and r
            % based on AREAXIALDATA, the average-angle and r lines are axes in the
            % histogram, or lines with the average angle as direction
            rL = rlim(polarAxs);
            avgAngRad = deg2rad(avgAng);
            if ~isempty(avgAng) && drawAvgAng %plot average angle
                if areAxialData
                    thetaAvgAng = [avgAngRad,avgAngRad+pi];
                    rhoAvgAng = [circR,circR];
                else
                    thetaAvgAng = [avgAngRad,avgAngRad];
                    rhoAvgAng = rL;
                end
                self.avgAngH = polarplot(self.polarAxs,thetaAvgAng,rhoAvgAng,lineSp...
                    ,'lineWidth',lineWAvgAng,'color',colorAvgAng,'Tag','avgAng');
            end
            if ~isempty(r) && ~isempty(avgAng) && drawR
                % make vector length relative to plot radius (after shift)
                rNorm = r * range(rL) + baseLineOffset;
                if areAxialData
                    thetaR = [avgAngRad,avgAngRad+pi];
                    rohR = [rNorm,rNorm];
                else
                    thetaR = [avgAngRad,avgAngRad];
                    rohR = [min(rL),rNorm];
                end
                self.rH = polarplot(self.polarAxs,thetaR,rohR,lineSp...
                    ,'lineWidth',lineWR,'color',colorR,'Tag','r');
            end
            
            % distance between confidence-interval line and plot border in percent
            avgAngCiWhiskOffset = 0.02;
            self.polarAxs.UserData.avgAngCiWhiskOffset = avgAngCiWhiskOffset;
            avgAngCiLW = 2; % line width

            avgAngCiRad = deg2rad(avgAngCi);
            % ISNAN(AVGANGCI) == TRUE if the requirements for confidence levels are not
            % met, see CIRC_CONFMEAN line 70 -> do not plot
            if ~isempty(avgAngCi) && ~isnan(avgAngCi) && drawAvgAngCi
                self.avgAngCiH = gobjects(0);
                avgAngCiPlotArgs = {'lineStyle',lineSp,'color',colorAvgAngCi ...
                    ,'Clipping','off'};
                
                whiskWidthEnd = self.stdWidth * 0.7; % width of whisker-endings
                thetaStepN = ceil(avgAngCi / 0.2); % plot line in 0.2-deg-steps
                thetaAvgAngCi = linspace( ...
                    avgAngRad - avgAngCiRad,avgAngRad + avgAngCiRad,thetaStepN);
                
                rhoAvgAngCi = max(rL) + avgAngCiWhiskOffset * range(rL);
                avgAngCiWhiskLen = avgAngCiWhiskOffset*2/3 * range(rL);
                rhoAvgAngCiWhiskEnd = ...
                    [rhoAvgAngCi + avgAngCiWhiskLen, rhoAvgAngCi - avgAngCiWhiskLen];
                rhoAvgAngCi = repmat(rhoAvgAngCi,thetaStepN,1);
                
                self.avgAngCiH(1,1) = polarplot(polarAxs,thetaAvgAngCi,rhoAvgAngCi,avgAngCiPlotArgs{:} ...
                    ,'lineWidth',avgAngCiLW,'Tag','avgAngCiWhisk');
                self.avgAngCiH(1,2) = polarplot(polarAxs,[thetaAvgAngCi(1),thetaAvgAngCi(1)] ...
                    ,rhoAvgAngCiWhiskEnd,avgAngCiPlotArgs{:} ...
                    ,'lineWidth',whiskWidthEnd,'Tag','avgAngCiWhiskEnd');
                self.avgAngCiH(1,3) = polarplot(polarAxs,[thetaAvgAngCi(end),thetaAvgAngCi(end)] ...
                    ,rhoAvgAngCiWhiskEnd,avgAngCiPlotArgs{:} ...
                    ,'lineWidth',whiskWidthEnd,'Tag','avgAngCiWhiskEnd');
                
                if areAxialData % plot again with mirrored angles
                    thetaAvgAngCi = thetaAvgAngCi + pi;
                    self.avgAngCiH(2,1) = polarplot(polarAxs,thetaAvgAngCi,rhoAvgAngCi ...
                        ,'lineWidth',avgAngCiLW,avgAngCiPlotArgs{:},'Tag','avgAngCiWhisk');
                    self.avgAngCiH(2,2) = polarplot(polarAxs,[thetaAvgAngCi(1),thetaAvgAngCi(1)] ...
                        ,rhoAvgAngCiWhiskEnd,avgAngCiPlotArgs{:} ...
                        ,'lineWidth',whiskWidthEnd,'Tag','avgAngCiWhiskEnd');
                    self.avgAngCiH(2,3) = polarplot(polarAxs,[thetaAvgAngCi(end),thetaAvgAngCi(end)] ...
                        ,rhoAvgAngCiWhiskEnd,avgAngCiPlotArgs{:} ...
                        ,'lineWidth',whiskWidthEnd,'Tag','avgAngCiWhiskEnd');
                end
            end
            
            %% edit axes
            polarAxs.Color = 'none'; % no background
                        
            polarAxs.ThetaAxis.FontSize = fontSize;
            polarAxs.LineWidth = 1;
            polarAxs.GridColor = 'k';
            polarAxs.GridAlpha = 0.5;
            
            % nifty way to dynamically create minor grid lines in 10 deg spacing, skipping
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
            details = sprintf(['N = %u , avgAng = %.2f\\circ\\pm%.2f\\circ, r = %.4f\n' ...
                'p_{Rayl} = %.3f, Z_{Rayl} = %.4f, p_{C. an.} = %.3f, R^{2}_{C. an.} = %.3f']...
                ,nSamples,avgAng,avgAngCi,r,pRayl,zRayl,corrAnP,corrAnR*corrAnR);
            title(polarAxs,details,'FontSize',9);
            %%
            colormap(polarAxs,white); % for axis appearance
            polarAxs.RTickLabel = [];
            
            if isempty(ax) % adjust figure-window size
                figH.OuterPosition(1) = figH.OuterPosition(1)/2; 
                figH.OuterPosition(2) = figH.OuterPosition(2)/2;
                figH.OuterPosition(3) = figH.OuterPosition(3)*1.2;
                figH.OuterPosition(4) = figH.OuterPosition(4)*1.5;
            end
            
            self.drawScale;
            figH.SizeChangedFcn = @self.redrawScale;
            % automatically DELETE object if the polaraxes is deleted
            polarAxs.DeleteFcn = @(~,~) delete(self);
            
            % this preserves the visiblity-state if the histogram is plotted into an
            % existing axes in an invisible figure.
            figH.Visible = figVisibleState;
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
            for iB = 1:nBars
                currAng = deg2rad(binCentersDeg(iB)); %angle of bar in rad
                currVal = self.histData(iB,1);
                if currVal == 0 || currVal <= lowerLim % skip bar-drawing
                    continue, end
                currStd = self.histData(iB,2);
                
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
            clr = self.colorBar; % baseline color, same as bar-color
            
            binSizeDeg = abs(self.edges(1) - self.edges(2));
            thetaBase = 0:deg2rad(binSizeDeg):2*pi; % degree-values between bars
            rhoBase = zeros(numel(thetaBase),1);
            polarplot(self.polarAxs,thetaBase,rhoBase,lineSp,'lineWidth',lineW...
                ,'color',clr,'Tag',tag);
            uistack(findobj(self.polarAxs,'Tag',tag),'bottom');
            uistack(findobj(self.polarAxs,'Tag',tag),'up');
        end
        %% drawScale; draw scale bar
        function drawScale(self)
            %drawScale Draws the scale bar. Used as SizeChangedFcn for figure so it is
            % drawn each time the figure size is changed. The scale bar is actually a
            % colorbar-object, thus it does not behave as neat as a conventional axis.
            
            scl = self.scaleBar;
            % non-empty and non-valid if the scale bar has been deleted
            if ~isempty(scl) && ~isvalid(scl), return; end
            
            figVisibleState = self.figH.Visible;
            self.figH.Visible = 'off'; % as recommended for SizeChangedFcn operations
            pAx = self.polarAxs;
                        
            sclLeft = strcmp(self.scaleBarSide,'left');

            % default font size. This line should be completely useless, but for some
            % reason, there are strange effects without it. The fact that the bar-drawing
            % is executed via SIZECHANGEDFCN seems to interrupt the link between the font
            % properties of the bar and the font properties of the parent axes, messing
            % something up.
            fontsz = 13; %#ok

            initialDraw = isempty(scl);
            if initialDraw % if the bar is drawn for the first time, use predefined label and font size
                hasLabel = true;
                label = self.axisLabel;
                fontsz = self.fontSize;
            else % read out (assumingly predefined) properties and use these for the next drawing
                hasLabel = ~isempty(scl.Label);
                if hasLabel, label  = scl.Label.String;
                             fontsz = scl.Label.FontSize;
                else,        label  = '';
                             fontsz = self.fontSize;
                end
                oldLineWidth = scl.LineWidth;
            end
            
            delete(scl); % deletes the previous bar
            
            scl = colorbar(pAx,'Location','east');
            self.scaleBar = scl; % re-assign
            
            % with this link, the label font-name is changed with the corresponding
            % axes-property using set(gca,'FontName',fontName), which is the default
            % behavior of regular colorbar labels.
            pAx.UserData.fontNameLink = linkprop([pAx,scl],'FontName');
            
            sclUnitsOld = scl.Units;
            scl.Units = 'pixels';
            polarAxsUnitsOld = pAx.Units;
            pAx.Units = 'pixels';
            
            sclWidth = 0.015;
            
            polarPos = pAx.Position; % position property == [left,bottom,width,height]
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
            % scale-bar must be located at the border of the labels. A side margin == 0
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
            sclTicks = pAx.RTick; % use same ticks as present in the POLARAXES
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
            % restore initial state of visibility
            self.figH.Visible = figVisibleState;
        end
        %% drawWhiteDisk
        function drawWhiteDisk(self)
            % draw white disk in the middle to obscure the axis-lines
            if ~isempty(self.whiteDiskH), delete(self.whiteDiskH); end
            
            offset = self.polarAxs.RLim(1);
            if offset >= 0, return, end % nothing to obscure
            
            pAx = self.polarAxs;
            rhoWhite = [offset,0];
            thetaLims = thetalim(pAx);
            if strcmp(pAx.ThetaAxisUnits,'radians'), thetaLims = rad2deg(thetaLims); end
            % 2 degree steps between white lines
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
            rlim(self.polarAxs,limits); %change limits
            
            % update line data
            self.drawBars
            if ~isempty(self.avgAngH) && isvalid(self.avgAngH)
                if self.areAxialData, self.avgAngH.RData(:) = limits(2);
                else,                 self.avgAngH.RData = limits;       end
            end
            if ~isempty(self.rH) && isvalid(self.rH)
                rNorm = self.r * range(limits) + limits(1);
                if self.areAxialData, self.rH.RData(:) = rNorm;
                else,                 self.rH.RData = [limits(1),rNorm]; end
            end
            if ~isempty(self.avgAngCiH) && all(all(isvalid(self.avgAngCiH)))
                oldRho = self.avgAngCiH(1,1).RData(1); % update CI circle-segment
                newRho = limits(2) + self.polarAxs.UserData.avgAngCiWhiskOffset * range(limits);
                avgAngCiRhoOffset = newRho - oldRho;
                whiskH = findobj(self.polarAxs,'Tag','avgAngCiWhisk');
                arrayfun(@(h)set(h,'RData',h.RData + avgAngCiRhoOffset),whiskH);
                
                avgAngCiWhiskLen = ... % update whisker-endings
                    self.polarAxs.UserData.avgAngCiWhiskOffset*2/3 * range(limits);
                newRhoAvgAngCiWhiskEnd = ...
                    [newRho + avgAngCiWhiskLen, newRho - avgAngCiWhiskLen];
                set(findobj(self.polarAxs,'Tag','avgAngCiWhiskEnd') ...
                    ,'RData',newRhoAvgAngCiWhiskEnd); % RDATA common for all whisker-ends
            end
            
            if isvalid(self.thetaLabel) && ~isempty(self.thetaLabel) ...
                    && ~isempty(self.thetaLabel.String)
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
            
            self.drawBaseLine;
            self.drawWhiteDisk; % adjust background
        end
        %% thetalim, enables calling obj.thetalim(limits)
        function limitsOut = thetalim(self,limits)
            %thetalim Returns the limits of the theta-axis if called without input
            %arguments, else sets the limits.
            %
            %   limits = obj.thetalim;
            %   obj.thetalim(newLimits); % where NEWLIMITS == [lower,upper]
            if nargin < 2, limitsOut = thetalim(self.polarAxs);
            else,          thetalim(self.polarAxs,limits);
            end
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
            self.thetaLabel = text(self.polarAxs,deg2rad(th),rlims(2) + range(rlims)*0.2 ...
                ,txt,'HorizontalAlignment','center','VerticalAlignment','cap' ...
                ,'Rotation',txtRot,'FontSize',self.fontSize);
            self.thetaLabel.UserData.location = location;
            % make label-color be linked to theta-axis color (default behavior for
            % labels). For some reason, the link only works when THETALABEL.COLOR is
            % changed, not the other way around ...
            self.polarAxs.UserData.thetaLabelColorLink = ...
                linkprop([self.polarAxs.ThetaAxis,self.thetaLabel],'Color');
        end
        %% change scale bar side
        function set.scaleBarSide(self,side)
            %scaleBarSide Determines the side of the scale bar.
            %
            %   obj.scaleBarSide(side); where SIDE is either 'left' or 'right'
            assert(any(strcmpi(side,{'left','right'})) ...
                ,'Property SCALEBARSIDE must be ''left'' or ''right''.');
            if strcmpi(side,self.scaleBarSide), return; end % nothing to do
            
            pAx = self.polarAxs;
            isLeft = strcmpi(side,'left');
            % initially, the plot is not positioned in the center of the axes, but
            % slightly right, and there seems to be no way to center it, so the first
            % adjustment needs to be different to subsequent adjustments
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
        %% change font size
        function set.fontSize(self,newFontSize)
            self.fontSize = newFontSize;
            if ~isempty(self.polarAxs) && isvalid(self.polarAxs)
                self.polarAxs.ThetaAxis.FontSize = newFontSize;
            end
        end
        %% change bar color
        function set.colorBar(self,color)
            self.colorBar = color;
            lineObjArr = findobj(self.polarAxs,'Tag','histBar','-or','Tag','baseLine');
            set(lineObjArr,'color',color);
        end
        %% change r-line color
        function set.colorR(self,color)
            self.colorR = color;
            set(self.rH,'Color',color);
        end
        %% change average-angle line color
        function set.colorAvgAng(self,color)
            self.colorAvgAng = color;
            set(self.avgAngH,'Color',color);
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
            %   circHistObj.stdWidth(newWidth);
            
            oldWidth = self.stdWidth;
            self.stdWidth = width;
            
            % happens at object construction
            if isempty(oldWidth), return; end
            
            scalingFactor = width / oldWidth;
            
            % get handles of different line objects
            maskLines = strcmp({self.stdH.Tag},'stdWhisk');
            stdMain = self.stdH(maskLines);
            stdEnding = self.stdH(~maskLines);
            
            oldEndingWidth = stdEnding(1).LineWidth; % scale proportionally
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
            % (default = 90 dpi). Specify RESOL as a string of the pattern '-r90'.
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
            % optionally specified by additional parameters which must be Name-Value pairs
            % as accepted by POLARPLOT. Optionally returns the graphics-object handle.
            %
            %   obj.drawCirc(rho);
            %   obj.drawCirc(rho,Name,Value);
            %   h = obj.drawCirc(___);
            %
            % See also polarplot
            
            assert(isnumeric(rho) && isscalar(rho) ...
                ,'First input RHO must be a numeric scalar.');
            
            theta = linspace(0,2*pi,100);
            h = polarplot(self.polarAxs,theta,repmat(rho,size(theta)),varargin{:});
            
            if isvalid(self.avgAngH), uistack(self.avgAngH,'top'); end
            if isvalid(self.rH),      uistack(self.rH,'top'); end
            
            if nargout > 0 
                hOut = h; end
        end
        %% drawArrow
        function hOut = drawArrow(self,theta,rho,varargin)
            %drawArrow Draws arrows (ANNOTATION objects) from the center of the plot at
            % the specified angles (in degrees) with the specified lengths (in data
            % units). To change the arrow properties, pass additional Name-Value pairs as
            % accepted by annotation('arrow',_,Name,Value). To change properties after
            % drawing, use obj.arrowH(N).PROPERTY, where N is the index of the arrow
            % object and PROPERTY is the respective property; alternatively, call this
            % function with an output argument and use this variable to access the
            % properties. Note that the arrows are ANNOTATION objects and that their
            % positional information is anchored in the figure coordinates, not the axes
            % coordinates, so the position might get messed up if the figure's or axes'
            % size/orientation/... is changed. In this case, call OBJ.UPDATEARROWS which
            % should re-position all arrows correctly. Tip: If an arrow is supposed to end
            % at a specific point on the scale, specify RHO as R+abs(min(rlim)), where R
            % is the desired point on the scale.
            %
            %   Input:
            %       theta       Scalar or vector of angles in degree.
            %       rho         Either a scalar or a vector of the same size as THETA. If
            %                   omitted or specified empty, the arrows end at the edge of
            %                   the plot. If THETA is a non-scalar vector and RHO is a
            %                   scalar, all arrows have the same lenght RHO. Note that
            %                   drawArrow(theta,[]) produces the same arrows as
            %                   drawArrow(theta,range(rlim)), but that, in the former
            %                   case, the arrows stick to the plot edge if the scale
            %                   limits are changed, whereas in the latter, the arrows keep
            %                   their lengths in data units.
            %
            %   obj.drawArrow(30); % arrow at 30 deg from the center to the edge of the plot
            %   obj.drawArrow(30,10); % same as above but the arrow is 10 data units long
            %   obj.drawArrow(30,10,Name,Value); % additional parameters as accepted by
            %                                      ANNOTATION
            %   obj.drawArrow(30,[],Name,Value); % arrow ends at edge of the plot
            %   h = obj.drawArrow(_); % returns the handle of the arrow objects
            %
            % See also annotation CircHist.updateArrows
                        
            validateattributes(theta,{'numeric'},{'vector'});
            theta = theta(:);
            nAr = numel(theta);
            
            drawnow
            pax = self.polarAxs;
            axUnitsOld = pax.Units; % save initial units to preserve them
            pax.Units = 'pixels';
            
            % empty means that the arrows end at the plot edge (default behavior)
            if nargin < 3, rho = []; end 
            
            if    isscalar(rho) || isempty(rho), rho = repmat(rho,size(theta));
            else, assert(numel(rho) == nAr,['If RHO is not a scalar, it must' ...
                    ,' have the same number of elements as THETA.']);
            end
            
            % get center coordinates and radius of polar plot
            cenX = pax.Position(1) + pax.Position(3)/2;
            cenY = pax.Position(2) + pax.Position(4)/2;
            radiusFull = min([pax.Position(3),pax.Position(4)])/2;
            rlimRange = range(rlim);
            
            switch pax.ThetaZeroLocation % angle-offset depending on THETAZEROLOCATION
                case 'right' , thetaOffset = 0;
                case 'top'   , thetaOffset = 90;
                case 'left'  , thetaOffset = 180;
                case 'bottom', thetaOffset = 270;
            end
            
            arrwH = gobjects(nAr,1); % array to store all created arrows
            
            for iAr = 1:nAr
                cTheta = theta(iAr);
                if isempty(rho), cRho = rlimRange;
                else,            cRho = rho(iAr);
                end
                
                % initialize zero-length arrow because else the UNIT change will not work
                cArrwH = annotation(self.figH,'arrow',[0,0],[0,0],varargin{:});
                arrwH(iAr) = cArrwH;
                
                arrwUnitsOld = cArrwH.Units;
                cArrwH.Units = 'pixels';
                
                cArrwH.UserData.rho = cRho; % save length and angle for use in UPDATEARROWS
                cArrwH.UserData.theta = cTheta;
                
                % calculate arrow length in pixels
                arrwLen = cRho * radiusFull / rlimRange;
                
                % calculate pixel coordinates of arrow-end
                [x,y] = pol2cart(deg2rad(cTheta + thetaOffset),arrwLen);
                cArrwH.X = [cenX,cenX + x]; % adjust arrow coordinates
                cArrwH.Y = [cenY,cenY + y];
                
                cArrwH.Units = arrwUnitsOld; % restore units
            end
            
            self.arrowH = [self.arrowH;arrwH]; % add to CircHist-object property
            
            self.polarAxs.Units = axUnitsOld; % restore units
            
            if nargout > 0, hOut = arrwH; end
        end
        %% update arrows
        function updateArrows(self)
            %updateArrows Updates all drawn arrows. Call this manually if they are messed
            %up.
            %
            %   obj.updateArrows;
            
            pax = self.polarAxs;
            % delete invalid arrows and RETURN if there are none
            mValid = isvalid(self.arrowH);
            self.arrowH = self.arrowH(mValid);
            if nnz(mValid) < 1, return, end
            
            arrwH = self.arrowH;
            axUnitsOld = pax.Units;
            pax.Units = 'pixels';
            arrwUnitsOld = arrwH(1).Units;
            set(arrwH,'Units','pixels');
            
            switch pax.ThetaZeroLocation % angle-offset depending on THETAZEROLOCATION
                case 'right' , thetaOffset = 0;
                case 'top'   , thetaOffset = 90;
                case 'left'  , thetaOffset = 180;
                case 'bottom', thetaOffset = 270;
            end
            
            % get center coordinates and radius of polar plot
            cenX = pax.Position(1) + pax.Position(3)/2;
            cenY = pax.Position(2) + pax.Position(4)/2;
            radiusFull = min([pax.Position(3),pax.Position(4)])/2;
            rlimRange = range(rlim);
            
            for iAr = 1:numel(arrwH)
                cArrwH = arrwH(iAr);
                cRho = cArrwH.UserData.rho;
                cTheta = cArrwH.UserData.theta;
                if isempty(cRho), cRho = rlimRange; end
                arrwLen = cRho * radiusFull / rlimRange;
                [x,y] = pol2cart(deg2rad(cTheta + thetaOffset),arrwLen);
                cArrwH.X = [cenX,cenX + x];
                cArrwH.Y = [cenY,cenY + y];
            end
            
            self.polarAxs.Units = axUnitsOld; % restore units
            set(arrwH,'Units',arrwUnitsOld);
        end
    end
    methods (Static)
        %% redrawScale
        function redrawScale(~,~)
            % Called when the figure window-size changes, should not be called manually
            fH = gcbf;
            set(0,'currentfigure',fH);
            % get all CIRCHIST-objects in this figure via the polaraxes' USERDATA entry
            pAxH = findobj(fH,'-function',@(o)isfield(o.UserData,'circHistObj'));
            if ~isempty(pAxH) % loop through objects, call DRAWSCALE and UPDATEARROWS
                arrayfun(@(ax) ax.UserData.circHistObj.drawScale,pAxH);
                arrayfun(@(ax) ax.UserData.circHistObj.updateArrows,pAxH);
            end
        end
        %% entry for separate usage-example script
        exampleCircHist
    end
end