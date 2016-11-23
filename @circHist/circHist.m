classdef CircHist < handle
    %CircHist   Class representing a figure with a circular histogram. Constructing an
    %object creates a new figure with polar coordinates containing a histogram. Circular
    %statistics (average angle, mean resultant vector, Rayleigh test of uniformity and
    %correlation analysis) are automatically calculated using the CircStat toolbox and
    %saved as object properties. Note that this is a handle class.
    %
    %   Requirements: CircStat toolbox (mathworks.com/matlabcentral/fileexchange/10676)
    %                 Ephys class
    %
    %   To change the axis limits, use CircHistObj.setRLim([lower,upper]).
    %   To change visual properties, either use the name-value pairs for the constructor
    %   as specified below or access the graphics-objects, namely CircHistObj.polarAxs for
    %   the coordinate system (font size, line width, etc.) and CircHistObj.scaleBar for
    %   the scale bar (change the axis label here). The scale bar is drawn anew each time
    %   the figure window changes size; note that this is not working flawlessly and check
    %   that the scale bar matches the coordinate-grid after changing the size.
    %   To adjust the bars, standard-deviation whiskers, phimax line and r line after
    %   plotting, use their handles which are saved as properties. Access them using dot
    %   notation, e.g., obj.phimaxH.lineWidth and delete them using delete(handle).
    %
    %   Usage:  CircHist(data,edges);
    %           CircHist(data,edges,Name,Value);
    %           obj = CircHist(___);
    %
    %   Methods:
    %       setRLim([lower,upper])      Change axis limits (get current limits by calling
    %                                   rlim)
    %       colorBar                    Change bar color (usage: obj.colorBar = newcolor)
    %       barWidth                    Change bar width (usage: obj.barWidth = newWidth)
    %       colorStd                    Change standard-deviation-whisker color
    %       stdWidth                    Change standard-deviation-whisker width
    %       toPdf(fileName)             Save as pdf-file with the specified file name
    %
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
    %   ---Optional input:
    %   dataType        'distribution'(default)/'psth'. Type of input data: 'distribution'
    %                   treats the input data as distributions of degree values; the data
    %                   of each input vector are binned inside the specified edges, each
    %                   bin is averaged and the standard deviation is calculated. 'psth'
    %                   treats the data as already binned data (converted to frequency),
    %                   e.g. from a peri-stimulus-time-histogram; if the input matrix has
    %                   a second column, it is taken as the standard-deviation values for
    %                   each bar.
    %
    %   psthType        'frequency'(default)/'count'. Based on this, the axis is labelled.
    %
    %   binSizeSec      Width of bins in seconds. Note that this needs to be specified if
    %                   dataType is 'distribution' and psthType is 'frequency', as well as
    %                   when dataType is 'psth' and psthType is 'count'.
    %
    %   includePhimax   'on'(default)/'off', plots average angle
    %
    %   phimax          Numeric value of the average angle. Should be specified if
    %                   dataType is 'psth' and includePhimax is 'on'.
    %
    %   includeR        'on'(default)/'off', plots resultant vector length as a black
    %                   overlay bar on top of the average angle. The black bar's length
    %                   equals the r-value in percent of the coordinate-system diameter.
    %
    %   r               Numeric value of the resultant vector length. Should be specified
    %                   if dataType is 'psth' and includeR is 'on'.
    %
    %   baseLineShift   Numeric value (default = 2) specifying the factor that scales the
    %                   size of the offset in the middle of the plot. The default value
    %                   produces nice results; lower the value (negative values allowed)
    %                   to increase the size, increase the value to decrease the size.
    %                   Specify it as NaN to produce an offset of 0. (This functionality
    %                   might need improvement, maybe just implement it as a percent value
    %                   of the plot diameter)
    %
    %   adjustSlope     Slope that defines how strong optical properties such as bar width
    %                   scale with bin-size; default = 0.3.
    %
    %   axialData       True(default)/false, specifying whether or not input data are
    %                   axial; else they are considered circular. This is taken into
    %                   account for statistical computations (axial data are are
    %                   multiplied by 2 before calculation).
    %
    %   ax              Axes handle to plot diagram into; becomes polarAxs property
    %                   (experimental feature, working in principle).
    %
    %   colorBar        Color of bars (default = [0 .45 .74]; (Matlab blue)).
    %   colorStd        Color of standard-deviation lines (default = 'k').
    %   colorPhimax     Color of phimax line (default = [.85 .33 .1]; (orange)).
    %   colorR          Color of r line (default = 'k').
    %   fontSize        Font size of axis labels (default = 13).
    %
    %
    %  ---Notes:
    %       Change axis label: obj.scaleBar.Label.String = 'your label';
    %
    %       Each standard-deviation whisker consists of a long line representing the
    %           magnitude of the deviation and a very short, thick line that marks the tip
    %           of the deviation-line. Both line types are comprised in the handle-array
    %           obj.stdH; however, they can be separately accessed by using
    %           findobj(obj.figH,'Tag',typeTag) where typeTag is either 'stdWhisk' for the
    %           "main" line or 'stdWhiskEnd' for the tips.
    %
    %
    %  ---Author: Frederick Zittrell
    %
    %
    properties (SetAccess = immutable)
        data            %Required input: Data.
        edges           %Required input: Histogram edges;
                
        dataType        %Optional input; 'distribution'(default)/'psth'
        psthType        %Optional input; 'frequency'(default)/'count'
        binSizeSec      %Optional input; Width of bins (s)
        includePhimax   %Optional input; 'on'(default)/'off'
        phimax          %Optional input; Numeric value of the average angle
        includeR        %Optional input; 'on'(default)/'off'
        r               %Optional input; Numeric value of the resultant vector length
        baseLineShift   %Optional input; scaling factor for the size of the offset in the middle of the diagram
        adjustSlope     %Optional input; Slope for scaling of visual properties with bin size
        axialData       %Optional input; True(default)/false
        
        polarAxs        %Polaraxes handle. Change visual properties such as line width of the axes here.
        figH            %Handle to figure where diagram is plotted
        phimaxH         %Handle to the phimax line
        rH              %Handle to the r line
        barH            %Array of handles to the bars (which are line objects)
        stdH            %Array of handles to the standard-deviation whiskers (which are line objects)
        
        psthData        %PSTH data as plotted; 1st column average counts, 2nd column standard deviations
        rayleighP       %P-value of Rayleigh test of uniformity
        rayleighZ       %Z-value of Rayleigh test of uniformity
        corrAnP         %P-value of correlation analysis
        corrAnR         %R-value of correlation analysis (square this to get the coefficient of determination)
    end
    
    properties
        scaleBar        %Handle of scale bar. Use to access visual properties.
        axisLabel       %Label of scale bar as originally set (change the scaleBar.Label.String property to adjust)
        
        colorBar        %Optional input; Color of bars (default = [0 .45 .74]; (Matlab blue))
        colorStd        %Optional input; Color of standard-deviation lines (default = 'k')
        colorPhimax     %Optional input; Color of phimax line (default = [.85 .33 .1]; (orange))
        colorR          %Optional input; Color of r line (default = 'k')
        fontSize        %Optional input; Font size of axis labels (default = 13)
        
        barWidth        %Width of bars. Change this property to adjust the bar width after plotting.
        stdWidth        %Width of standard-deviation whiskers. For adjustment after plotting.
    end
    
    methods
        %% constructor
        function self = CircHist(data,edges,varargin)   
            
            %% validate and parse input
            validateattributes(data,{'numeric','cell'},{'nonempty'});
            validateattributes(edges,{'numeric'},{'vector'});
            self.data  = data;
            self.edges = edges;
            
            validColor = @(x)validateattributes(x,{'numeric','char'},{'vector'});
            validScalarNum = @(N) isnumeric(N) & isscalar(N);
            
            %default values
            def.dataType      = 'distribution';
            def.psthType      = 'frequency';
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
            addOptional(pr,'psthType',  def.psthType,@(str) any(strcmpi(str,{'frequency','count'})));
            addOptional(pr,'binSizeSec',def.binSizeSec,...
                @(x) validateattributes(x,{'numeric'},{'scalar'}));
            
            addParameter(pr,'includePhimax',def.includePhimax);
            addParameter(pr,'phimax'       ,def.phimax,validScalarNum);
            addParameter(pr,'includeR'     ,def.includeR);
            addParameter(pr,'r'            ,def.r,validScalarNum);
            addParameter(pr,'baseLineShift',def.baseLineShift,validScalarNum);
            addParameter(pr,'adjustSlope'  ,def.adjustSlope,validScalarNum);
            addParameter(pr,'axialData'    ,def.axialData,...
                @(x) validateattributes(x,{'logical'},{'scalar'}));
            addParameter(pr,'ax'           ,def.axes...
                ,@(x) validateattributes(x,{'matlab.graphics.axis.Axes'},{'scalar'}));
            
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
            self.axialData     = pr.Results.axialData;
            axialData          = self.axialData;
            ax                 = pr.Results.ax;
            
            self.colorBar    = pr.Results.colorBar;
            colorBar         = self.colorBar;
            self.colorStd    = pr.Results.colorStd;
            colorStd         = self.colorStd;
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
                self.axisLabel = 'Activity /imp/s';
            else
                self.axisLabel = 'Spikes per bin';
            end
            % deduce histogram data from edges
            binSizeDeg = abs(edges(1) - edges(2));
            binCentersDeg = edges(1:end-1) + binSizeDeg/2;
            %% operations on input data based on dataType
            if areDistribData
                if isnumeric(data) %if it is only a vector, pack it into a cell
                    validateattributes(data,{'numeric'},{'vector'});
                    data = {data};
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
                self.phimax  = Ephys.phimax(degPool,axialData);
                phimax       = self.phimax;
                self.r       = Ephys.r(degPool,axialData);
                r            = self.r;
                [self.rayleighP,self.rayleighZ] = Ephys.rayl(degPool,axialData);
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
                Ephys.corrAn(binCentersDeg(:),psthData(:,1),axialData);
            corrAnR = self.corrAnR;
            corrAnP = self.corrAnP;
            
            %% initialize figure, set visual properties
            if ~isempty(ax) && isvalid(ax)
                currFig = gcf;
                currFig.CurrentAxes = ax;
            else
                currFig = figure;
            end
            currFig.Visible = 'off';
            self.figH = currFig;
            set(0,'currentfigure',currFig);
            polarplot(0);hold on
            self.polarAxs = gca;
            polarAxs = self.polarAxs;
            polarAxs.ThetaZeroLocation = 'top';
            polarAxs.Tag = 'Polar';
            
            lineSp = '-'; %continuous lines
            colorBase = colorBar; %baseline color
            lineWBar = adjustSlope * binSizeDeg + 4; %bar width
            self.barWidth = lineWBar;
            lineWStd = lineWBar / 3;
            self.stdWidth = lineWStd;
            %whisker-endings consist of a short line starting at (whisker-end - whiskLen)
            %and ending at the whisker-end
            lineWStdWhisk = lineWBar * 0.7; %width
            whiskLen = .15;
            lineWPhimax = adjustSlope * binSizeDeg/10 + 2;
            lineWR = lineWPhimax * 1.5;
            lineWBase = 1;
            
            %% plot bars
            nBars = numel(binCentersDeg);
            for n = 1:nBars
                currAng = deg2rad(binCentersDeg(n)); %angle of bar in rad
                currVal = psthData(n,1);
                currStd = psthData(n,2);
                
                %specify the bar: a line between the base line of the plot (at theta = angle and rho =
                %0) and the point specified by the angle and the radius-value (at theta = angle and
                %rho = bar height)
                thetaVal = [currAng,currAng];
                rhoVal = [0,currVal];
                polarplot(thetaVal,rhoVal,lineSp,'lineWidth',lineWBar,'color',colorBar...
                    ,'Tag','histBar');
                
                if currStd > 0
                    %specify the standard deviation whisker: a line between the outer end of the bar
                    %(at theta = angle and rho = bar height) and the same end extended by the standard
                    %deviation (at theta = angle and rho = bar height + standard deviation)
                    thetaStd = [currAng,currAng];
                    rhoStd = [currVal,currVal+currStd];
                    polarplot(thetaStd,rhoStd,lineSp,'lineWidth',lineWStd...
                        ,'color',colorStd,'Tag','stdWhisk');
                    
                    %whisker-endings: short lines thicker than the whisker, starting at (whisker-end -
                    %whiskLen) and ending at the whisker-end
                    thetaStdWhisk = thetaStd;
                    rhoStdWhisk = [currVal+currStd-whiskLen,currVal+currStd];
                    polarplot(thetaStdWhisk,rhoStdWhisk,lineSp,'lineWidth',...
                        lineWStdWhisk,'color',colorStd,'Tag','stdWhiskEnd');
                end
            end
            self.barH = findobj(self.figH,'Tag','histBar');
            self.stdH = findobj(self.figH,'Tag','stdWhisk','-or','Tag','stdWhiskEnd');
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
            
            todo ='';
            %TODO: plot white lines with the same width as the 'bars' between baseline and center to
            %make the center-hole white without ticks and lines. The below crap does not work yet.
            % rhoWhite = [baseLineOffset,0-lineWBase];
            % for n = 1:nBars
            %     currAng = deg2rad(binCentersDeg(n)); %angle of bar in rad
            %     polarplot([currAng,currAng],rhoWhite,'linestyle','-','color','w','linewidth',lineWBar);
            % end
            
            %draw baseline
            thetaBase = 0:deg2rad(binSizeDeg):2*pi; %degree-values between bars
            rhoBase = zeros(numel(thetaBase),1);
            polarplot(thetaBase,rhoBase,lineSp,'lineWidth',lineWBase,'color',colorBase...
                ,'Tag','baseLine');
            
            %shift radius baseline
            polarAxs.RLim = [baseLineOffset,circR];
            
            %% phimax and r
            phimaxRad = deg2rad(phimax);
            if ~isempty(phimax) && includePhimax %plot long diagonal phimax
                thetaPhimax = [phimaxRad,phimaxRad+pi];
                rhoPhimax = [circR,circR];
                self.phimaxH = polarplot(thetaPhimax,rhoPhimax,lineSp...
                    ,'lineWidth',lineWPhimax,'color',colorPhimax,'Tag','phimax');
            end
            if ~isempty(r) && ~isempty(phimax) && includeR
                rNorm = r * range(rlim) + baseLineOffset; %make vector length relative to plot radius (after shift)
                thetaR = [phimaxRad,phimaxRad+pi];
                rohR = [rNorm,rNorm];
                self.rH = polarplot(thetaR,rohR,lineSp,'lineWidth',lineWR,'color',colorR...
                    ,'Tag','r');
            end
            
            %% edit axes
            polarAxs.Color = 'none'; %no background
            labels = polarAxs.ThetaTickLabel;
            % add °-sign to labels
            for n = 1:numel(labels)
                labels{n} = [labels{n},'°'];end
            polarAxs.ThetaTickLabel=labels;
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
            
            %% title
            details = sprintf(['N = %u , phimax = %.2f°, r = %.4f\np_{Rayl} = %.3f, ',...
                'Z_{Rayl} = %.4f, p_{C. an.} = %.3f, R²_{C. an.} = %.3f'],...
                nSamples,phimax,r,pRayl,zRayl,corrAnP,corrAnR*corrAnR);
            title(details,'FontSize',9)
            
            %%
            colormap white %for axis appearance
            polarAxs.RTickLabel = [];
            
            if isempty(ax)
                currFig.OuterPosition(1) = currFig.OuterPosition(1)/2; %figure-window size
                currFig.OuterPosition(2) = currFig.OuterPosition(2)/2;
                currFig.OuterPosition(3) = currFig.OuterPosition(3)*1.2;
                currFig.OuterPosition(4) = currFig.OuterPosition(4)*1.5;
            end
            
            self.drawScale;
            currFig.Visible = 'on';
            currFig.SizeChangedFcn = {@self.redrawScale,self};
        end
        %% drawScale; draw scale bar
        function drawScale(self)
            %drawScale Draws the scale bar. Used as SizeChangedFcn for figure so it is
            %drawn each time the figure size is changed. The scale bar is actually a
            %colorbar-object, thus it does not behave as neat as a conventional axis.
            
            pAx = self.polarAxs;
            scl = self.scaleBar;
            
            %default font size. This line should be completely useless, but for some
            %reason, there are strange effects without it.
            fontsz = 13;

            initialDraw = isempty(scl);
            if initialDraw % if the bar is drawn for the first time, use predefined label and font size
                label = self.axisLabel;
                fontsz = self.fontSize;
            else % read out (assumingly predefined properties) and use these for the next drawing
                label = scl.Label.String;
                fontsz = scl.Label.FontSize;
                oldLineWidth = scl.LineWidth;
            end
            
            delete(scl); %deletes the previous bar
            
            scl = colorbar('Location','westoutside');
            % plotted like this, the scale height equals the polar plot diameter
            self.scaleBar = scl; % re-assign
            
            sclUnitsOld = scl.Units;
            scl.Units = 'pixels';
            polarAxsUnitsOld = pAx.Units;
            pAx.Units = 'pixels';
            
            scl.Position(1) = pAx.Position(1);
            polarPos = pAx.Position; %position property = [left,bottom,width,height]
            polarBot = polarPos(2);
            polarHeight = polarPos(4);
            
            sclWidth = 0.015;
            sclHeight = polarHeight/2; %halve the height -> bar spans from center to top
            offset = sclHeight * abs(pAx.RLim(1))/range(pAx.RLim);
            sclHeight = sclHeight - offset; %adjust baseline shift
            sclBot = polarBot + polarHeight - sclHeight; %calculate baseline position
            scl.Position = [pAx.Position(1),sclBot,sclWidth,sclHeight];
            %move scale left if it overlays the polar-axis label (hard-coded, empirical values, sadly)
            if offset < 15
                moveLeft = 15;
            else
                moveLeft = 0;
            end
            scl.Position(1) = pAx.Position(1) - moveLeft; %move
            
            scl.Label.String = label;
            scl.Limits = [0,max(pAx.RLim)];
            sclTicks = pAx.RTick; %use ticks as produced by the polarplot function
            sclTicks = sclTicks(sclTicks >= 0);
            scl.Ticks = sclTicks;
            scl.Box = 'off';
            scl.TickLength = 0.04;
            scl.FontSize = fontsz;
            
            scl.Units = sclUnitsOld;
            pAx.Units = polarAxsUnitsOld;
            
            if ~initialDraw
                scl.LineWidth = oldLineWidth;
            end
        end
        %% change scale limits
        function setRLim(self,limits)
            %setRLim Change scale limits to new two-element vector with [lower,upper]. Get
            %the current limits by calling rlim.
            rlim(limits); %change limits
            
            if limits(1) > 0
                todo = ''; % scale-bar bug in this case---------------------------------------------------
                rDataOffset = limits(1);
            else
                rDataOffset = 0;
            end
            
            nBars = numel(self.barH);
            for n = 1:nBars
                currBar = self.barH(n);
                barHeight = currBar.RData(2);
                currBar.RData = [rDataOffset,barHeight];                
            end
            
            self.phimaxH.RData(:) = limits(2); %update line data            
            self.rH.RData(:) = self.r * range(limits) + limits(1);
            
            % Sometimes, a warning with this identifier is issued for obscure reasons. To
            % suppress this warning, it is converted to an error by calling
            % warning('error',_). This error can then be caught and ignored.
            wrn = warning('error','MATLAB:callback:error'); %#ok<CTPCT>
            try   self.drawScale; %update scale
            catch ME
                if ~strcmp(ME.identifier,wrn.identifier)
                    rethrow(ME);end % re-throw error if it is not this error
            end
            warning(wrn); % set error back to being a warning
        end
        %% change bar color
        function set.colorBar(self,color)
            self.colorBar = color;
            lineObjArr = findobj(self.figH,'Tag','histBar','-or','Tag','baseLine');
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
            self.colorStd = color;
            set(self.stdH,'color',color);
        end
        %% change whisker width
        function set.stdWidth(self,width)
            % Since each whisker consists of two line-objects with different widths (the
            % "main" line and the ending), the width of the ending is scaled
            % proportionally to the width-change of the "main" line.
            
            oldWidth = self.stdWidth;
            self.stdWidth = width;
            
            %happens at object construction; this is NOT elegant and I am sorry
            if isempty(oldWidth) 
                return; end
            
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
            %toPdf  Save histogram as (fileName).pdf.
            
            % allow robust input with and without pdf-file-extension
            if strcmpi(fileName(end-3:end),'.pdf')
                fileName(end-3:end) = ''; end
            
            saveas(self.figH,[fileName,'.pdf']);
        end
    end
    methods (Static)
        %%
        function redrawScale(~,~,circHistObj)
           circHistObj.drawScale;
       end
    end
end