classdef SaccadeLabeling_v2 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        EyeTraces                     matlab.ui.control.UIAxes
        Velocities                    matlab.ui.control.UIAxes
        RadialView                    matlab.ui.control.UIAxes
        ControlsPanel                 matlab.ui.container.Panel
        NextSaccadeButton             matlab.ui.control.Button
        PreviousSaccadeButton         matlab.ui.control.Button
        OnsetLeft                     matlab.ui.control.Button
        OnsetRight                    matlab.ui.control.Button
        SaccadeOnsetLabel             matlab.ui.control.Label
        OffsetRight                   matlab.ui.control.Button
        OffsetLeft                    matlab.ui.control.Button
        SaccadeOffsetLabel            matlab.ui.control.Label
        DeleteCurrentSaccadeButton    matlab.ui.control.Button
        AddPSOtocurrentSaccadeButton  matlab.ui.control.Button
        PSOOffsetLeft                 matlab.ui.control.Button
        PSOOffsetRight                matlab.ui.control.Button
        PSOOffsetLabel                matlab.ui.control.Label
        DeleteCurrentPSOButton        matlab.ui.control.Button
        Label                         matlab.ui.control.Label
        BinSacs                       matlab.ui.control.UIAxes
        ZoomVelocity                  matlab.ui.control.UIAxes
        ZoomTrace                     matlab.ui.control.UIAxes
        CreatenewSaccadesPanel        matlab.ui.container.Panel
        NewSaccadeButton              matlab.ui.control.Button
        OnsetEditFieldLabel           matlab.ui.control.Label
        OnsetEditField                matlab.ui.control.NumericEditField
        OffsetEditFieldLabel          matlab.ui.control.Label
        OffsetEditField               matlab.ui.control.NumericEditField
        MakeBrushedSaccadesButton     matlab.ui.control.Button
        ActivateBrushButton           matlab.ui.control.Button
        MakeBrushedPSOButton          matlab.ui.control.Button
        TrialNavigationPanel          matlab.ui.container.Panel
        PreviousTrialButton           matlab.ui.control.Button
        NextTrialButton               matlab.ui.control.Button
        SaveCurrentTrialButton        matlab.ui.control.Button
        FilterBankPanel               matlab.ui.container.Panel
        MovingMeanSwitchLabel         matlab.ui.control.Label
        MovingMeanSwitch              matlab.ui.control.Switch
        MovingMedianSwitchLabel       matlab.ui.control.Label
        MovingMedianSwitch            matlab.ui.control.Switch
        VelocityLabel                 matlab.ui.control.Label
        VelocitySwitch                matlab.ui.control.RockerSwitch
        SmoothingKernelSizeKnobLabel  matlab.ui.control.Label
        SmoothingKernelSizeKnob       matlab.ui.control.Knob
        KernelSizeLabel               matlab.ui.control.Label
        Label2                        matlab.ui.control.Label
        Acceleration                  matlab.ui.control.UIAxes
        SaccadeDataLabel              matlab.ui.control.Label
        ManualSaccadeLabelingGUI      matlab.ui.control.Label
    end

    
    methods (Access = private)
        
        function [] = updatePlot(app,currentSac)
            
            cla(app.EyeTraces)
            cla(app.Velocities)
            cla(app.RadialView)
            cla(app.ZoomTrace)
            cla(app.ZoomVelocity)
            cla(app.BinSacs)
            cla(app.Acceleration)
            
            
            eyevec        = getappdata(0,'eyevec');
            raweyevec     = getappdata(0,'raweyevec');
            sac           = getappdata(0,'sac');
            vel           = getappdata(0,'vel');
            binsac        = getappdata(0,'binsac');
            PSO           = getappdata(0,'PSO');
            trial         = getappdata(0,'trial');
            eyeparameters = getappdata(0,'eyeparameters');

            % set knob value
            value           = round(app.SmoothingKernelSizeKnob.Value);
            app.KernelSizeLabel.Text = ("Kernel Size: " + num2str(value));
            
            %%%%% SWITCHES %%%%%%
            
            % Moving Mean Switch
            switch app.MovingMeanSwitch.Value
                case 'On'
                    eyevec(1,:) = tn_rmvlinenoise(eyevec(1,:),50,2,eyeparameters.srate);
                    eyevec(2,:) = tn_rmvlinenoise(eyevec(2,:),50,2,eyeparameters.srate);
                    vel         = [];
                    vel(:,1)    = diff(smoothdata(eyevec(1,:),'movmean',5)); % here the original data is smoothed again, because it is also smoothed in the saccade detection
                    vel(:,2)    = diff(smoothdata(eyevec(2,:),'movmean',5));
                case 'Off'
                    vel           = getappdata(0,'vel');
                    eyevec        = getappdata(0,'raweyevec');
                    
            end
            
            % Moving Median
            switch app.MovingMedianSwitch.Value
                case 'On'
                    eyevec(1,:) = tn_rmvlinenoise(eyevec(1,:),50,2,eyeparameters.srate);
                    eyevec(2,:) = tn_rmvlinenoise(eyevec(2,:),50,2,eyeparameters.srate);
                    vel         = [];
                    vel(:,1)    = diff(smoothdata(eyevec(1,:),'movmean',5)); % here the original data is smoothed again, because it is also smoothed in the saccade detection
                    vel(:,2)    = diff(smoothdata(eyevec(2,:),'movmean',5));
                case 'Off'
                    vel           = getappdata(0,'vel');
                    eyevec        = getappdata(0,'raweyevec');
                    
            end
            
            % Abs(Velocity)
            switch app.VelocitySwitch.Value
                case 'On'
                    vel = abs(vel);
                case 'Off'
                    vel           = getappdata(0,'vel');
                    eyevec        = getappdata(0,'raweyevec');
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%   KEEP TRACK OF SACCADE ONSET/OFFSET CONSTRAINTS   %%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %let current sac never get too small or too large
            try
                if  ~any(all(sac)) & ~any(all(PSO))
                    if sac(currentSac,1) < 2
                        sac(currentSac,1) = 2;
                    end
                    if sac(currentSac,2) > size(eyevec(1,:),2)
                        sac(currentSac,2) = size(eyevec(1,:),2);
                    end
                    
                    let saccade end never come before saccade beginning
                    if sac(currentSac,1) == sac(currentSac,2)
                        sac(currentSac,:) = [];
                    end

                    % Make sure to delete PSOs when the duration is zero
                    if PSO(currentSac,1) == PSO(currentSac,2) || PSO(currentSac,2) == sac(currentSac,2)
                        PSO(currentSac,:) = [0 0];
                    end
  
                end
                %%%%%%%%     BINARY VECTOR HOUSEKEEPING     %%%%%%%%%
                
                % Plot the binary vector to show all saccades %%%%
                binsac = zeros(1,size(eyevec,2));
                
                % update the PSO so that the beginning is always the end of
                % the saccade it belongs to
                try
                    if any(PSO(currentSac,:))
                        PSO(currentSac,1) = sac(currentSac,2);
                    end
                    
                    for s = 1:size(sac,1)
                        binsac(sac(s,1):sac(s,2)) = 1;
                        if any(PSO(s,:))
                            binsac(PSO(s,1):PSO(s,2)) = 2;
                        end
                    end
                end
            end
            
            setappdata(0,'sac',sac);
            setappdata(0,'PSO',PSO)
            
            
            
            
            
            try
                % plot the eye traces
                plot(app.EyeTraces,(eyevec(1,:)-mean(eyevec(1,:))),'LineWidth',1.5  ); hold(app.EyeTraces,'on')
                plot(app.EyeTraces,(eyevec(2,:)-mean(eyevec(2,:))),'LineWidth',1.5  );
                
                if ~isempty(sac)
                    yl = app.EyeTraces.YLim; ha = area(app.EyeTraces,[(sac(currentSac,1)) (sac(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'f';
                end
                % plot shaded PSO
                if any(PSO(currentSac,:))
                    yl = app.EyeTraces.YLim; ha = area(app.EyeTraces,[(PSO(currentSac,1)) (PSO(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'g';
                end
                
                
                
                % plot the velocity
                plot(app.Velocities,vel-mean(vel),'LineWidth',1.5 ); hold(app.Velocities,'on')
                if ~isempty(sac)
                    yl = app.Velocities.YLim; ha = area(app.Velocities,[(sac(currentSac,1)) (sac(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'f';
                end
                % plot shaded PSO
                if any(PSO(currentSac,:))
                    yl = app.Velocities.YLim; ha = area(app.Velocities,[(PSO(currentSac,1)) (PSO(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'g';
                end
                
                
                
                %% plot the radial view
                plot(app.RadialView, eyevec(1,:), eyevec(2,:),'LineWidth',1.5 ), hold(app.RadialView,"on")
                plot(app.RadialView, eyevec(1,(sac(currentSac,1):sac(currentSac,2))), eyevec(2,(sac(currentSac,1):sac(currentSac,2))),'r','linewidth',2)
                if any(PSO(currentSac,:))
                    plot(app.RadialView, eyevec(1,(PSO(currentSac,1):PSO(currentSac,2))), eyevec(2,(PSO(currentSac,1):PSO(currentSac,2))),'g','linewidth',2)
                end
                
                % plot the acceleration
                
                plot(app.Acceleration,diff(vel)-mean(diff(vel)),'linewidth',2); hold(app.Acceleration,'on')
                app.Acceleration.XLim = [[(sac(currentSac,1))-20, (sac(currentSac,1))+70]];
                if ~isempty(sac)
                    xline(app.Acceleration,sac(currentSac,1),'r','LineWidth',1.5 );
                    xline(app.Acceleration,sac(currentSac,2),'b','LineWidth',1.5 )
                end
                
                if ~isempty(sac)
                    yl = app.Acceleration.YLim; ha = area(app.Acceleration,[(sac(currentSac,1)) (sac(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'f';
                end
                % plot shaded PSO
                if any(PSO(currentSac,:))
                    yl = app.Acceleration.YLim; ha = area(app.Acceleration,[(PSO(currentSac,1)) (PSO(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'g';
                end
                
                
                
                
                %%%%% PLOT ZOOMED EYE TRACES %%%%%%
                mymean1 = mean(eyevec(1,:)); mymean2 = mean(eyevec(2,:)); % mean over whole vector
                try
                    mymean1 = mean(eyevec(1,(sac(currentSac,1))-20:(sac(currentSac,1))+70));
                    mymean2 = mean(eyevec(1,(sac(currentSac,1))-20:(sac(currentSac,1))+70));
                end
                plot(app.ZoomTrace,(eyevec(1,:)-mymean1),'LineWidth',1.5  ); hold(app.ZoomTrace,'on')
                plot(app.ZoomTrace,(eyevec(2,:)-mymean2),'LineWidth',1.5  );
                app.ZoomTrace.XLim = [(sac(currentSac,1))-20, (sac(currentSac,1))+70];
                
                
                if ~isempty(sac)
                    xline(app.ZoomTrace,sac(currentSac,1),'r','LineWidth',1.5 );
                    xline(app.ZoomTrace,sac(currentSac,2),'b','LineWidth',1.5 )
                end
                if ~isempty(sac)
                    yl = app.ZoomTrace.YLim; ha = area(app.ZoomTrace,[(sac(currentSac,1)) (sac(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'f';
                end
                % plot shaded PSO
                if any(PSO(currentSac,:))
                    yl = app.ZoomTrace.YLim; ha = area(app.ZoomTrace,[(PSO(currentSac,1)) (PSO(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'g';
                end
                
                
                %%%%% PLOT ZOOMED IN VELOCITY %%%%
                
                plot(app.ZoomVelocity,vel,'LineWidth',1.5 ); hold(app.ZoomVelocity,'on')
                app.ZoomVelocity.XLim = [[(sac(currentSac,1))-20, (sac(currentSac,1))+70]];
                if ~isempty(sac)
                    xline(app.ZoomVelocity,sac(currentSac,1),'r','LineWidth',1.5 );
                    xline(app.ZoomVelocity,sac(currentSac,2),'b','LineWidth',1.5 )
                end
                if ~isempty(sac)
                    yl = app.ZoomVelocity.YLim; ha = area(app.ZoomVelocity,[(sac(currentSac,1)) (sac(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'f';
                end
                % plot shaded PSO
                if any(PSO(currentSac,:))
                    yl = app.ZoomVelocity.YLim; ha = area(app.ZoomVelocity,[(PSO(currentSac,1)) (PSO(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'g';
                end
                
                
                
                
                
                %%%% PLOT BINARY VEC %%%%%
                
                plot(app.BinSacs,binsac,'LineWidth',1.5 ); hold(app.BinSacs,"on")
                % plot shaded saccade
                if ~isempty(sac)
                    yl = app.BinSacs.YLim; ha = area(app.BinSacs,[(sac(currentSac,1)) (sac(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'f';
                end
                
                % plot shaded PSO
                if any(PSO(currentSac,:))
                    yl = app.BinSacs.YLim; ha = area(app.BinSacs,[(PSO(currentSac,1)) (PSO(currentSac,2))], [yl(2) yl(2)], yl(1));
                    ha.FaceAlpha = 0.1; ha.EdgeAlpha = 0.1; ha.FaceColor = 'g';
                end
            
            if ~isempty(sac)
                app.Label2.Text = ("Current saccade duration: " + num2str(sac(currentSac,2)-sac(currentSac,1)));
            else
                app.Label2.Text = "No Saccades detected!";
            end
            end
            
            if trial == 10
                uiwait(msgbox('Good job for labeling 10 trials! Keep going :)'));
                figure(app.UIFigure)
            end
           
        end
    end  
 

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)

            figure(app.UIFigure)
            selpath = uigetdir; % select folder where the X and Y data are stored
            figure(app.UIFigure)
            cd(selpath)
            
            if exist("outputdir", 'dir') == 0
                mkdir("outputdir")
            end
            
            cd('outputdir'); % create output directory
            dlmwrite('X_coords.csv',[],'delimiter',','); % create output files
            dlmwrite('Y_coords.csv',[],'delimiter',',');
            dlmwrite('labels.csv',  [],'delimiter',',');
            cd(selpath) % go back to working memory files

            
            % open dialog box
            prompt    = {'Enter sample rate [Hz]','Enter minimum duration for saccades:', 'Enter SD threshold (VFAC):',...
                        'Make saccade suggestions on smoothed data (1 = true, 0 = false)?:', 'Line Noise Frequency (Hz):'};
            dlgtitle  = 'Input';
            dims      = [1 55];
            definput  = {'1000','8','5','1','50'};
            answer    = inputdlg(prompt,dlgtitle,dims,definput);
            srate     = str2double(answer{1});
            MINDUR    = str2double(answer{2});
            VFAC      = str2double(answer{3});
            SMOOTHDAT = str2double(answer{4});
            LNFREQ    = str2double(answer{5});
            
            
            eyeparameters.rmvlinenoise = true;
            eyeparameters.VFAC         = VFAC;
            eyeparameters.MINDUR       = MINDUR;
            eyeparameters.srate        = srate;
            eyeparameters.mergeint     = 50; % saccades separater less than mergeint are merged into 1 event
            eyeparameters.slength      = 1; % kernel for velocities
            eyeparameters.smoothdata   = SMOOTHDAT;
            eyeparameters.LNFREQ       = LNFREQ;
            
            % get filenames in folder
            D      = selpath;
            T      = dir(fullfile(D,'*.csv'));
            C      = {T(~[T.isdir]).name};
            xname  = C{find(contains(C,'X'))}; % find the files that contain the strings X and Y
            yname  = C{find(contains(C,'Y'))};
            xdata  = load(xname);
            ydata  = load(yname); % load the files
            
            % detect MSs here
            trial       = 1;
            eyevec      = [xdata(trial,:) ; ydata(trial,:)];
            raweyevec   = [xdata(trial,:) ; ydata(trial,:)];
            [sac, vel]  = msdetect(raweyevec,eyeparameters);
            
            
            binsac = zeros(1,size(eyevec,2));
            for s = 1:size(sac,1)
                binsac(sac(s,1):sac(s,2)) = 1;
            end
            
            %% set the app data
            currentSac = 1;
            currentPSO = 1;
            PSO        = zeros(size(sac(:,1:2))); %init PSO times
            setappdata(0,'currentSac',currentSac)
            setappdata(0,'currentPSO',currentPSO)
            setappdata(0,'eyeparameters',eyeparameters)
            setappdata(0,'trial',trial);
            setappdata(0,'sac',sac);
            setappdata(0,'vel',vel);
            setappdata(0,'xdata',xdata);
            setappdata(0,'ydata',ydata);
            setappdata(0,'eyevec',eyevec);
            setappdata(0,'raweyevec',raweyevec);
            setappdata(0,'binsac',binsac);
            setappdata(0,'PSO',PSO)
            
            
            % Set labs
            value                     = round(app.SmoothingKernelSizeKnob.Value);
            app.KernelSizeLabel.Text  = ("Kernel Size: " + num2str(value));
            app.Label.Text            = ("Current trial: " + num2str(trial)); % show current trial
            app.Label2.Text           = ("Current saccade duration: " + num2str(sac(currentSac,2)-sac(currentSac,1)));
            
            updatePlot(app,currentSac)
        end

        % Button pushed function: NextSaccadeButton
        function NextSaccadeButtonPushed(app, event)
            
            % get shared variables
            currentSac = getappdata(0,'currentSac');
            sac        = getappdata(0,'sac');
            
            currentSac = currentSac + 1;
            if currentSac > size(sac,1); currentSac = size(sac,1) ; end % sets max for current sac
            updatePlot(app,currentSac)
            
            % update the global data
            setappdata(0,'currentSac',currentSac)
        end

        % Button pushed function: PreviousSaccadeButton
        function PreviousSaccadeButtonPushed(app, event)
            
            % get shared variables
            currentSac = getappdata(0,'currentSac');
            currentSac = currentSac - 1;
            
            if currentSac == 0; currentSac = 1 ; end % sets min for current sac
            updatePlot(app,currentSac)
            
            % update the global data
            setappdata(0,'currentSac',currentSac);
        end

        % Button pushed function: NextTrialButton
        function NextTrialButtonPushed(app, event)
            
            % get shared variables
            eyeparameters = getappdata(0,'eyeparameters');
            trial         = getappdata(0,'trial');
            xdata         = getappdata(0,'xdata');
            ydata         = getappdata(0,'ydata');

            trial      = trial + 1; % count another trial
            currentSac = 1; % reset current saccade to 1;
            
            % update the eyevec and raw eyevec vectors
            eyevec      = [xdata(trial,:) ; ydata(trial,:)];
            raweyevec   = [xdata(trial,:) ; ydata(trial,:)];
            [sac, vel]  = msdetect(raweyevec,eyeparameters);
            
            % init PSO
            if isempty(sac)
                PSO = zeros(1,2);
            else
                PSO = zeros(size(sac(:,1:2))); %init PSO times
            end

            
            % update the global data
            setappdata(0,'currentSac',currentSac);
            setappdata(0,'sac',sac);
            setappdata(0,'vel',vel);
            setappdata(0,'eyevec',eyevec);
            setappdata(0,'raweyevec',raweyevec);
            setappdata(0,'trial',trial);
            setappdata(0,'PSO',PSO)
            
            % Reset the switches
            app.MovingMeanSwitch.Value = 'Off';
            app.MovingMedianSwitch.Value    = 'Off';
            app.VelocitySwitch.Value        = 'Off';

            updatePlot(app,currentSac)
            app.Label.Text = ("Current trial: " + num2str(trial)); % show current trial
            
        end

        % Button pushed function: PreviousTrialButton
        function PreviousTrialButtonPushed(app, event)
            % get shared variables
            eyeparameters = getappdata(0,'eyeparameters');
            trial         = getappdata(0,'trial');
            xdata         = getappdata(0,'xdata');
            ydata         = getappdata(0,'ydata');
            
            trial      = trial - 1; % count another trial
            if trial == 0; trial = 1; end % make sure current trial cant go to zero
            currentSac = 1; % reset current saccade to 1;
            
            % update the eyevec and raw eyevec vectors
            eyevec    = [xdata(trial,:) ; ydata(trial,:)];
            raweyevec = [xdata(trial,:) ; ydata(trial,:)];
            [sac, vel]  = msdetect(raweyevec,eyeparameters);
            
            % init PSO
            if isempty(sac)
                PSO = zeros(1,2);
            else
                PSO = zeros(size(sac(:,1:2))); %init PSO times
            end
            
            % update the global data
            setappdata(0,'currentSac',currentSac);
            setappdata(0,'sac',sac);
            setappdata(0,'vel',vel);
            setappdata(0,'eyevec',eyevec);
            setappdata(0,'raweyevec',raweyevec);
            setappdata(0,'trial',trial);
            setappdata(0,'PSO',PSO)
            
            % Reset the switches
            app.MovingMeanSwitch.Value = 'Off';
            app.MovingMedianSwitch.Value    = 'Off';
            app.VelocitySwitch.Value        = 'Off';
            
            updatePlot(app,currentSac)
            app.Label.Text = ("Current trial: " + num2str(trial)); % show current trial
        end

        % Button pushed function: OnsetRight
        function OnsetRightButtonPushed(app, event)
            currentSac    = getappdata(0,'currentSac');
            sac           = getappdata(0,'sac');
            
            % move the saccade to the right
            sac(currentSac,1) = sac(currentSac,1) + 1;
            setappdata(0,'sac',sac);
            
            updatePlot(app,currentSac)
        end

        % Button pushed function: OnsetLeft
        function OnsetLeftButtonPushed(app, event)
            currentSac    = getappdata(0,'currentSac');
            sac           = getappdata(0,'sac');
            
            % move the saccade to the right
            sac(currentSac,1) = sac(currentSac,1) - 1;
            setappdata(0,'sac',sac);
            
            updatePlot(app,currentSac)
            
        end

        % Button pushed function: OffsetRight
        function OffsetRightButtonPushed(app, event)
            currentSac    = getappdata(0,'currentSac');
            sac           = getappdata(0,'sac');
            
            % move the saccade to the right
            sac(currentSac,2) = sac(currentSac,2) + 1;
            setappdata(0,'sac',sac);
            
            updatePlot(app,currentSac)
        end

        % Button pushed function: OffsetLeft
        function OffsetLeftButtonPushed(app, event)
            currentSac    = getappdata(0,'currentSac');
            sac           = getappdata(0,'sac');
            
            % move the saccade to the right
            sac(currentSac,2) = sac(currentSac,2) - 1;
            setappdata(0,'sac',sac);
            
            updatePlot(app,currentSac)
        end

        % Button pushed function: NewSaccadeButton
        function NewSaccadeButtonPushed(app, event)
            currentSac    = getappdata(0,'currentSac');
            sac           = getappdata(0,'sac');
            PSO           = getappdata(0,'PSO');
            eyevec        = getappdata(0,'eyevec');
            
            if (app.OnsetEditField.Value <  app.OffsetEditField.Value) && (app.OnsetEditField.Value > 1) && (app.OffsetEditField.Value < size(eyevec,2))

                new     = [app.OnsetEditField.Value app.OffsetEditField.Value  0 0 0 0 0 ];
                newPSO  = [0 0];
                
                sac     = vertcat(sac, new);
                sac     = sort(sac,1);
                PSO     = vertcat(PSO,newPSO);
                PSO     = sort(PSO,1);

                setappdata(0,'sac',sac);
                setappdata(0,'PSO',PSO);

            else
                 uialert(app.UIFigure,'Please select valid start and end points!','Warning')
            end
            
            updatePlot(app,currentSac)
            
        end

        % Button pushed function: DeleteCurrentSaccadeButton
        function DeleteCurrentSaccadeButtonPushed(app, event)

            currentSac    = getappdata(0,'currentSac');
            sac           = getappdata(0,'sac');
            PSO           = getappdata(0,'PSO');
            
            sac(currentSac,:) = [ ];
            PSO(currentSac,:) = [ ];
            currentSac        = 1;
            setappdata(0,'sac',sac);
            setappdata(0,'PSO',PSO);
            setappdata(0,'currentSac',currentSac);
            
            updatePlot(app,currentSac)
        end

        % Button pushed function: SaveCurrentTrialButton
        function SaveCurrentTrialButtonPushed(app, event)
            sac           = getappdata(0,'sac');
            eyevec        = getappdata(0,'eyevec');
            trial         = getappdata(0,'trial');
            xdata         = getappdata(0,'xdata');
            ydata         = getappdata(0,'ydata');
            
            
            % create the binary output vector again
            binsac = zeros(1,size(eyevec,2));
            for s = 1:size(sac,1)
                binsac(sac(s,1):sac(s,2)) = 1;
            end
            
            cd('outputdir') % cd to the variable directory
            
            % write trials in X and Y coordinates to the matricies. Also
            % append the binary labels to the matrix labels
            dlmwrite('X_coords.csv',eyevec(1,:),'delimiter',',','-append');
            dlmwrite('Y_coords.csv',eyevec(2,:),'delimiter',',','-append');
            dlmwrite('labels.csv',  binsac,     'delimiter',',','-append');
            
            
            % write individual files if necessary
            xfname = ['X_trial_' num2str(trial) '.csv'];
            yfname = ['Y_trial_' num2str(trial) '.csv'];
            lfname = ['L_trial_' num2str(trial) '.csv'];
            
            dlmwrite(xfname,eyevec(1,:));
            dlmwrite(yfname,eyevec(2,:));
            dlmwrite(lfname,binsac);
            
            % create to do matricies that include only the trials which
            % have not been labeled yet
            X_todo = xdata(trial:end,:);
            Y_todo = ydata(trial:end,:);
            
            save('X_todo.mat','X_todo');
            save('Y_todo.mat','Y_todo');
            
            uiwait(msgbox('Trial saved and appended'));
            
            
            cd .. % cd back to directory
        end

        % Button pushed function: MakeBrushedSaccadesButton
        function MakeBrushedSaccadesButtonPushed(app, event)
            currentSac    = getappdata(0,'currentSac');
            sac           = getappdata(0,'sac');
            
            figure_handles  = findall(app.ZoomTrace);
            Handles         = figure_handles( arrayfun(@(H) isprop(H, 'BrushData'), figure_handles) );
            info            = [num2cell(Handles),arrayfun(@(H) find(H.BrushData), Handles, 'uniform', 0)];
          
            
            beg_indx           = min([info{2,2}  info{3,2}]); % smallest index
            end_indx           = max([info{2,2}  info{3,2}]); % largest index
            sac(currentSac,:)  = zeros(1,7);
            sac(~any(sac')',:) = [];
            
            
            try
                new = [beg_indx end_indx  0 0 0 0 0 ];
                sac = vertcat(sac, new);
                sac = sort(sac,1);
                
                
                % check if 2 saccades overlap and potentially merge
                
                differences = diff(reshape(sac(:,1:2)',1,[]));
                differences = differences(2:2:end);
                
                for j = find(differences < 0)
                    sac(j,1)   = min([sac(j,1) sac(j+1,1)]);
                    sac(j,2)   = max([sac(j,2) sac(j+1,2)]);
                    sac(j+1,:) = [];
                end
                
                
                setappdata(0,'sac',sac);
            catch
                uialert(app.UIFigure,'Please brush some data','Warning')
            end
            
            updatePlot(app,currentSac)
        end

        % Button pushed function: ActivateBrushButton
        function ActivateBrushButtonPushed(app, event)
            brush(app.ZoomTrace,'on')
        end

        % Button pushed function: MakeBrushedPSOButton
        function MakeBrushedPSOButtonPushed(app, event)
            currentSac    = getappdata(0,'currentSac');
            PSO           = getappdata(0,'PSO');
            sac           = getappdata(0,'sac');
            
            figure_handles = findall(app.ZoomTrace);
            Handles        = figure_handles( arrayfun(@(H) isprop(H, 'BrushData'), figure_handles) );
            info           = [num2cell(Handles),arrayfun(@(H) find(H.BrushData), Handles, 'uniform', 0)];
            
            beg_indx           = min([info{2,2}  info{3,2}]); % smallest index
            end_indx           = max([info{2,2}  info{3,2}]); % largest index
            
            try
                PSO(currentSac,:)  = zeros(1,2);  % delete current PSO/overwrite it
                PSO(currentSac,:) = [beg_indx end_indx];
                % in this case, where the PSO is manually overwritten,
                % overwrite the saccade end
                sac(currentSac,2) = PSO(currentSac,1);
                
                
                
                setappdata(0,'PSO',PSO); % update the data
                setappdata(0,'sac',sac); % update the data
            catch
                uialert(app.UIFigure,'Please brush some data','Warning')
            end
            
            
            
            updatePlot(app,currentSac)
        end

        % Callback function
        function PreviousPSOButtonPushed(app, event)
            
        end

        % Callback function
        function NextPSOButtonPushed(app, event)
            
        end

        % Button pushed function: AddPSOtocurrentSaccadeButton
        function AddPSOtocurrentSaccadeButtonPushed(app, event)
            currentSac    = getappdata(0,'currentSac');
            sac           = getappdata(0,'sac');
            PSO           = getappdata(0,'PSO');

            
            PSObeg       = sac(currentSac,2)+1;
            PSOend       = sac(currentSac,2)+11;
            
            PSO(currentSac,1) = PSObeg; % create the label 2 as PSO labels in the binsac vector
            PSO(currentSac,2) = PSOend; % create the label 2 as PSO labels in the binsac vector
            
            

            setappdata(0,'PSO',PSO); % update the data
           
            updatePlot(app,currentSac)
            
        end

        % Button pushed function: PSOOffsetRight
        function PSOOffsetRightButtonPushed(app, event)
            currentSac    = getappdata(0,'currentSac');
            PSO           = getappdata(0,'PSO');
            
            if any(PSO(currentSac,:))
                PSO(currentSac,2) = PSO(currentSac,2) + 1;
            end

            
            setappdata(0,'PSO',PSO); % update the data
            
            updatePlot(app,currentSac)
       
        end

        % Button pushed function: PSOOffsetLeft
        function PSOOffsetLeftButtonPushed(app, event)
            currentSac    = getappdata(0,'currentSac');
            PSO           = getappdata(0,'PSO');
            
            if any(PSO(currentSac,:))
                PSO(currentSac,2) = PSO(currentSac,2) - 1;
            end
   
            setappdata(0,'PSO',PSO); % update the data
            
            updatePlot(app,currentSac)
            
        end

        % Button pushed function: DeleteCurrentPSOButton
        function DeleteCurrentPSOButtonPushed(app, event)
            currentSac    = getappdata(0,'currentSac');
            PSO           = getappdata(0,'PSO');
            
        
            PSO(currentSac,:) = [0 0];
            
            setappdata(0,'PSO',PSO);
            setappdata(0,'currentSac',1);
            
            updatePlot(app,currentSac)
        end

        % Value changed function: MovingMeanSwitch
        function MovingMeanSwitchValueChanged(app, event)
            eyevec        = getappdata(0,'eyevec');
            raweyevec     = getappdata(0,'raweyevec');
            currentSac    = getappdata(0,'currentSac');
            knobValue     = app.SmoothingKernelSizeKnob.Value;

            
%             switch app.MovingMeanSwitch.Value
%                 case 'On'
%                     eyevec(1,:) = movmean(eyevec(1,:),uint8(knobValue));
%                     eyevec(2,:) = movmean(eyevec(2,:),uint8(knobValue));
%                 case 'Off'
%                     eyevec = raweyevec;
%             end
            
            
            %setappdata(0,'eyevec',eyevec);
            updatePlot(app,currentSac)
        end

        % Value changed function: MovingMedianSwitch
        function MovingMedianSwitchValueChanged(app, event)
            eyevec        = getappdata(0,'eyevec');
            raweyevec     = getappdata(0,'raweyevec');
            currentSac    = getappdata(0,'currentSac');
            knobValue     = app.SmoothingKernelSizeKnob.Value;

            
%             switch app.MovingMedianSwitch.Value
%                 case 'On'
%                     eyevec(1,:) = movmedian(eyevec(1,:),uint8(knobValue));
%                     eyevec(2,:) = movmedian(eyevec(2,:),uint8(knobValue));
%                 case 'Off'
%                     eyevec = raweyevec;
%             end
            
            
            setappdata(0,'eyevec',eyevec);
            updatePlot(app,currentSac)            
        end

        % Value changed function: VelocitySwitch
        function VelocitySwitchValueChanged(app, event)
            currentSac    = getappdata(0,'currentSac');
            
            updatePlot(app,currentSac);
        end

        % Value changed function: SmoothingKernelSizeKnob
        function SmoothingKernelSizeKnobValueChanged(app, event)
 
            currentSac                      = getappdata(0,'currentSac');
            app.MovingMeanSwitch.Value      = 'Off';
            app.MovingMedianSwitch.Value    = 'Off';
            app.VelocitySwitch.Value        = 'Off';
            
            updatePlot(app,currentSac)
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1406 884];
            app.UIFigure.Name = 'MATLAB App';

            % Create EyeTraces
            app.EyeTraces = uiaxes(app.UIFigure);
            title(app.EyeTraces, 'Eye Traces')
            xlabel(app.EyeTraces, 'time [s]')
            ylabel(app.EyeTraces, '')
            app.EyeTraces.PlotBoxAspectRatio = [21.5555555555556 1 1];
            app.EyeTraces.YTick = [];
            app.EyeTraces.YTickLabel = '';
            app.EyeTraces.BackgroundColor = [1 1 1];
            app.EyeTraces.Clipping = 'off';
            app.EyeTraces.Position = [1 703 940 145];

            % Create Velocities
            app.Velocities = uiaxes(app.UIFigure);
            title(app.Velocities, 'Velocity')
            xlabel(app.Velocities, '')
            ylabel(app.Velocities, '')
            app.Velocities.PlotBoxAspectRatio = [22.3846153846154 1 1];
            app.Velocities.XTick = [];
            app.Velocities.XTickLabel = '';
            app.Velocities.YTick = [];
            app.Velocities.ZTick = [];
            app.Velocities.BackgroundColor = [1 1 1];
            app.Velocities.Position = [1 591 940 113];

            % Create RadialView
            app.RadialView = uiaxes(app.UIFigure);
            title(app.RadialView, '2-D View')
            xlabel(app.RadialView, '')
            ylabel(app.RadialView, '')
            app.RadialView.PlotBoxAspectRatio = [1.16666666666667 1 1];
            app.RadialView.Box = 'on';
            app.RadialView.XTick = [];
            app.RadialView.YTick = [];
            app.RadialView.BackgroundColor = [1 1 1];
            app.RadialView.Position = [947 495 443 357];

            % Create ControlsPanel
            app.ControlsPanel = uipanel(app.UIFigure);
            app.ControlsPanel.TitlePosition = 'centertop';
            app.ControlsPanel.Title = 'Controls';
            app.ControlsPanel.Position = [19 12 267 484];

            % Create NextSaccadeButton
            app.NextSaccadeButton = uibutton(app.ControlsPanel, 'push');
            app.NextSaccadeButton.ButtonPushedFcn = createCallbackFcn(app, @NextSaccadeButtonPushed, true);
            app.NextSaccadeButton.Position = [144 422 111 22];
            app.NextSaccadeButton.Text = 'Next Saccade';

            % Create PreviousSaccadeButton
            app.PreviousSaccadeButton = uibutton(app.ControlsPanel, 'push');
            app.PreviousSaccadeButton.ButtonPushedFcn = createCallbackFcn(app, @PreviousSaccadeButtonPushed, true);
            app.PreviousSaccadeButton.Position = [22 422 111 22];
            app.PreviousSaccadeButton.Text = 'Previous Saccade';

            % Create OnsetLeft
            app.OnsetLeft = uibutton(app.ControlsPanel, 'push');
            app.OnsetLeft.ButtonPushedFcn = createCallbackFcn(app, @OnsetLeftButtonPushed, true);
            app.OnsetLeft.Position = [19 363 100 22];
            app.OnsetLeft.Text = '<';

            % Create OnsetRight
            app.OnsetRight = uibutton(app.ControlsPanel, 'push');
            app.OnsetRight.ButtonPushedFcn = createCallbackFcn(app, @OnsetRightButtonPushed, true);
            app.OnsetRight.Position = [144 360 100 22];
            app.OnsetRight.Text = '>';

            % Create SaccadeOnsetLabel
            app.SaccadeOnsetLabel = uilabel(app.ControlsPanel);
            app.SaccadeOnsetLabel.Position = [90 384 88 22];
            app.SaccadeOnsetLabel.Text = 'Saccade Onset';

            % Create OffsetRight
            app.OffsetRight = uibutton(app.ControlsPanel, 'push');
            app.OffsetRight.ButtonPushedFcn = createCallbackFcn(app, @OffsetRightButtonPushed, true);
            app.OffsetRight.Position = [144 292 100 22];
            app.OffsetRight.Text = '>';

            % Create OffsetLeft
            app.OffsetLeft = uibutton(app.ControlsPanel, 'push');
            app.OffsetLeft.ButtonPushedFcn = createCallbackFcn(app, @OffsetLeftButtonPushed, true);
            app.OffsetLeft.Position = [19 294 100 22];
            app.OffsetLeft.Text = '<';

            % Create SaccadeOffsetLabel
            app.SaccadeOffsetLabel = uilabel(app.ControlsPanel);
            app.SaccadeOffsetLabel.Position = [90 317 88 22];
            app.SaccadeOffsetLabel.Text = 'Saccade Offset';

            % Create DeleteCurrentSaccadeButton
            app.DeleteCurrentSaccadeButton = uibutton(app.ControlsPanel, 'push');
            app.DeleteCurrentSaccadeButton.ButtonPushedFcn = createCallbackFcn(app, @DeleteCurrentSaccadeButtonPushed, true);
            app.DeleteCurrentSaccadeButton.Position = [58 78 143 22];
            app.DeleteCurrentSaccadeButton.Text = {'Delete Current Saccade'; ''};

            % Create AddPSOtocurrentSaccadeButton
            app.AddPSOtocurrentSaccadeButton = uibutton(app.ControlsPanel, 'push');
            app.AddPSOtocurrentSaccadeButton.ButtonPushedFcn = createCallbackFcn(app, @AddPSOtocurrentSaccadeButtonPushed, true);
            app.AddPSOtocurrentSaccadeButton.Position = [44 26 170 22];
            app.AddPSOtocurrentSaccadeButton.Text = 'Add PSO to current Saccade';

            % Create PSOOffsetLeft
            app.PSOOffsetLeft = uibutton(app.ControlsPanel, 'push');
            app.PSOOffsetLeft.ButtonPushedFcn = createCallbackFcn(app, @PSOOffsetLeftButtonPushed, true);
            app.PSOOffsetLeft.Position = [19 209 100 22];
            app.PSOOffsetLeft.Text = '<';

            % Create PSOOffsetRight
            app.PSOOffsetRight = uibutton(app.ControlsPanel, 'push');
            app.PSOOffsetRight.ButtonPushedFcn = createCallbackFcn(app, @PSOOffsetRightButtonPushed, true);
            app.PSOOffsetRight.Position = [144 209 100 22];
            app.PSOOffsetRight.Text = '>';

            % Create PSOOffsetLabel
            app.PSOOffsetLabel = uilabel(app.ControlsPanel);
            app.PSOOffsetLabel.Position = [101 240 66 22];
            app.PSOOffsetLabel.Text = 'PSO Offset';

            % Create DeleteCurrentPSOButton
            app.DeleteCurrentPSOButton = uibutton(app.ControlsPanel, 'push');
            app.DeleteCurrentPSOButton.ButtonPushedFcn = createCallbackFcn(app, @DeleteCurrentPSOButtonPushed, true);
            app.DeleteCurrentPSOButton.Position = [69 125 120 22];
            app.DeleteCurrentPSOButton.Text = 'Delete Current PSO';

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.Position = [1096 428 258 33];

            % Create BinSacs
            app.BinSacs = uiaxes(app.UIFigure);
            title(app.BinSacs, '')
            xlabel(app.BinSacs, '')
            ylabel(app.BinSacs, '')
            app.BinSacs.PlotBoxAspectRatio = [18.4471544715447 1 1];
            app.BinSacs.XTick = [];
            app.BinSacs.YTick = [];
            app.BinSacs.BackgroundColor = [1 1 1];
            app.BinSacs.Position = [1 505 940 87];

            % Create ZoomVelocity
            app.ZoomVelocity = uiaxes(app.UIFigure);
            title(app.ZoomVelocity, 'Zoom Velocity')
            xlabel(app.ZoomVelocity, '')
            ylabel(app.ZoomVelocity, '')
            app.ZoomVelocity.Box = 'on';
            app.ZoomVelocity.XTick = [];
            app.ZoomVelocity.YTick = [];
            app.ZoomVelocity.BackgroundColor = [1 1 1];
            app.ZoomVelocity.Position = [303 12 371 218];

            % Create ZoomTrace
            app.ZoomTrace = uiaxes(app.UIFigure);
            title(app.ZoomTrace, 'Zoom Trace')
            xlabel(app.ZoomTrace, '')
            ylabel(app.ZoomTrace, '')
            app.ZoomTrace.Box = 'on';
            app.ZoomTrace.XTick = [];
            app.ZoomTrace.YTick = [];
            app.ZoomTrace.BackgroundColor = [1 1 1];
            app.ZoomTrace.Position = [303 303 371 208];

            % Create CreatenewSaccadesPanel
            app.CreatenewSaccadesPanel = uipanel(app.UIFigure);
            app.CreatenewSaccadesPanel.TitlePosition = 'centertop';
            app.CreatenewSaccadesPanel.Title = 'Create new Saccades';
            app.CreatenewSaccadesPanel.Position = [704 371 363 118];

            % Create NewSaccadeButton
            app.NewSaccadeButton = uibutton(app.CreatenewSaccadesPanel, 'push');
            app.NewSaccadeButton.ButtonPushedFcn = createCallbackFcn(app, @NewSaccadeButtonPushed, true);
            app.NewSaccadeButton.Position = [237 27 77 60];
            app.NewSaccadeButton.Text = 'Create!';

            % Create OnsetEditFieldLabel
            app.OnsetEditFieldLabel = uilabel(app.CreatenewSaccadesPanel);
            app.OnsetEditFieldLabel.HorizontalAlignment = 'right';
            app.OnsetEditFieldLabel.Position = [32 65 37 22];
            app.OnsetEditFieldLabel.Text = 'Onset';

            % Create OnsetEditField
            app.OnsetEditField = uieditfield(app.CreatenewSaccadesPanel, 'numeric');
            app.OnsetEditField.Position = [84 65 100 22];

            % Create OffsetEditFieldLabel
            app.OffsetEditFieldLabel = uilabel(app.CreatenewSaccadesPanel);
            app.OffsetEditFieldLabel.HorizontalAlignment = 'right';
            app.OffsetEditFieldLabel.Position = [32 25 37 22];
            app.OffsetEditFieldLabel.Text = 'Offset';

            % Create OffsetEditField
            app.OffsetEditField = uieditfield(app.CreatenewSaccadesPanel, 'numeric');
            app.OffsetEditField.Position = [84 25 100 22];

            % Create MakeBrushedSaccadesButton
            app.MakeBrushedSaccadesButton = uibutton(app.UIFigure, 'push');
            app.MakeBrushedSaccadesButton.ButtonPushedFcn = createCallbackFcn(app, @MakeBrushedSaccadesButtonPushed, true);
            app.MakeBrushedSaccadesButton.Position = [524 275 142 22];
            app.MakeBrushedSaccadesButton.Text = 'Make Brushed Saccades';

            % Create ActivateBrushButton
            app.ActivateBrushButton = uibutton(app.UIFigure, 'push');
            app.ActivateBrushButton.ButtonPushedFcn = createCallbackFcn(app, @ActivateBrushButtonPushed, true);
            app.ActivateBrushButton.Position = [323 275 99.99 22];
            app.ActivateBrushButton.Text = 'Activate Brush';

            % Create MakeBrushedPSOButton
            app.MakeBrushedPSOButton = uibutton(app.UIFigure, 'push');
            app.MakeBrushedPSOButton.ButtonPushedFcn = createCallbackFcn(app, @MakeBrushedPSOButtonPushed, true);
            app.MakeBrushedPSOButton.Position = [524 244 142 22];
            app.MakeBrushedPSOButton.Text = 'Make Brushed PSO';

            % Create TrialNavigationPanel
            app.TrialNavigationPanel = uipanel(app.UIFigure);
            app.TrialNavigationPanel.TitlePosition = 'centertop';
            app.TrialNavigationPanel.Title = 'Trial Navigation';
            app.TrialNavigationPanel.Position = [704 241 363 109];

            % Create PreviousTrialButton
            app.PreviousTrialButton = uibutton(app.TrialNavigationPanel, 'push');
            app.PreviousTrialButton.ButtonPushedFcn = createCallbackFcn(app, @PreviousTrialButtonPushed, true);
            app.PreviousTrialButton.Position = [52 13 111 22];
            app.PreviousTrialButton.Text = 'Previous Trial';

            % Create NextTrialButton
            app.NextTrialButton = uibutton(app.TrialNavigationPanel, 'push');
            app.NextTrialButton.ButtonPushedFcn = createCallbackFcn(app, @NextTrialButtonPushed, true);
            app.NextTrialButton.Position = [203 13 111 22];
            app.NextTrialButton.Text = 'Next Trial';

            % Create SaveCurrentTrialButton
            app.SaveCurrentTrialButton = uibutton(app.TrialNavigationPanel, 'push');
            app.SaveCurrentTrialButton.ButtonPushedFcn = createCallbackFcn(app, @SaveCurrentTrialButtonPushed, true);
            app.SaveCurrentTrialButton.Position = [125 55 110 22];
            app.SaveCurrentTrialButton.Text = 'Save Current Trial';

            % Create FilterBankPanel
            app.FilterBankPanel = uipanel(app.UIFigure);
            app.FilterBankPanel.TitlePosition = 'centertop';
            app.FilterBankPanel.Title = 'Filter Bank';
            app.FilterBankPanel.Position = [1094 27 296 348];

            % Create MovingMeanSwitchLabel
            app.MovingMeanSwitchLabel = uilabel(app.FilterBankPanel);
            app.MovingMeanSwitchLabel.HorizontalAlignment = 'center';
            app.MovingMeanSwitchLabel.Position = [3 266 141 22];
            app.MovingMeanSwitchLabel.Text = 'Moving Mean';

            % Create MovingMeanSwitch
            app.MovingMeanSwitch = uiswitch(app.FilterBankPanel, 'slider');
            app.MovingMeanSwitch.ValueChangedFcn = createCallbackFcn(app, @MovingMeanSwitchValueChanged, true);
            app.MovingMeanSwitch.Position = [51 303 45 20];

            % Create MovingMedianSwitchLabel
            app.MovingMedianSwitchLabel = uilabel(app.FilterBankPanel);
            app.MovingMedianSwitchLabel.HorizontalAlignment = 'center';
            app.MovingMedianSwitchLabel.Position = [5 195 128 22];
            app.MovingMedianSwitchLabel.Text = 'Moving Median';

            % Create MovingMedianSwitch
            app.MovingMedianSwitch = uiswitch(app.FilterBankPanel, 'slider');
            app.MovingMedianSwitch.ValueChangedFcn = createCallbackFcn(app, @MovingMedianSwitchValueChanged, true);
            app.MovingMedianSwitch.Position = [40 232 58 25];

            % Create VelocityLabel
            app.VelocityLabel = uilabel(app.FilterBankPanel);
            app.VelocityLabel.HorizontalAlignment = 'center';
            app.VelocityLabel.Position = [201 193 59 22];
            app.VelocityLabel.Text = '| Velocity |';

            % Create VelocitySwitch
            app.VelocitySwitch = uiswitch(app.FilterBankPanel, 'rocker');
            app.VelocitySwitch.ValueChangedFcn = createCallbackFcn(app, @VelocitySwitchValueChanged, true);
            app.VelocitySwitch.Position = [221 251 20 45];

            % Create SmoothingKernelSizeKnobLabel
            app.SmoothingKernelSizeKnobLabel = uilabel(app.FilterBankPanel);
            app.SmoothingKernelSizeKnobLabel.HorizontalAlignment = 'center';
            app.SmoothingKernelSizeKnobLabel.Position = [143 1 130 22];
            app.SmoothingKernelSizeKnobLabel.Text = 'Smoothing Kernel Size';

            % Create SmoothingKernelSizeKnob
            app.SmoothingKernelSizeKnob = uiknob(app.FilterBankPanel, 'continuous');
            app.SmoothingKernelSizeKnob.Limits = [1 39];
            app.SmoothingKernelSizeKnob.ValueChangedFcn = createCallbackFcn(app, @SmoothingKernelSizeKnobValueChanged, true);
            app.SmoothingKernelSizeKnob.Position = [164 50 89 89];
            app.SmoothingKernelSizeKnob.Value = 1;

            % Create KernelSizeLabel
            app.KernelSizeLabel = uilabel(app.FilterBankPanel);
            app.KernelSizeLabel.Position = [160 166 130 22];
            app.KernelSizeLabel.Text = 'Kernel Size';

            % Create Label2
            app.Label2 = uilabel(app.UIFigure);
            app.Label2.Position = [1096 406 259 22];
            app.Label2.Text = 'Label2';

            % Create Acceleration
            app.Acceleration = uiaxes(app.UIFigure);
            title(app.Acceleration, 'Acceleration')
            xlabel(app.Acceleration, '')
            ylabel(app.Acceleration, '')
            app.Acceleration.Box = 'on';
            app.Acceleration.XTick = [];
            app.Acceleration.YTick = [];
            app.Acceleration.BackgroundColor = [1 1 1];
            app.Acceleration.Position = [690 12 387 218];

            % Create SaccadeDataLabel
            app.SaccadeDataLabel = uilabel(app.UIFigure);
            app.SaccadeDataLabel.FontWeight = 'bold';
            app.SaccadeDataLabel.Position = [1096 467 90 22];
            app.SaccadeDataLabel.Text = 'Saccade Data:';

            % Create ManualSaccadeLabelingGUI
            app.ManualSaccadeLabelingGUI = uilabel(app.UIFigure);
            app.ManualSaccadeLabelingGUI.FontSize = 24;
            app.ManualSaccadeLabelingGUI.FontWeight = 'bold';
            app.ManualSaccadeLabelingGUI.Position = [593 847 354 31];
            app.ManualSaccadeLabelingGUI.Text = {'Manual Saccade Labeling GUI'; ''};

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SaccadeLabeling_v2

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end