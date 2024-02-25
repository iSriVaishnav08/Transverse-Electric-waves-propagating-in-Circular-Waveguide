classdef new < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                   matlab.ui.Figure
        TabGroup                   matlab.ui.container.TabGroup
        DesignTab                  matlab.ui.container.Tab
        H_zH_yButton               matlab.ui.control.Button
        H_xH_yButton               matlab.ui.control.Button
        E_xE_yButton               matlab.ui.control.Button
        MHzLabel_2                 matlab.ui.control.Label
        LCFreqLabel                matlab.ui.control.Label
        LeastCutOffFrequencyLabel  matlab.ui.control.Label
        MHzLabel                   matlab.ui.control.Label
        CFreqLabel                 matlab.ui.control.Label
        CutOffFrquencyLabel        matlab.ui.control.Label
        HzLabel                    matlab.ui.control.Label
        FrequencyEditField         matlab.ui.control.NumericEditField
        FrequencyEditFieldLabel    matlab.ui.control.Label
        n0Label                    matlab.ui.control.Label
        m0Label                    matlab.ui.control.Label
        mLabel                     matlab.ui.control.Label
        rEditField_2               matlab.ui.control.NumericEditField
        rEditField_2Label          matlab.ui.control.Label
        modevaluemEditField        matlab.ui.control.NumericEditField
        modevaluemEditFieldLabel   matlab.ui.control.Label
        modevaluenEditField        matlab.ui.control.NumericEditField
        modevaluenEditFieldLabel   matlab.ui.control.Label
        RadiusEditField            matlab.ui.control.NumericEditField
        RadiusEditFieldLabel       matlab.ui.control.Label
        TransverseelectricfieldsinsideaCircularWaveguideLabel  matlab.ui.control.Label
        InputParametersPanel       matlab.ui.container.Panel
        rEditField                 matlab.ui.control.NumericEditField
        rEditFieldLabel            matlab.ui.control.Label
        EnterButton                matlab.ui.control.Button
        UIAxes                     matlab.ui.control.UIAxes
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: EnterButton
        function EnterButtonPushed(app, event)
            a = app.RadiusEditField.Value;
            m = app.modevaluemEditField.Value;
            n = app.modevaluenEditField.Value;
            epsilon_r = app.rEditField.Value;
            mu_r = app.rEditField_2.Value;
            
            epsilon=8.85*(10^-12)*epsilon_r;
            mu=4*pi*(10^-7)*mu_r;
            
            % Cutoff frequency calculation
            h=BessDerivZerosBisect2(m,n,3)/a;
            f_c=h/(2*(10^6)*pi*sqrt(mu*epsilon));
            l_f_c=BessDerivZerosBisect2(1,1,3)/(2*(10^6)*a*pi*sqrt(mu*epsilon));


            app.CFreqLabel.Text = num2str(f_c);
            app.LCFreqLabel.Text = num2str(l_f_c);
        end

        % Button pushed function: E_xE_yButton
        function E_xE_yButtonPushed(app, event)
            a = app.RadiusEditField.Value;
            m = app.modevaluemEditField.Value;
            n = app.modevaluenEditField.Value;
            epsilon_r = app.rEditField.Value;
            mu_r = app.rEditField_2.Value;
            f = app.FrequencyEditField.Value;
            epsilon = 8.85e-12 * epsilon_r;
            mu = 4 * pi * 1e-7 * mu_r;
            syms x y
            h = BessDerivZerosBisect2(m, n, 3) / a;
            f_c = h / (2e6 * pi * sqrt(mu * epsilon));
            l_f_c = BessDerivZerosBisect2(1, 1, 3) / (2e6 * a * pi * sqrt(mu * epsilon));
            
            % Variables to write electric and magnetic fields
            omega = 2 * pi * f;
            k = omega * sqrt(mu * epsilon);
            beta = k * sqrt(1 - (l_f_c / f)^2);
            r = sqrt(x^2 + y^2);
            c_phi = x / r;
            s_phi = y / r;
            c_mphi = cos(m * acos(c_phi));
            s_mphi = sin(m * asin(s_phi));
            
            % Electric field equations
            E_x = ((omega * mu) / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * c_phi) - (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * s_phi));
            E_y = ((omega * mu) / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * s_phi) + (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * c_phi));
            E_z = 0;
            
            % Magnetic field equations
            H_x = (beta / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * s_phi) + (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * c_phi));
            H_y = -(beta / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * c_phi) - (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * s_phi));
            H_z = besselj(m, h * r) * c_mphi;
            
            % Define a grid for plotting
            [xGrid, yGrid] = meshgrid(linspace(-a, a, 30)); % Adjust the number of points as needed
            
            % Substitute numeric values for symbolic variables
            E_x_numeric = double(subs(E_x, {x, y}, {xGrid(:), yGrid(:)}));
            E_y_numeric = double(subs(E_y, {x, y}, {xGrid(:), yGrid(:)}));
            
            % Reshape to match the grid size
            E_x_numeric = reshape(E_x_numeric, size(xGrid));
            E_y_numeric = reshape(E_y_numeric, size(yGrid));
            
            % Plotting electric field vectors
            quiver(app.UIAxes, xGrid, yGrid, E_x_numeric, E_y_numeric, 'r');
            title(app.UIAxes, 'Electric field (E_x and E_y)');
            xlabel(app.UIAxes, 'X-axis');
            ylabel(app.UIAxes, 'Y-axis');
            grid(app.UIAxes, 'on');

        end

        % Button pushed function: H_xH_yButton
        function H_xH_yButtonPushed(app, event)
            a = app.RadiusEditField.Value;
            m = app.modevaluemEditField.Value;
            n = app.modevaluenEditField.Value;
            epsilon_r = app.rEditField.Value;
            mu_r = app.rEditField_2.Value;
            f = app.FrequencyEditField.Value;
            epsilon = 8.85e-12 * epsilon_r;
            mu = 4 * pi * 1e-7 * mu_r;
            syms x y
            h = BessDerivZerosBisect2(m, n, 3) / a;
            f_c = h / (2e6 * pi * sqrt(mu * epsilon));
            l_f_c = BessDerivZerosBisect2(1, 1, 3) / (2e6 * a * pi * sqrt(mu * epsilon));
            
            % Variables to write electric and magnetic fields
            omega = 2 * pi * f;
            k = omega * sqrt(mu * epsilon);
            beta = k * sqrt(1 - (l_f_c / f)^2);
            r = sqrt(x^2 + y^2);
            c_phi = x / r;
            s_phi = y / r;
            c_mphi = cos(m * acos(c_phi));
            s_mphi = sin(m * asin(s_phi));
            
            % Electric field equations
            E_x = ((omega * mu) / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * c_phi) - (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * s_phi));
            E_y = ((omega * mu) / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * s_phi) + (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * c_phi));
            E_z = 0;
            
            % Magnetic field equations
            H_x = (beta / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * s_phi) + (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * c_phi));
            H_y = -(beta / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * c_phi) - (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * s_phi));
            H_z = besselj(m, h * r) * c_mphi;
            
            % Define a grid for plotting
            [xGrid, yGrid] = meshgrid(linspace(-a, a, 30)); % Adjust the number of points as needed
            
            % Substitute numeric values for symbolic variables
H_x_numeric = double(subs(H_x, {x, y}, {xGrid(:), yGrid(:)}));
H_y_numeric = double(subs(H_y, {x, y}, {xGrid(:), yGrid(:)}));

% Reshape to match the grid size
H_x_numeric = reshape(H_x_numeric, size(xGrid));
H_y_numeric = reshape(H_y_numeric, size(yGrid));

% Plotting magnetic field vectors
quiver(app.UIAxes, xGrid, yGrid, H_x_numeric, H_y_numeric, 'b');  % Use UIAxes2 for magnetic field
title(app.UIAxes, 'Magnetic field (H_x and H_y)');
xlabel(app.UIAxes, 'X-axis');
ylabel(app.UIAxes, 'Y-axis');
grid(app.UIAxes, 'on');
        end

        % Button pushed function: H_zH_yButton
        function H_zH_yButtonPushed(app, event)
            a = app.RadiusEditField.Value;
            m = app.modevaluemEditField.Value;
            n = app.modevaluenEditField.Value;
            epsilon_r = app.rEditField.Value;
            mu_r = app.rEditField_2.Value;
            f = app.FrequencyEditField.Value;
            epsilon = 8.85e-12 * epsilon_r;
            mu = 4 * pi * 1e-7 * mu_r;
            syms x y
            h = BessDerivZerosBisect2(m, n, 3) / a;
            f_c = h / (2e6 * pi * sqrt(mu * epsilon));
            l_f_c = BessDerivZerosBisect2(1, 1, 3) / (2e6 * a * pi * sqrt(mu * epsilon));
            
            % Variables to write electric and magnetic fields
            omega = 2 * pi * f;
            k = omega * sqrt(mu * epsilon);
            beta = k * sqrt(1 - (l_f_c / f)^2);
            r = sqrt(x^2 + y^2);
            c_phi = x / r;
            s_phi = y / r;
            c_mphi = cos(m * acos(c_phi));
            s_mphi = sin(m * asin(s_phi));
            
            % Electric field equations
            E_x = ((omega * mu) / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * c_phi) - (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * s_phi));
            E_y = ((omega * mu) / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * s_phi) + (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * c_phi));
            E_z = 0;
            
            % Magnetic field equations
            H_x = (beta / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * s_phi) + (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * c_phi));
            H_y = -(beta / h) * (((m / (h * r)) * besselj(m, h * r) * s_mphi * c_phi) - (0.5 * (besselj(m - 1, h * r) - besselj(m + 1, h * r)) * c_mphi * s_phi));
            H_z = besselj(m, h * r) * c_mphi;
            
            % Define a grid for plotting
            [xGrid, yGrid] = meshgrid(linspace(-a, a, 30)); % Adjust the number of points as needed

            % Substitute numeric values for symbolic variables
	    H_y_numeric = double(subs(H_y, {x, y}, {xGrid(:), yGrid(:)}));
	    H_z_numeric = double(subs(H_z, {x, y}, {xGrid(:), yGrid(:)}));

	    % Reshape to match the grid size
	    H_y_numeric = reshape(H_y_numeric, size(yGrid));
	    H_z_numeric = reshape(H_z_numeric, size(xGrid));

	    % Plotting magnetic field vectors
	    quiver(app.UIAxes, xGrid, yGrid, H_z_numeric, H_y_numeric, 'b');  % Use UIAxes2 for magnetic field
	    title(app.UIAxes, 'Magnetic field (H_y and H_z)');
	    xlabel(app.UIAxes, 'X-axis');
	    ylabel(app.UIAxes, 'Y-axis');
	    grid(app.UIAxes, 'on');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 858 715];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [2 1 857 715];

            % Create DesignTab
            app.DesignTab = uitab(app.TabGroup);
            app.DesignTab.Title = 'Design';

            % Create UIAxes
            app.UIAxes = uiaxes(app.DesignTab);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [351 165 477 431];

            % Create InputParametersPanel
            app.InputParametersPanel = uipanel(app.DesignTab);
            app.InputParametersPanel.Title = 'Input Parameters';
            app.InputParametersPanel.Position = [45 157 278 439];

            % Create EnterButton
            app.EnterButton = uibutton(app.InputParametersPanel, 'push');
            app.EnterButton.ButtonPushedFcn = createCallbackFcn(app, @EnterButtonPushed, true);
            app.EnterButton.Position = [87 21 100 22];
            app.EnterButton.Text = 'Enter';

            % Create rEditFieldLabel
            app.rEditFieldLabel = uilabel(app.InputParametersPanel);
            app.rEditFieldLabel.HorizontalAlignment = 'right';
            app.rEditFieldLabel.FontSize = 20;
            app.rEditFieldLabel.Position = [66 199 27 26];
            app.rEditFieldLabel.Text = ' εr';

            % Create rEditField
            app.rEditField = uieditfield(app.InputParametersPanel, 'numeric');
            app.rEditField.Position = [108 201 100 22];
            app.rEditField.Value = 1;

            % Create TransverseelectricfieldsinsideaCircularWaveguideLabel
            app.TransverseelectricfieldsinsideaCircularWaveguideLabel = uilabel(app.DesignTab);
            app.TransverseelectricfieldsinsideaCircularWaveguideLabel.FontSize = 20;
            app.TransverseelectricfieldsinsideaCircularWaveguideLabel.Position = [162 633 489 26];
            app.TransverseelectricfieldsinsideaCircularWaveguideLabel.Text = 'Transverse electric fields inside a Circular Waveguide';

            % Create RadiusEditFieldLabel
            app.RadiusEditFieldLabel = uilabel(app.DesignTab);
            app.RadiusEditFieldLabel.HorizontalAlignment = 'right';
            app.RadiusEditFieldLabel.Position = [96 540 43 22];
            app.RadiusEditFieldLabel.Text = 'Radius';

            % Create RadiusEditField
            app.RadiusEditField = uieditfield(app.DesignTab, 'numeric');
            app.RadiusEditField.Position = [154 540 100 22];
            app.RadiusEditField.Value = 1;

            % Create modevaluenEditFieldLabel
            app.modevaluenEditFieldLabel = uilabel(app.DesignTab);
            app.modevaluenEditFieldLabel.HorizontalAlignment = 'right';
            app.modevaluenEditFieldLabel.Position = [59 422 81 22];
            app.modevaluenEditFieldLabel.Text = ' mode value n';

            % Create modevaluenEditField
            app.modevaluenEditField = uieditfield(app.DesignTab, 'numeric');
            app.modevaluenEditField.Position = [155 422 100 22];
            app.modevaluenEditField.Value = 1;

            % Create modevaluemEditFieldLabel
            app.modevaluemEditFieldLabel = uilabel(app.DesignTab);
            app.modevaluemEditFieldLabel.HorizontalAlignment = 'right';
            app.modevaluemEditFieldLabel.Position = [59 482 81 22];
            app.modevaluemEditFieldLabel.Text = 'mode value m';

            % Create modevaluemEditField
            app.modevaluemEditField = uieditfield(app.DesignTab, 'numeric');
            app.modevaluemEditField.Position = [155 482 100 22];
            app.modevaluemEditField.Value = 1;

            % Create rEditField_2Label
            app.rEditField_2Label = uilabel(app.DesignTab);
            app.rEditField_2Label.HorizontalAlignment = 'right';
            app.rEditField_2Label.FontSize = 20;
            app.rEditField_2Label.Position = [115 294 25 26];
            app.rEditField_2Label.Text = 'μr';

            % Create rEditField_2
            app.rEditField_2 = uieditfield(app.DesignTab, 'numeric');
            app.rEditField_2.Position = [155 298 100 22];
            app.rEditField_2.Value = 1;

            % Create mLabel
            app.mLabel = uilabel(app.DesignTab);
            app.mLabel.Position = [263 540 25 22];
            app.mLabel.Text = 'm';

            % Create m0Label
            app.m0Label = uilabel(app.DesignTab);
            app.m0Label.Position = [263 482 36 22];
            app.m0Label.Text = 'm>=0';

            % Create n0Label
            app.n0Label = uilabel(app.DesignTab);
            app.n0Label.Position = [263 422 26 22];
            app.n0Label.Text = 'n>0';

            % Create FrequencyEditFieldLabel
            app.FrequencyEditFieldLabel = uilabel(app.DesignTab);
            app.FrequencyEditFieldLabel.HorizontalAlignment = 'right';
            app.FrequencyEditFieldLabel.Position = [78 238 62 22];
            app.FrequencyEditFieldLabel.Text = 'Frequency';

            % Create FrequencyEditField
            app.FrequencyEditField = uieditfield(app.DesignTab, 'numeric');
            app.FrequencyEditField.Position = [155 238 100 22];
            app.FrequencyEditField.Value = 100000000;

            % Create HzLabel
            app.HzLabel = uilabel(app.DesignTab);
            app.HzLabel.Position = [259 238 25 22];
            app.HzLabel.Text = 'Hz';

            % Create CutOffFrquencyLabel
            app.CutOffFrquencyLabel = uilabel(app.DesignTab);
            app.CutOffFrquencyLabel.Position = [65 95 97 22];
            app.CutOffFrquencyLabel.Text = 'Cut-Off Frquency';

            % Create CFreqLabel
            app.CFreqLabel = uilabel(app.DesignTab);
            app.CFreqLabel.Position = [205 95 39 22];
            app.CFreqLabel.Text = 'CFreq';

            % Create MHzLabel
            app.MHzLabel = uilabel(app.DesignTab);
            app.MHzLabel.Position = [262 95 30 22];
            app.MHzLabel.Text = 'MHz';

            % Create LeastCutOffFrequencyLabel
            app.LeastCutOffFrequencyLabel = uilabel(app.DesignTab);
            app.LeastCutOffFrequencyLabel.Position = [20 48 137 22];
            app.LeastCutOffFrequencyLabel.Text = 'Least Cut-Off Frequency';

            % Create LCFreqLabel
            app.LCFreqLabel = uilabel(app.DesignTab);
            app.LCFreqLabel.Position = [203 48 46 22];
            app.LCFreqLabel.Text = 'LCFreq';

            % Create MHzLabel_2
            app.MHzLabel_2 = uilabel(app.DesignTab);
            app.MHzLabel_2.Position = [263 48 30 22];
            app.MHzLabel_2.Text = 'MHz';

            % Create E_xE_yButton
            app.E_xE_yButton = uibutton(app.DesignTab, 'push');
            app.E_xE_yButton.ButtonPushedFcn = createCallbackFcn(app, @E_xE_yButtonPushed, true);
            app.E_xE_yButton.Position = [400 95 100 22];
            app.E_xE_yButton.Text = 'E_x,E_y';

            % Create H_xH_yButton
            app.H_xH_yButton = uibutton(app.DesignTab, 'push');
            app.H_xH_yButton.ButtonPushedFcn = createCallbackFcn(app, @H_xH_yButtonPushed, true);
            app.H_xH_yButton.Position = [571 95 100 22];
            app.H_xH_yButton.Text = 'H_x,H_y';

            % Create H_zH_yButton
            app.H_zH_yButton = uibutton(app.DesignTab, 'push');
            app.H_zH_yButton.ButtonPushedFcn = createCallbackFcn(app, @H_zH_yButtonPushed, true);
            app.H_zH_yButton.Position = [729 95 100 22];
            app.H_zH_yButton.Text = 'H_z,H_y';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = new

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

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