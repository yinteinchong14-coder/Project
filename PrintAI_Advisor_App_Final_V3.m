classdef PrintAI_Advisor_App_Final_V3 < matlab.apps.AppBase
    % =====================================================================
    % PrintAI Advisor — Final V3 (Clean, Full)
    % ---------------------------------------------------------------------
    % - Retrain shows progress dialog.
    % - No auto-prediction; only on button press.
    % - Material 2 condition greys out & disables when "None".
    % - Notes only suggest support when support is NOT enabled.
    % - Log to CSV only; if support not enabled => support cols are "0".
    % - No XLSX / highlight logic anywhere.
    % - Compatible with Train_PrintAI_Master_V3() that overwrites a single
    %   PrintAI_Master_Models.mat file.
    % =====================================================================

    %% ===================== UI =====================
    properties (Access = public)
        UIFigure               matlab.ui.Figure

        % Step labels
        Step0Label             matlab.ui.control.Label
        Step1Label             matlab.ui.control.Label
        Step2Label             matlab.ui.control.Label
        Step3Label             matlab.ui.control.Label
        Step4Label             matlab.ui.control.Label
        Step5Label             matlab.ui.control.Label

        % Left viewer
        ViewerPanel            matlab.ui.container.Panel
        UploadSTLButton        matlab.ui.control.Button
        UIAxes                 matlab.ui.control.UIAxes

        % Top-left retrain
        RetrainModelButton     matlab.ui.control.Button

        % Middle controls
        ControlsPanel          matlab.ui.container.Panel
        ControlsGrid           matlab.ui.container.GridLayout
        Material1Label         matlab.ui.control.Label
        Material1DropDown      matlab.ui.control.DropDown
        Material2Label         matlab.ui.control.Label
        Material2DropDown      matlab.ui.control.DropDown
        Cond1Label             matlab.ui.control.Label
        Cond1DropDown          matlab.ui.control.DropDown
        Cond2Label             matlab.ui.control.Label
        Cond2DropDown          matlab.ui.control.DropDown
        ProfileLabel           matlab.ui.control.Label
        ProfileDropDown        matlab.ui.control.DropDown
        ResolutionLabel        matlab.ui.control.Label
        ResolutionDropDown     matlab.ui.control.DropDown

        % Structure panel
        StructurePanel         matlab.ui.container.Panel
        StructureGrid          matlab.ui.container.GridLayout
        InfillPatternLabel     matlab.ui.control.Label
        InfillPatternDropDown  matlab.ui.control.DropDown
        InfillPercentLabel     matlab.ui.control.Label
        InfillPercentSpinner   matlab.ui.control.Spinner
        WallThickLabel         matlab.ui.control.Label
        WallThickSpinner       matlab.ui.control.Spinner
        TopBotThickLabel       matlab.ui.control.Label
        TopBotThickSpinner     matlab.ui.control.Spinner

        % Action buttons
        PredictButton          matlab.ui.control.Button
        SaveParamsButton       matlab.ui.control.Button
        GenerateSWReportButton matlab.ui.control.Button
        LogCSVButton           matlab.ui.control.Button

        % Right results
        ResultsPanel           matlab.ui.container.Panel
        PrintabilityGauge      matlab.ui.control.SemicircularGauge
        SuccessNumLabel        matlab.ui.control.Label

        % Suggestions panel
        SuggestPanel           matlab.ui.container.Panel
        Nozzle1TempLabel       matlab.ui.control.Label
        Nozzle2TempLabel       matlab.ui.control.Label
        BedTempLabel           matlab.ui.control.Label
        FanLabel               matlab.ui.control.Label
        SpeedLabel             matlab.ui.control.Label
        AdhesionOutLabel       matlab.ui.control.Label
        SupportNeededLabel     matlab.ui.control.Label
        SupportStructLabel     matlab.ui.control.Label
        SupportPatternLabel    matlab.ui.control.Label
        SupportDensityLabel    matlab.ui.control.Label
        OrientationValueLabel  matlab.ui.control.Label
        PrimeTowerLabel        matlab.ui.control.Label
        NotesLabel             matlab.ui.control.Label
    end

    %% ===================== STATE =====================
    properties (Access = private)
        CurrentSTLPath string
        Models struct
        Metadata struct
        LastTriSurf
        LastInputs struct
        LastPredictions struct
    end

    %% ===================== CALLBACKS =====================
    methods (Access = private)
        function RetrainModelButtonPushed(app, ~)
            d = uiprogressdlg(app.UIFigure, ...
                'Title','Training Models', ...
                'Message','Please wait while models are being trained...', ...
                'Indeterminate','on', ...
                'Cancelable','off');
            drawnow;
            try
                pause(0.2);
                if exist('Train_PrintAI_Master_V3', 'file') == 2
                    Train_PrintAI_Master_V3(); % overwrites single MAT
                else
                    app.internalQuickTrain('print_trainset.csv');
                end
                if isvalid(d); close(d); end
                uialert(app.UIFigure,'Models retrained and saved to PrintAI_Master_Models.mat','Retrain complete');
                app.Models = []; % force reload on next predict
            catch ME
                if isvalid(d); close(d); end
                uialert(app.UIFigure, "Retrain failed: " + ME.message, 'Retrain Error');
            end
        end

        function UploadSTLButtonPushed(app, ~)
            [f,p] = uigetfile({'*.stl','STL Files (*.stl)'}, 'Select STL file');
            if isequal(f,0), return; end
            fp = fullfile(p,f);
            [F,V] = app.safeReadSTL(fp);
            if isempty(F)
                uialert(app.UIFigure,'Unsupported STL format. Try re-saving as ASCII or binary STL.','STL Error');
                return;
            end
            app.CurrentSTLPath = fp;

            cla(app.UIAxes);
            app.drawPrintBed();
            app.LastTriSurf = trisurf(F,V(:,1),V(:,2),V(:,3), ...
                'Parent',app.UIAxes,'FaceColor',[0.8 0.86 1.0],'EdgeColor','none');
            axis(app.UIAxes,'equal'); camlight(app.UIAxes,'headlight'); lighting(app.UIAxes,'gouraud');
            title(app.UIAxes, sprintf('Loaded: %s', f), 'Interpreter','none');
            % No auto predict
        end

        function ProfileDropDownValueChanged(app, ~)
            p = string(app.ProfileDropDown.Value);
            names = app.resolutionNames();
            switch p
                case "Balanced",    items = names;
                case "Visual",      items = names([1 2 3 4]);
                case "Draft",       items = names([4 5]);
                case "Engineering", items = names([2 3 4]);
                otherwise,          items = names(3);
            end
            app.ResolutionDropDown.Items = items;
            if ~ismember(app.ResolutionDropDown.Value, items)
                app.ResolutionDropDown.Value = items{1};
            end
            % No auto predict
        end

        function Material1Changed(app, ~)
            app.updateConditionOptions(1);
        end

        function Material2Changed(app, ~)
            app.updateConditionOptions(2);
        end

        function InputOnlyChanged(~, ~)
            % No auto predict
        end

        function PredictButtonPushed(app, ~)
            try
                app.doPredictAll();
            catch ME
                uialert(app.UIFigure, ME.message, 'Prediction Error');
            end
        end

        function SaveParametersButtonPushed(app, ~)
            if isempty(app.LastInputs) || isempty(app.LastPredictions)
                uialert(app.UIFigure, 'Please run prediction first.', 'No Data');
                return;
            end
            [file, path] = uiputfile('parameters_summary.txt', 'Save Parameters As');
            if isequal(file,0), return; end
            fp = fullfile(path,file);
            fid = fopen(fp,'w');
            if fid < 0
                uialert(app.UIFigure,'Unable to create file.','File Error'); return;
            end

            fprintf(fid, '=== PRINTAI PARAMETERS SUMMARY ===\n\n');
            fprintf(fid, 'Timestamp: %s\n', datestr(now,31));
            if ~isempty(app.CurrentSTLPath)
                [~, fname, ext] = fileparts(app.CurrentSTLPath);
                fprintf(fid, 'Model: %s%s\n\n', fname, ext);
            else
                fprintf(fid, 'Model: (none loaded)\n\n');
            end

            fprintf(fid, '--- User Inputs ---\n');
            flds = fieldnames(app.LastInputs);
            for i = 1:numel(flds)
                val = app.LastInputs.(flds{i});
                if isstruct(val), continue; end
                fprintf(fid, '%-28s : %s\n', flds{i}, string(val));
            end

            if isfield(app.LastInputs,'geometry')
                g = app.LastInputs.geometry;
                fprintf(fid, '\n--- Geometry ---\n');
                fprintf(fid, 'Length (mm)                : %.3f\n', g.L);
                fprintf(fid, 'Depth (mm)                 : %.3f\n', g.D);
                fprintf(fid, 'Height (mm)                : %.3f\n', g.H);
                fprintf(fid, 'Volume (mm^3)              : %.3f\n', g.Vol);
                fprintf(fid, 'Surface Area (mm^2)        : %.3f\n', g.SA);
                fprintf(fid, 'Aspect Ratio               : %.4f\n', g.aspect);
                fprintf(fid, 'Compactness                : %.6f\n', g.compact);
                fprintf(fid, 'Overhang Ratio             : %.4f\n', g.overhang_ratio);
            end

            fprintf(fid, '\n--- AI Suggestions ---\n');
            flds = fieldnames(app.LastPredictions);
            for i = 1:numel(flds)
                val = app.LastPredictions.(flds{i});
                fprintf(fid, '%-28s : %s\n', flds{i}, string(val));
            end

            fclose(fid);
            uialert(app.UIFigure, sprintf('Parameters saved to:\n%s', fp), 'Saved');
        end

        function GenerateSWReportButtonPushed(app, ~)
            if isempty(app.LastInputs) || isempty(app.LastPredictions)
                uialert(app.UIFigure, 'Please run prediction first.', 'No Data');
                return;
            end
            [file, path] = uiputfile('strength_weight_report.txt', 'Save Report As');
            if isequal(file,0), return; end
            fp = fullfile(path,file);
            fid = fopen(fp,'w');
            if fid < 0
                uialert(app.UIFigure,'Unable to create file.','File Error'); return;
            end

            geom   = app.LastInputs.geometry;
            wtot   = app.LastInputs.total_weight_g;
            fsol   = app.LastInputs.solid_fraction;
            s2w    = app.LastInputs.strength_to_weight;
            wall_t = app.LastInputs.wall_thickness_mm;
            m1     = string(app.LastInputs.material_1);
            m2     = string(app.LastInputs.material_2);
            infpat = string(app.LastInputs.infill_pattern);
            infpct = app.LastInputs.infill_density_pct;
            tbt    = app.LastInputs.top_bottom_thickness_mm;

            [d1,~] = app.materialDensities(m1, m2); % g/mm^3
            d1_cm3 = d1 * 1e3;                       % g/cm^3
            density = wtot / max(geom.Vol * 1e-3, 1e-6);

            switch m1
                case "PLA",        sigma_t = 60;
                case "Tough PLA",  sigma_t = 65;
                case "ABS",        sigma_t = 40;
                case "PETG",       sigma_t = 50;
                otherwise,         sigma_t = 55;
            end
            spec_strength = (sigma_t * s2w) / max(density, 1e-6);
            meanDim = mean([geom.L, geom.D, geom.H]);
            eff_index = (sigma_t / max(density,1e-6)) * (wall_t / max(meanDim,1e-6));

            oh = geom.overhang_ratio;
            if oh < 0.25
                oh_level = 'Low';      oh_note  = 'Low overhang risk — stable without additional support.';
            elseif oh < 0.40
                oh_level = 'Moderate'; oh_note  = 'Moderate risk — partial support recommended.';
            else
                oh_level = 'High';     oh_note  = 'High overhang risk — support required for successful print.';
            end

            fprintf(fid, '=== STRENGTH-TO-WEIGHT ENGINEERING REPORT ===\n\n');
            fprintf(fid, 'Date: %s\n', datestr(now,31));
            if ~isempty(app.CurrentSTLPath)
                fprintf(fid, 'Model File: %s\n\n', app.CurrentSTLPath);
            else
                fprintf(fid, 'Model File: (none loaded)\n\n');
            end

            fprintf(fid, '--- GEOMETRY ---\n');
            fprintf(fid, '  Length (mm): %.2f\n', geom.L);
            fprintf(fid, '  Depth  (mm): %.2f\n', geom.D);
            fprintf(fid, '  Height (mm): %.2f\n', geom.H);
            fprintf(fid, '  Volume (mm^3): %.2f\n', geom.Vol);
            fprintf(fid, '  Surface Area (mm^2): %.2f\n', geom.SA);
            fprintf(fid, '  Aspect Ratio: %.3f\n\n', geom.aspect);

            fprintf(fid, '--- MATERIAL & PROCESS INPUTS ---\n');
            fprintf(fid, '  Material 1: %s\n', m1);
            fprintf(fid, '  Material 2: %s\n', m2);
            fprintf(fid, '  Infill Pattern: %s\n', infpat);
            fprintf(fid, '  Infill Density (%%): %g\n', infpct);
            fprintf(fid, '  Wall Thickness (mm): %g\n', wall_t);
            fprintf(fid, '  Top/Bottom Thickness (mm): %g\n\n', tbt);

            fprintf(fid, '--- STRUCTURAL PERFORMANCE ---\n');
            fprintf(fid, '  Total Weight: %.2f g\n', wtot);
            fprintf(fid, '  Solid Fraction: %.3f\n', fsol);
            fprintf(fid, '  Strength-to-Weight Ratio: %.4f (heuristic)\n', s2w);
            fprintf(fid, '  Specific Strength: %.3f MPa·cm^3/g\n', spec_strength);
            fprintf(fid, '  Efficiency Index: %.6f\n\n', eff_index);

            fprintf(fid, '--- OVERHANG ANALYSIS ---\n');
            fprintf(fid, '  Overhang Ratio: %.3f\n', oh);
            fprintf(fid, '  Risk Level: %s\n', oh_level);
            fprintf(fid, '  Note: %s\n\n', oh_note);

            fprintf(fid, 'Overall Summary:\n');
            fprintf(fid, '  %s\n', oh_note);
            fclose(fid);

            uialert(app.UIFigure, sprintf('Detailed S/W report saved to:\n%s', fp), 'Report Saved');
        end

        function LogCSVButtonPushed(app, ~)
            try
                if isempty(app.LastInputs) || isempty(app.LastPredictions)
                    uialert(app.UIFigure,'Run Predict first to log AI suggestions.','No Data');
                    return;
                end
                csvName  = 'print_trainset.csv';

                if isfile(csvName)
                    T = readtable(csvName);
                else
                    T = table();
                end

                geom = app.LastInputs.geometry;
                m1 = string(app.LastInputs.material_1);
                m2 = string(app.LastInputs.material_2);
                infpat = string(app.LastInputs.infill_pattern);
                infpct = app.LastInputs.infill_density_pct;
                wall_t = app.LastInputs.wall_thickness_mm;
                tbt    = app.LastInputs.top_bottom_thickness_mm;

                % Build one-row record R
                R = table();
                R.job_id = height(T)+1;
                R.material_1 = categorical(m1);
                R.material_2 = categorical(m2);
                R.material_compatibility = categorical(app.computeCompatibility(m1,m2));
                R.filament_condition_1 = categorical("good");
                R.filament_condition_2 = categorical("good");
                R.profile = categorical(string(app.ProfileDropDown.Value));
                R.layer_height_mm = app.layerHeightFromName(string(app.ResolutionDropDown.Value));
                R.infill_pattern = categorical(infpat);
                R.infill_density_pct = infpct;
                R.wall_thickness_mm = wall_t;
                R.top_bottom_thickness_mm = tbt;
                R.model_length_mm = geom.L;
                R.model_depth_mm = geom.D;
                R.model_height_mm = geom.H;
                R.surface_area_mm2 = geom.SA;
                R.volume_mm3 = geom.Vol;
                R.compactness = geom.compact;
                R.aspect_ratio = geom.aspect;
                R.overhang_ratio = geom.overhang_ratio;
                R.overhang_area_mm2 = geom.overhang_area;
                R.geometry_type = categorical("box");

                if ~isempty(T)
                    maybeCols = {'orientation_ex_deg','orientation_ey_deg','orientation_ez_deg'};
                    for k=1:numel(maybeCols)
                        if ismember(maybeCols{k}, string(T.Properties.VariableNames))
                            R.(maybeCols{k}) = 0;
                        end
                    end
                end

                % Support
                R.support_needed = logical(app.LastPredictions.support_enabled);
                if ~R.support_needed
                    R.support_structure   = categorical("0");
                    R.support_pattern     = categorical("0");
                    R.support_density_pct = 0;
                else
                    R.support_structure = categorical(lower(string(app.LastPredictions.support_structure)));
                    R.support_pattern   = categorical(lower(string(app.LastPredictions.support_pattern)));
                    if isnumeric(app.LastPredictions.support_density)
                        supDen = app.LastPredictions.support_density;
                    else
                        supDen = 0;
                    end
                    R.support_density_pct = supDen;
                end

                R.prime_tower = logical(string(app.PrimeTowerLabel.Text) ~= "Prime Tower: Off (Single-Material)");
                R.nozzle1_temp_c = app.LastPredictions.nozzle1_temp_c;
                if isfield(app.LastPredictions,'nozzle2_temp_c')
                    R.nozzle2_temp_c = app.LastPredictions.nozzle2_temp_c;
                else
                    R.nozzle2_temp_c = 0;
                end
                R.bed_temp_c = app.LastPredictions.bed_temp_c;
                R.fan_speed_pct = app.LastPredictions.fan_speed_pct;
                R.print_speed_mms = app.LastPredictions.print_speed_mms;
                R.adhesion_type = categorical(lower(string(app.LastPredictions.adhesion_type)));
                R.print_time_hr = NaN;

                R.orientation_class = categorical(string(app.LastPredictions.orientation_class));
                R.weight_mat1_g = app.LastInputs.total_weight_g * 0.9;
                R.weight_mat2_g = app.LastInputs.total_weight_g * 0.1;
                R.total_weight_g = app.LastInputs.total_weight_g;
                R.strength_to_weight = app.LastInputs.strength_to_weight;
                R.pass_fail = categorical(""); % for manual fill later
                if isfield(app.LastPredictions,'success_probability')
                    R.print_success_prob = app.LastPredictions.success_probability/100;
                else
                    R.print_success_prob = NaN;
                end

                % Stitch/schema align
                if ~isempty(T)
                    need = string(T.Properties.VariableNames);
                    have = string(R.Properties.VariableNames);
                    miss = setdiff(need, have);
                    for m = miss, R.(m) = app.defaultValueForColumn(m, R); end
                    R = R(:, need);
                    missT = setdiff(string(R.Properties.VariableNames), string(T.Properties.VariableNames));
                    for m = missT, T.(m) = app.defaultValueForColumn(m, T); end
                    T = [T; R];
                else
                    T = R;
                end

                writetable(T, csvName);
                uialert(app.UIFigure, sprintf('Logged one entry to %s', csvName), "Logged");
            catch ME
                uialert(app.UIFigure, "Log failed: " + ME.message, "Log Error");
            end
        end
    end

    %% ===================== CORE PREDICTION =====================
    methods (Access = private)
        function doPredictAll(app)
            if isempty(app.CurrentSTLPath)
                uialert(app.UIFigure,'Please upload an STL first.','Missing Input'); return;
            end

            if isempty(app.Models)
                if isfile('PrintAI_Master_Models.mat')
                    S = load('PrintAI_Master_Models.mat');
                    if isfield(S,'models')
                        app.Models = S.models;
                        if isfield(S,'metadata'), app.Metadata = S.metadata; end
                    end
                end
            end

            [F,V] = app.safeReadSTL(app.CurrentSTLPath);
            geom = app.extractGeometry(V,F);

            m1 = string(app.Material1DropDown.Value);
            m2 = string(app.Material2DropDown.Value);
            c1 = lower(string(app.Cond1DropDown.Value));
            c2 = lower(string(app.Cond2DropDown.Value));
            prof = string(app.ProfileDropDown.Value);
            layer_h = app.layerHeightFromName(string(app.ResolutionDropDown.Value));
            infpat = string(app.InfillPatternDropDown.Value);
            infpct = app.InfillPercentSpinner.Value;
            wall_t = app.WallThickSpinner.Value;
            topbot_t = app.TopBotThickSpinner.Value;

            fsol = app.solidFraction(infpct,wall_t,topbot_t,[geom.L,geom.D,geom.H],0.4);
            [d1,d2] = app.materialDensities(m1,m2);
            Vsol = geom.Vol * fsol;

            if geom.overhang_ratio > 0.22 || any(m2 == ["PVA Natural","Breakaway"])
                V1 = 0.9*Vsol; V2 = 0.1*Vsol;
            else
                V1 = Vsol; V2 = 0;
            end
            w1 = d1*V1; w2 = d2*V2; wtot = min(800,w1+w2);

            s2w = (0.35*(infpct/100) + 0.25*(wall_t/(3*0.4)) + 0.25*(topbot_t/(4*0.4))) / max(wtot, 1e-6);
            compat = app.computeCompatibility(m1,m2);

            T = table;
            T.material_1 = categorical(m1);
            T.material_2 = categorical(m2);
            T.material_compatibility = categorical(compat);
            T.filament_condition_1 = categorical(c1);
            T.filament_condition_2 = categorical(c2);
            T.profile = categorical(prof);
            T.layer_height_mm = layer_h;
            T.infill_pattern = categorical(infpat);
            T.infill_density_pct = infpct;
            T.wall_thickness_mm = wall_t;
            T.top_bottom_thickness_mm = topbot_t;

            T.nozzle1_temp_c = 0; T.nozzle2_temp_c = 0; T.bed_temp_c = 0;
            T.fan_speed_pct  = 0; T.print_speed_mms = 0;

            T.model_length_mm = geom.L; T.model_depth_mm  = geom.D; T.model_height_mm = geom.H;
            T.surface_area_mm2= geom.SA; T.volume_mm3 = geom.Vol; T.aspect_ratio = geom.aspect;
            T.compactness = geom.compact; T.overhang_ratio  = geom.overhang_ratio; T.overhang_area_mm2 = geom.overhang_area;

            T.support_needed       = false;
            T.support_structure    = categorical("none");
            T.support_pattern      = categorical("none");
            T.support_density_pct  = 0;
            T.adhesion_type        = categorical("none");

            T.weight_mat1_g = w1; T.weight_mat2_g = w2; T.total_weight_g= wtot;
            T.strength_to_weight = s2w; T.geometry_type = categorical("box");
            T.prime_tower = m2 ~= "None";
            if ~ismember("print_success_prob", string(T.Properties.VariableNames))
                T.print_success_prob = 0;
            end
            T.orientation_class = categorical("Upright"); % default

            % Predictions
            nozzle1 = app.predictRegIfPresent('nozzle1_temp_c', T, 210);
            if m2=="None", nozzle2 = NaN; else, nozzle2 = app.predictRegIfPresent('nozzle2_temp_c', T, 220); end
            bedT    = app.predictRegIfPresent('bed_temp_c',     T, 60);
            fan1    = app.predictRegIfPresent('fan_speed_pct',  T, 40);
            speed   = app.predictRegIfPresent('print_speed_mms',T, 60);
            supDen  = app.predictRegIfPresent('support_density_pct', T, 20);

            if isfield(app.Models, "print_success_prob") || isfield(app.Models, "pass_fail")
                passProb = app.predictPassProb(T); % model if present
            else
                passProb = app.heuristicSuccessProb(T); % fallback
            end

            supportAI    = app.predictLabelIfPresent(T, 'support_needed', false);
            orientClass  = app.predictStringIfPresent(T, 'orientation_class', "Upright");
            supStruct    = app.predictStringIfPresent(T, 'support_structure', "normal");
            supPattern   = app.predictStringIfPresent(T, 'support_pattern', "gyroid");
            adhesion     = app.predictStringIfPresent(T, 'adhesion_type', "brim");

            % Force support if using soluble/breakaway and high overhang
            if any(m2 == ["PVA Natural","Breakaway"]) && geom.overhang_ratio > 0.25
                supportAI = true;
                if string(supStruct)=="none",   supStruct  = "normal"; end
                if string(supPattern)=="none",  supPattern = "grid";   end
                if isnan(supDen) || supDen<=0,  supDen = 20;           end
            end

            % Orientation heuristic for large flat face
            if geom.overhang_ratio < 0.25 && geom.aspect < 2.0
                orientClass = "Flat";
                noteOrient = "Orientation automatically set to Flat due to large flat surface.";
            else
                noteOrient = "";
            end

            [supportEnabled, dispStruct, dispPattern, dispDen] = app.sanitizeSupportOutputs(supportAI, supStruct, supPattern, supDen, geom);
            dispAdhesion = app.cap1(adhesion);
            primeTxt = app.ifelse(m2~="None","On (Dual-Material)","Off (Single-Material)");

            pct = max(0,min(100,100*passProb));
            app.PrintabilityGauge.Value = round(pct,1);
            if pct < 50
                clr = [1 0 0];
            elseif pct < 80
                clr = [0.85 0.6 0];
            else
                clr = [0 0.5 0];
            end
            app.SuccessNumLabel.Text = sprintf('Success Probability: %.1f%%', pct);
            app.SuccessNumLabel.FontColor = clr;

            cla(app.UIAxes);
            app.drawPrintBed();
            app.LastTriSurf = trisurf(F,V(:,1),V(:,2),V(:,3), ...
                'Parent',app.UIAxes,'FaceColor',[0.8 0.86 1.0],'EdgeColor','none');
            app.applyOrientationView(app.LastTriSurf, orientClass);

            app.Nozzle1TempLabel.Text   = sprintf('Nozzle 1 Temperature: %d °C', round(nozzle1));
            app.Nozzle2TempLabel.Text   = sprintf('Nozzle 2 Temperature: %s', app.numOrDash(nozzle2,'°C'));
            app.BedTempLabel.Text       = sprintf('Bed Temperature: %d °C', round(bedT));
            app.FanLabel.Text           = sprintf('Fan Speed: %d %%', round(fan1));
            app.SpeedLabel.Text         = sprintf('Print Speed: %d mm/s', round(speed));
            app.AdhesionOutLabel.Text   = "Adhesion: " + dispAdhesion;

            app.SupportNeededLabel.Text   = "Support Enabled: " + app.ifelse(supportEnabled,"Yes","No");
            app.SupportStructLabel.Text   = "Support Structure: " + dispStruct;
            app.SupportPatternLabel.Text  = "Support Pattern: " + dispPattern;
            app.SupportDensityLabel.Text  = "Support Density: " + dispDen;
            app.OrientationValueLabel.Text= "Orientation: " + orientClass;
            app.PrimeTowerLabel.Text      = "Prime Tower: " + primeTxt;

            notes = strings(0,1);
            if geom.Vol > 1.0e6 && infpct < 30
                notes(end+1) = "Large part with low infill; consider ≥ 40% infill or thicker walls/top/bottom.";
            end
            if strcmpi(string(infpat),'Gyroid')
                notes(end+1) = "Triangle/Grid infill may bridge better than Gyroid for heavy overhangs.";
            end
            % Only suggest support when it is NOT enabled
            if ~supportEnabled && geom.overhang_ratio > 0.35
                notes(end+1) = "Geometry indicates protruding features or high overhang ratio — support is recommended.";
            end
            if strlength(noteOrient) > 0
                notes(end+1) = noteOrient;
            end
            if isempty(notes)
                notes = "All parameters shown are AI-predicted from your training set. Tweak inputs and re-predict to see probability shift.";
            else
                notes = strjoin(notes, newline);
            end
            app.NotesLabel.Text = "Notes:" + newline + notes;

            app.LastInputs = struct( ...
                'material_1', m1, ...
                'material_2', m2, ...
                'infill_pattern', infpat, ...
                'infill_density_pct', infpct, ...
                'wall_thickness_mm', wall_t, ...
                'top_bottom_thickness_mm', topbot_t, ...
                'total_weight_g', wtot, ...
                'solid_fraction', fsol, ...
                'strength_to_weight', s2w, ...
                'geometry', geom);

            app.LastPredictions = struct( ...
                'nozzle1_temp_c', round(nozzle1), ...
                'nozzle2_temp_c', round(nozzle2), ...
                'bed_temp_c', round(bedT), ...
                'fan_speed_pct', round(fan1), ...
                'print_speed_mms', round(speed), ...
                'support_enabled', supportEnabled, ...
                'support_structure', dispStruct, ...
                'support_pattern', dispPattern, ...
                'support_density', dispDen, ...
                'adhesion_type', dispAdhesion, ...
                'orientation_class', orientClass, ...
                'success_probability', pct);
        end

        % ---------- AI helpers ----------
        function y = predictRegIfPresent(app, fieldName, Trow, defaultVal)
            y = defaultVal;
            if isfield(app.Models, fieldName)
                mdl = app.Models.(fieldName);
                X = app.ensurePredictorColumns(Trow, mdl);
                X = app.selectPredictorSubset(X, mdl);
                try
                    y = predict(mdl, X);
                    if ~isscalar(y), y = y(1); end
                catch
                    y = defaultVal;
                end
            end
        end

        function p = predictPassProb(app, Trow)
            if isfield(app.Models,'print_success_prob')
                try
                    mdl = app.Models.print_success_prob;
                    X = app.ensurePredictorColumns(Trow, mdl);
                    X = app.selectPredictorSubset(X, mdl);
                    p = predict(mdl, X); p = min(0.999, max(0.001, p(1)));
                    return;
                catch
                end
            end
            if isfield(app.Models,'pass_fail')
                try
                    mdl = app.Models.pass_fail;
                    X = app.ensurePredictorColumns(Trow, mdl);
                    X = app.selectPredictorSubset(X, mdl);
                    [~, scores] = predict(mdl, X);
                    p = app.posteriorForClass(scores, mdl, [], "pass"); p = min(0.999, max(0.001, p(1)));
                    return;
                catch
                end
            end
            p = app.heuristicSuccessProb(Trow);
        end

        function p = heuristicSuccessProb(~, Trow)
            oh   = double(Trow.overhang_ratio);
            asp  = double(Trow.aspect_ratio);
            inf  = double(Trow.infill_density_pct)/100;
            wall = double(Trow.wall_thickness_mm);
            tbt  = double(Trow.top_bottom_thickness_mm);
            lay  = double(Trow.layer_height_mm);
            m2   = string(Trow.material_2);
            c1   = string(Trow.filament_condition_1);
            c2   = string(Trow.filament_condition_2);

            penalty = 0.60*oh + 0.20*max(0, asp-2)/4;
            base    = 0.78 - penalty;

            base = base + 0.18*inf ...
                        + 0.09*min(wall/1.2,1) ...
                        + 0.07*min(tbt/1.2,1) ...
                        - 0.05*min(lay/0.30,1);

            if any(m2 == ["PVA Natural","Breakaway"]) && oh > 0.25
                base = base + 0.04;
            end

            if any(lower(c1)==["aged","moist"]),  base = base - 0.08; end
            if any(lower(c2)==["aged","moist"]),  base = base - 0.04; end

            p = min(0.98, max(0.05, base));
        end

        function tf = predictLabelIfPresent(app, Trow, modelName, defaultTF)
            tf = defaultTF;
            if isfield(app.Models, modelName)
                mdl = app.Models.(modelName);
                X = app.ensurePredictorColumns(Trow, mdl);
                X = app.selectPredictorSubset(X, mdl);
                try
                    y = predict(mdl, X);
                    ys = string(y);
                    tf = any(ys == ["Yes","1","true","True"]);
                catch
                end
            end
        end

        function s = predictStringIfPresent(app, Trow, modelName, defaultStr)
            s = string(defaultStr);
            if isfield(app.Models, modelName)
                mdl = app.Models.(modelName);
                X = app.ensurePredictorColumns(Trow, mdl);
                X = app.selectPredictorSubset(X, mdl);
                try
                    y = predict(mdl, X);
                    s = string(y);
                    if isempty(s), s = string(defaultStr); end
                catch
                    s = string(defaultStr);
                end
            end
        end

        function [enabled, dispStruct, dispPattern, dispDen] = sanitizeSupportOutputs(app, supportAI, supStruct, supPattern, supDen, geom)
            enabled = logical(supportAI);
            normNone = @(s) any(strcmpi(string(s),["none","no","false","0","—","-",""]));

            if ~enabled
                dispStruct = "0"; dispPattern = "0"; dispDen = "0"; return;
            end

            if normNone(supStruct)
                if geom.overhang_ratio > 0.25, supStruct = "tree"; else, supStruct = "normal"; end
            end
            if normNone(supPattern)
                if      geom.overhang_ratio > 0.35, supPattern = "grid";
                elseif  geom.overhang_ratio > 0.25, supPattern = "gyroid";
                else,                               supPattern = "zigzag";
                end
            end

            if isnan(supDen) || supDen<=0, dispDen = "—"; else, dispDen = sprintf('%d %%', round(supDen)); end
            dispStruct  = app.cap1(supStruct);
            dispPattern = app.cap1(supPattern);
        end
    end

    %% ===================== UI CREATION =====================
    methods (Access = private)
        function createComponents(app)
            app.UIFigure = uifigure('Visible','off','Name','PrintAI Advisor (Final)');
            app.UIFigure.Position = [40 40 1320 720];

            % Base layout rectangles
            viewerPos   = [20 120 520 560];
            inputsPos   = [560 430 360 250];
            structPos   = [560 250 360 170];
            predictPos  = [560 170 360 30];
            resultsPos  = [940 580 360 100];
            suggestPos  = [940 95 360 470];

            stepFont = {'FontSize',14,'FontWeight','bold','FontColor',[0.25 0.25 0.25]};

            % STEP 1
            app.Step0Label = uilabel(app.UIFigure, 'Text','STEP 1: UPDATE TO LATEST DATABASE', ...
                stepFont{:}, 'HorizontalAlignment','left', ...
                'Position',[viewerPos(1) viewerPos(2)+viewerPos(4)+62 520 22]);
            app.RetrainModelButton = uibutton(app.UIFigure,'push','Text','Update Model', ...
                'Position',[viewerPos(1) viewerPos(2)+viewerPos(4)+32 520 28], ...
                'ButtonPushedFcn',@(~,~)app.RetrainModelButtonPushed());

            % STEP 2
            app.Step1Label = uilabel(app.UIFigure, 'Text','STEP 2: UPLOAD STL', stepFont{:}, ...
                'HorizontalAlignment','left', ...
                'Position',[viewerPos(1) viewerPos(2)+viewerPos(4)+2 520 22]);

            % Left viewer
            app.ViewerPanel = uipanel(app.UIFigure,'Title','3D Viewer','Position',viewerPos);
            app.UploadSTLButton = uibutton(app.ViewerPanel,'push','Text','Upload STL', ...
                'Position',[200 515 120 30],'ButtonPushedFcn',@(~,~)app.UploadSTLButtonPushed());
            app.UIAxes = uiaxes(app.ViewerPanel,'Position',[20 40 480 460]);
            title(app.UIAxes,'No model loaded');
            app.drawPrintBed();

            % STEP 3
            app.Step2Label = uilabel(app.UIFigure, 'Text','STEP 3: INPUT PARAMETERS', ...
                stepFont{:}, 'HorizontalAlignment','left', ...
                'Position',[inputsPos(1) inputsPos(2)+inputsPos(4)+12 360 22]);
            app.ControlsPanel = uipanel(app.UIFigure,'Title','Inputs','Position',inputsPos);
            app.ControlsGrid = uigridlayout(app.ControlsPanel,[6 2]);
            app.ControlsGrid.RowHeight   = {28,28,28,28,28,28};
            app.ControlsGrid.ColumnWidth = {160,160};
            app.ControlsGrid.Padding     = [10 8 10 8];
            app.ControlsGrid.RowSpacing  = 8; app.ControlsGrid.ColumnSpacing = 12;

            app.Material1Label = uilabel(app.ControlsGrid,'Text','Material 1'); app.Material1Label.Layout.Row=1; app.Material1Label.Layout.Column=1;
            app.Material1DropDown = uidropdown(app.ControlsGrid,'Items',{'PLA','Tough PLA','ABS','PETG'},'Value','PLA',...
                'ValueChangedFcn',@(~,~)app.Material1Changed());
            app.Material1DropDown.Layout.Row=1; app.Material1DropDown.Layout.Column=2;

            app.Material2Label = uilabel(app.ControlsGrid,'Text','Material 2'); app.Material2Label.Layout.Row=2; app.Material2Label.Layout.Column=1;
            app.Material2DropDown = uidropdown(app.ControlsGrid,'Items',{'None','PLA','Tough PLA','ABS','PETG','PVA Natural','Breakaway'}, ...
                'Value','None','ValueChangedFcn',@(~,~)app.Material2Changed());
            app.Material2DropDown.Layout.Row=2; app.Material2DropDown.Layout.Column=2;

            app.Cond1Label = uilabel(app.ControlsGrid,'Text','Condition (M1)'); app.Cond1Label.Layout.Row=3; app.Cond1Label.Layout.Column=1;
            app.Cond1DropDown = uidropdown(app.ControlsGrid,'Items',{'Good','Aged'},'Value','Good',...
                'ValueChangedFcn',@(~,~)app.InputOnlyChanged());
            app.Cond1DropDown.Layout.Row=3; app.Cond1DropDown.Layout.Column=2;

            app.Cond2Label = uilabel(app.ControlsGrid,'Text','Condition (M2)'); app.Cond2Label.Layout.Row=4; app.Cond2Label.Layout.Column=1;
            app.Cond2DropDown = uidropdown(app.ControlsGrid,'Items',{'Good'},'Value','Good','Enable','off', ...
                'ValueChangedFcn',@(~,~)app.InputOnlyChanged());
            app.Cond2DropDown.Layout.Row=4; app.Cond2DropDown.Layout.Column=2;

            app.ProfileLabel = uilabel(app.ControlsGrid,'Text','Profile'); app.ProfileLabel.Layout.Row=5; app.ProfileLabel.Layout.Column=1;
            app.ProfileDropDown = uidropdown(app.ControlsGrid,'Items',{'Balanced','Visual','Draft','Engineering'}, ...
                'Value','Balanced','ValueChangedFcn',@(~,~)app.ProfileDropDownValueChanged());
            app.ProfileDropDown.Layout.Row=5; app.ProfileDropDown.Layout.Column=2;

            app.ResolutionLabel = uilabel(app.ControlsGrid,'Text','Resolution'); app.ResolutionLabel.Layout.Row=6; app.ResolutionLabel.Layout.Column=1;
            app.ResolutionDropDown = uidropdown(app.ControlsGrid,'Items',app.resolutionNames(),'Value','Normal (0.15 mm)',...
                'ValueChangedFcn',@(~,~)app.InputOnlyChanged());
            app.ResolutionDropDown.Layout.Row=6; app.ResolutionDropDown.Layout.Column=2;

            % Structure
            app.StructurePanel = uipanel(app.UIFigure,'Title','Structure','Position',structPos);
            app.StructureGrid = uigridlayout(app.StructurePanel,[4 2]);
            app.StructureGrid.RowHeight = {28,28,28,28};
            app.StructureGrid.ColumnWidth = {160,160};
            app.StructureGrid.Padding = [10 8 10 8];
            app.StructureGrid.RowSpacing = 8; app.StructureGrid.ColumnSpacing = 12;

            app.InfillPatternLabel = uilabel(app.StructureGrid,'Text','Infill Pattern'); app.InfillPatternLabel.Layout.Row=1; app.InfillPatternLabel.Layout.Column=1;
            app.InfillPatternDropDown = uidropdown(app.StructureGrid,'Items',{'Triangle','Zig zag','Gyroid','Cubic','Grid'},'Value','Gyroid',...
                'ValueChangedFcn',@(~,~)app.InputOnlyChanged());
            app.InfillPatternDropDown.Layout.Row=1; app.InfillPatternDropDown.Layout.Column=2;

            app.InfillPercentLabel = uilabel(app.StructureGrid,'Text','Infill Density'); app.InfillPercentLabel.Layout.Row=2; app.InfillPercentLabel.Layout.Column=1;
            app.InfillPercentSpinner = uispinner(app.StructureGrid,'Limits',[10 100],'Step',10,'Value',40,...
                'ValueChangedFcn',@(~,~)app.InputOnlyChanged());
            app.InfillPercentSpinner.Layout.Row=2; app.InfillPercentSpinner.Layout.Column=2;

            app.WallThickLabel = uilabel(app.StructureGrid,'Text','Wall Thickness (mm)'); app.WallThickLabel.Layout.Row=3; app.WallThickLabel.Layout.Column=1;
            app.WallThickSpinner = uispinner(app.StructureGrid,'Limits',[0.4 2.4],'Step',0.4,'Value',0.8,...
                'ValueChangedFcn',@(~,~)app.InputOnlyChanged());
            app.WallThickSpinner.Layout.Row=3; app.WallThickSpinner.Layout.Column=2;

            app.TopBotThickLabel = uilabel(app.StructureGrid,'Text','Top / Bottom Thickness (mm)'); app.TopBotThickLabel.Layout.Row=4; app.TopBotThickLabel.Layout.Column=1;
            app.TopBotThickSpinner = uispinner(app.StructureGrid,'Limits',[0.4 3.2],'Step',0.2,'Value',0.8,...
                'ValueChangedFcn',@(~,~)app.InputOnlyChanged());
            app.TopBotThickSpinner.Layout.Row=4; app.TopBotThickSpinner.Layout.Column=2;

            % STEP 4 Predict
            app.Step3Label = uilabel(app.UIFigure, 'Text','STEP 4: PREDICT (AI)', stepFont{:}, ...
                'HorizontalAlignment','left', ...
                'Position',[predictPos(1) predictPos(2)+predictPos(4)+5 360 22]);
            app.PredictButton = uibutton(app.UIFigure,'push','Text','Predict & Suggest (AI)', ...
                'Position',predictPos,'ButtonPushedFcn',@(~,~)app.PredictButtonPushed());

            % STEP 5
            app.Step4Label = uilabel(app.UIFigure, 'Text','STEP 5: ADJUST INPUTS & RE-PREDICT', ...
                stepFont{:}, 'HorizontalAlignment','left', ...
                'Position',[resultsPos(1) resultsPos(2)+resultsPos(4)+12 360 22]);

            % Right results - AI Predictability
            app.ResultsPanel = uipanel(app.UIFigure,'Title','AI Predictability','Position',resultsPos);
            panelW = 360; panelH = 100;
            gaugeW = 160; gaugeH = 100;
            centerY = panelH/2;
            gaugeY  = centerY + 50;
            gaugeX  = (panelW - gaugeW)/2;
            app.PrintabilityGauge = uigauge(app.ResultsPanel,'semicircular','Position',[gaugeX gaugeY gaugeW gaugeH]);
            app.PrintabilityGauge.Limits = [0 100];
            app.PrintabilityGauge.MajorTicks = 0:20:100;
            app.PrintabilityGauge.ScaleColors = [1 0 0; 1 1 0; 0 0.8 0];
            app.PrintabilityGauge.ScaleColorLimits = [0 50; 50 80; 80 100];

            app.SuccessNumLabel = uilabel(app.ResultsPanel, ...
                'Text','Success Probability: —', ...
                'FontWeight','bold', ...
                'HorizontalAlignment','center', ...
                'Position',[0 10 panelW 20]);

            % Suggested Parameters
            app.SuggestPanel = uipanel(app.UIFigure,'Title','Suggested Parameters & Recommendations','Position',suggestPos);
            y0 = 420; dy = 24;
            app.Nozzle1TempLabel    = uilabel(app.SuggestPanel,'Text','Nozzle 1 Temperature: —','Position',[18 y0 320 22]); y0=y0-dy;
            app.Nozzle2TempLabel    = uilabel(app.SuggestPanel,'Text','Nozzle 2 Temperature: —','Position',[18 y0 320 22]); y0=y0-dy;
            app.BedTempLabel        = uilabel(app.SuggestPanel,'Text','Bed Temperature: —','Position',[18 y0 320 22]); y0=y0-dy;
            app.FanLabel            = uilabel(app.SuggestPanel,'Text','Fan Speed: —','Position',[18 y0 320 22]); y0=y0-dy;
            app.SpeedLabel          = uilabel(app.SuggestPanel,'Text','Print Speed: —','Position',[18 y0 320 22]); y0=y0-dy;
            app.AdhesionOutLabel    = uilabel(app.SuggestPanel,'Text','Adhesion: —','Position',[18 y0 320 22]); y0=y0-dy;
            app.SupportNeededLabel  = uilabel(app.SuggestPanel,'Text','Support Enabled: —','FontWeight','bold','Position',[18 y0 320 22]); y0=y0-dy;
            app.SupportStructLabel  = uilabel(app.SuggestPanel,'Text','Support Structure: —','Position',[18 y0 320 22]); y0=y0-dy;
            app.SupportPatternLabel = uilabel(app.SuggestPanel,'Text','Support Pattern: —','Position',[18 y0 320 22]); y0=y0-dy;
            app.SupportDensityLabel = uilabel(app.SuggestPanel,'Text','Support Density: —','Position',[18 y0 320 22]); y0=y0-dy;
            app.OrientationValueLabel = uilabel(app.SuggestPanel,'Text','Orientation: —','FontWeight','bold','Position',[18 y0 320 22]); y0=y0-dy;
            app.PrimeTowerLabel     = uilabel(app.SuggestPanel,'Text','Prime Tower: —','Position',[18 y0 320 22]); y0=y0-dy-6;
            app.NotesLabel          = uilabel(app.SuggestPanel,'Text','Notes: —','WordWrap','on','Position',[18 10 320 120]);

            % STEP 6 Save & Analyze
            app.Step5Label = uilabel(app.UIFigure, 'Text','STEP 6: SAVE & ANALYZE', stepFont{:}, ...
                'HorizontalAlignment','left', ...
                'Position',[560 140 360 22]);
            app.SaveParamsButton = uibutton(app.UIFigure, 'push', 'Text', 'Save Parameters', ...
                'Position', [560 110 360 30], 'ButtonPushedFcn', @(~,~)app.SaveParametersButtonPushed());
            app.GenerateSWReportButton = uibutton(app.UIFigure, 'push', 'Text', 'Generate Strength/Weight Report', ...
                'Position', [560 80 360 30], 'ButtonPushedFcn', @(~,~)app.GenerateSWReportButtonPushed());
            app.LogCSVButton = uibutton(app.UIFigure,'push','Text','Log Entry to CSV', ...
                'Position', [560 50 360 30], 'ButtonPushedFcn', @(~,~)app.LogCSVButtonPushed());

            app.UIFigure.Visible = 'on';

            % Initialize defaults
            app.ProfileDropDownValueChanged();
            app.updateConditionOptions(1);
            app.updateConditionOptions(2);
        end
    end

    %% ===================== VISUAL HELPERS =====================
    methods (Access = private)
        function drawPrintBed(app)
            ax = app.UIAxes; hold(ax,'on');
            bedX = [-150 150 150 -150]; bedY = [-150 -150 150 150]; bedZ=[0 0 0 0];
            patch(ax, bedX, bedY, bedZ, [0.90 0.98 0.90], 'EdgeColor', [0.7 0.9 0.7]);
            grid(ax,'on'); view(ax,[45 20]); xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z'); hold(ax,'off');
        end

        function applyOrientationView(app, h, oc)
            ax = app.UIAxes;
            ctr = mean([h.XData(:),h.YData(:),h.ZData(:)],1);
            t = hgtransform('Parent', ax); set(h,'Parent', t);
            switch string(oc)
                case "Flat",       R = eye(4);
                case "LowTilt",    R = makehgtform('xrotate', deg2rad(15));
                case "MediumTilt", R = makehgtform('xrotate', deg2rad(45));
                case "HighTilt",   R = makehgtform('xrotate', deg2rad(60));
                case "Steep",      R = makehgtform('xrotate', deg2rad(75));
                case "Upright",    R = makehgtform('xrotate', deg2rad(90));
                otherwise,         R = eye(4);
            end
            t.Matrix = makehgtform('translate', -ctr) * R * makehgtform('translate', ctr);
            axis(ax,'equal'); camlight(ax,'headlight'); lighting(ax,'gouraud');
        end
    end

    %% ===================== LOGIC HELPERS =====================
    methods (Access = private)
        function names = resolutionNames(~)
            names = {'Extra Fine (0.06 mm)','Fine (0.10 mm)','Normal (0.15 mm)','Fast (0.20 mm)','Extra Fast (0.30 mm)'};
        end

        function h = layerHeightFromName(~, name)
            switch string(name)
                case "Extra Fine (0.06 mm)", h = 0.06;
                case "Fine (0.10 mm)",       h = 0.10;
                case "Normal (0.15 mm)",     h = 0.15;
                case "Fast (0.20 mm)",       h = 0.20;
                case "Extra Fast (0.30 mm)", h = 0.30;
                otherwise,                   h = 0.15;
            end
        end

        function out=ifelse(~,cond,a,b); if cond, out=a; else, out=b; end; end

        function s = cap1(~, txt)
            s = string(txt);
            if strlength(s)==0, s=""; return; end
            s = upper(extractBefore(s,2)) + extractAfter(s,1);
        end

        function updateConditionOptions(app, whichM)
            if whichM==1
                mat = string(app.Material1DropDown.Value);
                dd  = app.Cond1DropDown;
            else
                mat = string(app.Material2DropDown.Value);
                dd  = app.Cond2DropDown;
            end
            if whichM==2 && mat=="None"
                dd.Items = {'Good'}; dd.Value = 'Good'; dd.Enable = 'off';
                try, dd.BackgroundColor = [0.95 0.95 0.95]; catch, end
                return;
            else
                dd.Enable = 'on';
                try, dd.BackgroundColor = [1 1 1]; catch, end
            end
            if any(mat == ["PVA Natural","PETG"])
                dd.Items = {'Good','Moist'};
                if ~ismember(dd.Value, dd.Items), dd.Value = 'Good'; end
            elseif any(mat == ["PLA","Tough PLA","ABS","Breakaway"])
                dd.Items = {'Good','Aged'};
                if ~ismember(dd.Value, dd.Items), dd.Value = 'Good'; end
            else
                dd.Items = {'Good'}; dd.Value = 'Good';
            end
        end

        function [F,V]=safeReadSTL(~,fp)
            F=[]; V=[];
            try
                meshData = stlread(fp);
                if isa(meshData,"triangulation")
                    F = meshData.ConnectivityList; V = meshData.Points;
                elseif isstruct(meshData) && isfield(meshData,"Faces")
                    F = meshData.Faces; V = meshData.Vertices;
                end
            catch
                try
                    [V,F] = stlread_ascii_fallback(fp);
                catch
                    F=[]; V=[];
                end
            end
            function [V,F]=stlread_ascii_fallback(file)
                fid=fopen(file,'r'); if fid<0, error('Cannot open STL'); end
                V=[]; F=[]; verts=[];
                tline=fgetl(fid);
                while ischar(tline)
                    if contains(tline,'vertex')
                        nums=sscanf(tline,'%*s %f %f %f');
                        verts(end+1,:)=nums.'; %#ok<AGROW>
                    elseif contains(tline,'endfacet')
                        n=size(V,1);
                        V=[V; verts(end-2:end,:)]; %#ok<AGROW>
                        F=[F; n+1 n+2 n+3]; %#ok<AGROW>
                        verts=[];
                    end
                    tline=fgetl(fid);
                end
                fclose(fid);
                [V,~,ic]=unique(V,'rows'); F=reshape(ic(F),size(F));
            end
        end

        function geom=extractGeometry(app,V,F)
            V = double(V); F = double(F);
            minv=min(V,[],1); maxv=max(V,[],1);
            L=max(maxv(1)-minv(1),1e-12);
            D=max(maxv(2)-minv(2),1e-12);
            H=max(maxv(3)-minv(3),1e-12);

            p1=V(F(:,1),:); p2=V(F(:,2),:); p3=V(F(:,3),:);
            triN = cross(p2-p1,p3-p1,2);
            triArea = 0.5 * sqrt(sum(triN.^2,2));
            SA = sum(triArea);

            Vvol = sum(dot(p1, cross(p2, p3, 2), 2), 'omitnan') / 6;
            Vol  = abs(Vvol);

            dims    = sort([L D H]);
            aspect  = dims(3)/max(dims(1),1e-12);
            compact = Vol / max(SA,1e-12)^(1.5);

            n = app.faceNormals(V,F);
            steep = n(:,3) < cosd(45);
            overhang_ratio = mean(steep);
            overhang_area  = sum(triArea(steep));

            geom=struct('L',L,'D',D,'H',H,'SA',SA,'Vol',Vol,...
                        'aspect',aspect,'compact',compact,...
                        'overhang_ratio',overhang_ratio,'overhang_area',overhang_area);
        end

        function N=faceNormals(~,V,F)
            p1=V(F(:,1),:);p2=V(F(:,2),:);p3=V(F(:,3),:);
            n=cross(p2-p1,p3-p1,2);
            len = sqrt(sum(n.^2,2)); len(len<1e-12) = 1e-12;
            N = n ./ len;
        end

        function compat=computeCompatibility(~,m1,m2)
            if m2=="None"
                compat="Single Material";
            elseif m1=="ABS" && m2=="PVA Natural"
                compat="Incompatible";
            elseif (m1=="PLA" && (m2=="ABS"||m2=="Tough PLA")) || ...
                   (m1=="ABS" && (m2=="PLA"||m2=="Tough PLA"))
                compat="Poor";
            elseif (m1=="PETG" && m2=="PVA Natural") || (m1=="PVA Natural" && m2=="PETG")
                compat="Fair";
            elseif m1==m2
                compat="Excellent";
            else
                compat="Good";
            end
        end

        function val = defaultValueForColumn(~, colName, likeTable)
            lc = lower(string(colName));
            if contains(lc, ["material","profile","pattern","type","geometry","condition","orientation","adhesion","class","compat"])
                val = categorical(repmat("none", height(likeTable), 1));
            elseif contains(lc, ["needed","prime_tower","support_needed"])
                val = false(height(likeTable),1);
            else
                val = zeros(height(likeTable),1);
            end
        end

        function T = ensurePredictorColumns(~, T, model)
            try
                if isprop(model, 'PredictorNames')
                    need = string(model.PredictorNames);
                else
                    need = string(T.Properties.VariableNames);
                end
                have = string(T.Properties.VariableNames);
                logicalVars = ["prime_tower","support_needed"];

                missing = setdiff(need, have, 'stable');
                for m = missing
                    ml = lower(m);
                    if any(strcmpi(m, logicalVars))
                        T.(m) = false(height(T),1);
                    elseif contains(ml, {'material','profile','pattern','type','geometry','condition','orientation','adhesion','class','compat'})
                        T.(m) = categorical(repmat("none",height(T),1));
                    else
                        T.(m) = zeros(height(T),1);
                    end
                end

                for lv = logicalVars
                    if any(strcmpi(need, lv)) && ismember(lv, string(T.Properties.VariableNames))
                        if ~islogical(T.(lv))
                            if iscategorical(T.(lv)) || isstring(T.(lv)) || iscellstr(T.(lv))
                                T.(lv) = ismember(string(T.(lv)), ["true","yes","1","on"]);
                            else
                                T.(lv) = logical(T.(lv));
                            end
                        end
                    end
                end
            catch
            end
        end

        function X = selectPredictorSubset(~, T, model)
            try
                if isprop(model, 'PredictorNames')
                    need = string(model.PredictorNames);
                else
                    X = T; return;
                end
                have = string(T.Properties.VariableNames);
                keep = intersect(need, have, 'stable');

                if isempty(keep)
                    X = table();
                    for n = need
                        if contains(lower(n), {'material','profile','pattern','type','geometry','condition','orientation','adhesion','class','compat'})
                            X.(n) = categorical("none");
                        else
                            X.(n) = 0;
                        end
                    end
                else
                    X = T(:, keep);
                    [~, idx] = ismember(need, keep);
                    X = [X, table()]; % ensure table
                    X = X(:, idx(idx>0));
                    if width(X) ~= numel(need)
                        miss = need(idx==0);
                        for m = miss
                            if contains(lower(m), {'material','profile','pattern','type','geometry','condition','orientation','adhesion','class','compat'})
                                X.(m) = categorical("none");
                            else
                                X.(m) = 0;
                            end
                        end
                        X = X(:, need);
                    end
                end
            catch
                X = T;
            end
        end

        function p = posteriorForClass(~, scores, model, ~, targetClass)
            try
                classes = string(model.ClassNames);
                idx = find(strcmpi(classes,targetClass),1);
                if isempty(idx)
                    [~, mxi] = max(scores,[],2);
                    p = scores(sub2ind(size(scores),(1:numel(mxi))',mxi));
                else
                    p = scores(:,idx);
                end
            catch
                p = 0.5 + 0*scores(:,1);
            end
        end

        function txt = numOrDash(~, val, unit)
            if nargin<3, unit=""; end
            if isnan(val) || (strcmp(unit,'°C') && val==0)
                txt = "—";
            else
                switch unit
                    case '°C', txt = sprintf('%d °C', round(val));
                    case '%',  txt = sprintf('%d %%', round(val));
                    otherwise, txt = sprintf('%g', val);
                end
            end
        end

        function f = solidFraction(~, infillPct, wall_t, topbot_t, dims, nozzle)
            %#ok<INUSD>
            inf = infillPct / 100;
            L = dims(1); D = dims(2); H = dims(3);
            baseVol = L * D * H;
            shellVol = 3*(L*D*topbot_t + L*H*wall_t + D*H*wall_t);
            shellVol = min(shellVol, 0.9 * baseVol);
            infVol = baseVol * inf * 0.95;
            f = min(1, (shellVol + infVol) / baseVol);
            f = max(0.05, f);
        end

        function [d1,d2]=materialDensities(~,m1,m2)
            d1=density(m1); d2=density(m2);
            function d=density(m)
                switch string(m)
                    case "PLA",          d=1.24e-3;  % g/mm^3
                    case "Tough PLA",    d=1.24e-3;
                    case "ABS",          d=1.05e-3;
                    case "PETG",         d=1.27e-3;
                    case "PVA Natural",  d=1.19e-3;
                    case "Breakaway",    d=1.20e-3;
                    case "None",         d=0;
                    otherwise,           d=1.20e-3;
                end
            end
        end

        % ============== Internal quick training (fallback) ==============
        function internalQuickTrain(~, csvName)
            T = readtable(csvName);
            allTargets = ["pass_fail","support_needed","orientation_class", ...
                          "support_structure","support_pattern","adhesion_type", ...
                          "nozzle1_temp_c","nozzle2_temp_c","bed_temp_c", ...
                          "fan_speed_pct","print_speed_mms","support_density_pct", ...
                          "print_success_prob"];
            predictors = setdiff(T.Properties.VariableNames, allTargets);

            % Convert predictors strings to categorical
            for i = 1:numel(predictors)
                c = predictors{i};
                if iscellstr(T.(c)) || isstring(T.(c))
                    T.(c) = categorical(T.(c));
                end
            end

            models = struct();
            for k = 1:numel(allTargets)
                target = allTargets(k);
                if ~ismember(target, string(T.Properties.VariableNames))
                    continue;
                end
                Y = T.(target);
                X = T(:, predictors);
                try
                    if iscellstr(Y) || isstring(Y), Y = categorical(Y); end
                    if isnumeric(Y)
                        M = fitrensemble(X, Y, "Method","Bag","NumLearningCycles",120, ...
                                         "Learners",templateTree("MaxNumSplits",40,"MinLeafSize",4));
                    else
                        M = fitcensemble(X, Y, "Method","Bag","NumLearningCycles",120, ...
                                         "Learners",templateTree("MaxNumSplits",40,"MinLeafSize",4), ...
                                         "CategoricalPredictors","all");
                    end
                    models.(char(target)) = M;
                catch
                end
            end
            metadata.timestamp = datetime("now");
            metadata.features = predictors;
            save("PrintAI_Master_Models.mat","models","metadata");
        end
    end

    %% ===================== BOOT =====================
    methods (Access = public)
        function app = PrintAI_Advisor_App_Final_V3()
            createComponents(app);
        end
    end
end
