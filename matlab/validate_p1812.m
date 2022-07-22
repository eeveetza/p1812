% MATLAB script that is used to verify the implementation of
% Recommendation ITU-R P.1812-6 (as defined in the file tl_p1812.m and the
% functions called therefrom) using a set of test terrain profiles provided by the user.
%
% The script reads all the test profiles from the folder defined by
% the variable <test_profiles>, calculates path profile parameters,
% computes the field strength for the input parameters defined in those files,
% and (in case flag_debug is set to 1) logs all the intermediate results
% in *.csv files (placed in the folder defined by the variable: <test_results>).
% Additionally, the script plots terrain profiles in case flag_plot is set to 1.
%
% Author: Ivica Stevanovic (IS), Federal Office of Communications, Switzerland
% Revision History:
% Date            Revision
% 11FEB22         Aligned with ITU-R P.1812-6
% 28JUL20         Introduced alternative method to compute Lbulls w/o using
%                 terrain profile (Attachment 4 to Annex 1)
% 19MAR20         Modified to align to ITU-R P.1812-5
% 07JUL16         Initial version for ITU-R P.1812-4


%% Input variables

clear all;
close all;
fclose all;

try
    
    % add path to the folder where the functions are defined
    s = pwd;
    if ~exist('read_sg3_measurements.m','file')
        addpath([s '/private/'])
    end
    
    if (isOctave)
        page_screen_output(0);
        page_output_immediately(1);
    end
catch
    error('Folder ./private/ does not appear to be on Octave search path.');
    
end

% path to the folder containing test profiles
test_profiles = './validation_profiles/';

% path to the folder where the resulting log files will be saved
test_results = './validation_results/';

% format of the test profile (measurement) files
fileformat='Fryderyk_csv';

% Clutter code type

ClutterCode = 'GlobCover';

%     ClutterCode='default'; % default clutter code assumes land, rural area with R1 = R2 = 10;
%     ClutterCode='TBD'
%     ClutterCode='OFCOM';
%     ClutterCode='NLCD';
%     ClutterCode='LULC';
%     ClutterCode='GlobCover';
%     ClutterCode='DNR1812';

% set to 1 if the csv log files need to be produced
flag_debug = 1;

% set to 1 if the plots of the height profile are to be shown
flag_plot = 0;

% pathprofile is available (=1), not available (=0)
flag_path = 1;

% set to 1 if Attachment 4 to Annex 1 is to be used for computation of
% the spherical earth diffraction Lbs w/o terrain profile analysis

flag4 = 0;

% set variabilities to zero and location percentage to 50
pL = 50;
sigmaL = 0;

%% begin code
% Collect all the filenames .csv in the folder pathname that contain the profile data
filenames = dir(fullfile(test_profiles, '*.csv')); % filenames(i).name is the filename

% create the output directory
[s,mess,messid] = mkdir(test_results);

if (flag_debug==1)
    fid_all = fopen([test_results 'combined_results.csv'], 'w');
    if (fid_all == -1)
        error('The file combined_results.csv could not be opened');
    end
    fprintf(fid_all,'%s,%s,%s,%s,%s,%s\n','Folder','Filename','Dataset #','Reference','Predicted','Deviation: Predicted-Reference');
end

if (length(filenames) < 1)
    error(['There are no .csv files in the test profile folder ' test_profiles]);
end

for iname = 1 : length(filenames)
    
    filename1 = filenames(iname).name;
    fprintf(1,'***********************************************\n');
    fprintf(1,'Processing file %s%s ...\n', test_profiles, filename1);
    fprintf(1,'***********************************************\n');
    
    % read the file and populate sg3db input data structure
    
    sg3db=read_sg3_measurements([test_profiles filename1],fileformat);
    
    sg3db.debug = flag_debug;
    
    % update the data structure with the Tx Power (kW)
    for kindex=1:sg3db.Ndata
        PERP= sg3db.ERPMaxTotal(kindex);
        HRED= sg3db.HRPred(kindex);
        PkW=10^(PERP/10)*1e-3; %kW
        
        if(isnan(PkW))
            % use complementary information from Basic Transmission Loss and
            % received measured strength to compute the transmitter power +
            % gain
            E=sg3db.MeasuredFieldStrength(kindex);
            PL=sg3db.BasicTransmissionLoss(kindex);
            f=sg3db.frequency(kindex);
            PdBkW=-137.2217+E-20*log10(f)+PL;
            PkW=10^(PdBkW/10);
        end
        
        sg3db.TransmittedPower(kindex)=PkW;
    end
    
    sg3db.ClutterCode=[];
    
    % plot the height profile
    
    x=sg3db.x;
    h_gamsl=sg3db.h_gamsl;
    
    
    %% plot the profile
    if (flag_plot)
        figure;
        ax=axes;
        plot(x,h_gamsl,'LineWidth',2,'Color','k');
        set(ax,'XLim',[min(x) max(x)]);
        hTx=sg3db.hTx;
        hRx=sg3db.hRx;
        %area(ax,x,h_gamsl)
        set(0,'DefaulttextInterpreter','none');
        title(['Tx: ' sg3db.TxSiteName ', Rx: ' sg3db.RxSiteName ', ' sg3db.TxCountry sg3db.MeasurementFileName]);
        set(ax,'XGrid','on','YGrid','on');
        xlabel('distance [km]');
        ylabel('height [m]');
    end
    
    %% plot the position of transmitter/receiver
    
    x=sg3db.x;
    h_gamsl=sg3db.h_gamsl;
    
    hTx=sg3db.hTx;
    hRx=sg3db.hRx;
    
    if(flag_plot)
        hold('on');
        
    end
    for dataset = 1:length(hRx)
        if(~isempty(dataset))
            if (dataset > length(hRx) || dataset <1)
                error('The chosen dataset does not exist.')
            end
            fprintf(1,'Computing the fields for Dataset #%d\n', dataset);
            sg3db.userChoiceInt = dataset;
            hhRx=hRx(dataset);
            % this will be a separate function
            % Transmitter
            if(flag_plot)
                if (sg3db.first_point_transmitter)
                    plot([ x(1) x(1)], [h_gamsl(1),h_gamsl(1)+hTx(1)],'LineWidth',2,'Color','b');
                    plot( x(1), h_gamsl(1)+hTx(1), 'Marker','v','Color','b');
                    plot([ x(end) x(end)], [h_gamsl(end),h_gamsl(end)+hhRx],'LineWidth',2,'Color','r');
                    plot( x(end), h_gamsl(end)+hhRx, 'Marker','v','Color','r');
                else
                    plot([ x(end) x(end)], [h_gamsl(end),h_gamsl(end)+hTx(1)],'LineWidth',2,'Color','b');
                    plot( x(end), h_gamsl(1)+hTx(1), 'Marker','v','Color','b');
                    plot([ x(1) x(1)], [h_gamsl(1),h_gamsl(1)+hhRx],'LineWidth',2,'Color','r');
                    plot( x(1), h_gamsl(1)+hhRx, 'Marker','v','Color','r');
                end
            end
        end
        
        
        
        if(~isempty(sg3db.coveragecode))
            
            
            % fill in the  missing fields in Rx clutter
            i=sg3db.coveragecode(end);
            [RxClutterCode RxP1546Clutter R2external] = clutter(i, ClutterCode);
            i=sg3db.coveragecode(1);
            [TxClutterCode TxP1546Clutter R1external] = clutter(i, ClutterCode);
            
            
            sg3db.RxClutterCodeP1546 = RxP1546Clutter;
            
            if(~isempty(sg3db.h_ground_cover))
                if(~isnan(sg3db.h_ground_cover(end)))
                    if (sg3db.h_ground_cover(end) > 3)
                        sg3db.RxClutterHeight = sg3db.h_ground_cover(end);
                    else
                        sg3db.RxClutterHeight = R2external;
                    end
                else
                    sg3db.RxClutterHeight = R2external;
                end
                
                if(~isnan(sg3db.h_ground_cover(1)))
                    sg3db.TxClutterHeight = sg3db.h_ground_cover(1);
                    if (sg3db.h_ground_cover(1) > 3)
                        sg3db.TxClutterHeight = sg3db.h_ground_cover(1);
                    else
                        sg3db.TxClutterHeight = R1external;
                    end
                else
                    sg3db.TxClutterHeight = R1external;
                end
                
                
            else
                sg3db.RxClutterHeight = R2external;
                sg3db.TxClutterHeight = R1external;
            end
        end
        
        %         % Execute P.1812
        fid_log = -1;
        if (flag_debug)
            
            filename2 = [test_results filename1(1:end-4) '_'  num2str(dataset) '.csv'];
            fid_log = fopen(filename2, 'w');
            if (fid_log == -1)
                error_str = [filename2 ' cannot be opened.'];
                error(error_str);
            end
        end
        
        sg3db.fid_log = fid_log;
        %
        
        sg3db.dct = 500;
        sg3db.dcr = 500;
        
        if sg3db.radio_met_code(1) == 1 % Tx at sea
            sg3db.dct = 0;
        end
        
        if sg3db.radio_met_code(end) == 1 %Rx at sea
            sg3db.dcr = 0;
        end
        try
            
            [sg3db.Lb, sg3db.PredictedFieldStrength] = tl_p1812( ...
                sg3db.frequency(dataset)/1e3, ...  %sg3db.frequency is in MHz, tl_p1812 needs GHz
                sg3db.TimePercent(dataset), ...
                sg3db.x, ...
                sg3db.h_gamsl, ...
                sg3db.h_ground_cover, ...
                sg3db.coveragecode, ...
                sg3db.radio_met_code,...
                sg3db.hTx(dataset), ...
                sg3db.hRx(dataset), ...
                sg3db.polHVC(dataset), ...
                'phi_t', sg3db.TxLat, ...
                'phi_r', sg3db.RxLat, ...
                'lam_t', sg3db.TxLon, ...
                'lam_r', sg3db.RxLon, ...
                'pL', pL, ...
                'sigmaL', sigmaL, ...
                'Ptx', sg3db.TransmittedPower(dataset), ...
                'DN', sg3db.DN, ...
                'N0', sg3db.N0, ...
                'dct', sg3db.dct, ...
                'dcr', sg3db.dcr, ...
                'flag4', flag4, ...
                'debug', flag_debug, ...
                'fid_log', sg3db.fid_log);
            
            %fprintf(1,',,%.8f,%.8f\n',sg3db.PredictedFieldStrength,sg3db.Lb)
        catch message
            disp('Input parameters out of bounds');
            
            sg3db.Lb = NaN;
            sg3db.PredictedFieldStrength = NaN;
            
            rethrow(message);
            
            
        end
        if (flag_debug)
            fclose(fid_log);
            
            % print the deviation of the predicted from the measured value,
            % Measurement folder | Measurement File | Dataset | Measured Field Strength | Predicted Field Strength | Deviation from Measurement
            fprintf(fid_all,'%s,%s,%d,%.8f,%.8f,%.8f\n',sg3db.MeasurementFolder,sg3db.MeasurementFileName,dataset, sg3db.MeasuredFieldStrength(dataset), sg3db.PredictedFieldStrength, sg3db.PredictedFieldStrength - sg3db.MeasuredFieldStrength(dataset));
            
        end
        
    end
end % for all files in ./tests

if (flag_debug)
    fclose(fid_all);
end

fprintf(1,'***********************************************\n');
fprintf(1,'Results are written in folder: %s \n', test_results);
fprintf(1,'***********************************************\n');
