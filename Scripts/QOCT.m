%% Raster-scan acquisition of SPDC-SFG depth measurements

clearvars -except dev voa;
clc; close all;

targetDir = pwd;

% MEMS scan settings
iniFile = 'mtidevice.ini';
comPort = '';
dataScale = 1.0;
xs = linspace(-0.3, 0.3, 40);
ys = linspace(-0.3, 0.3, 40);
laserOnValue = 255;
mirrorSettle_s = 0.20;

% Delay scan settings
delayRange = 626.5:0.5:632;
odlSettle_s = 0.05;

% Optional power measurement
measureNSR = false;
eta = 0.0085;
NSR = NaN;
P1 = NaN;
P2 = NaN;

Nx = numel(xs);
Ny = numel(ys);
nDelay = numel(delayRange);

% SFG detector
pm = visa('ni', 'USB0::0x1313::0x8091::M00975323::0::INSTR');
pm.InputBufferSize = 2^10;
pm.Timeout = 10;
fopen(pm);

% MEMS scanner
mDev = MTIDevice;
if isempty(comPort)
    mDev.ConnectDevice();
else
    mDev.ConnectDevice(comPort);
end
if mDev.GetLastError() ~= MTIError.MTI_SUCCESS
    error('Could not connect to the MEMS scanning mirror.');
end
mDev.ResetDevicePosition();
devParams = mDev.LoadDeviceParams(iniFile);
mDev.SetDeviceParams(devParams);
mDev.SetDeviceParam(MTIParam.DataScale, dataScale);
mDev.SetDeviceParam(MTIParam.SampleRate, round(Nx * Ny / 0.1));
mDev.SetDeviceParam(MTIParam.MEMSDriverEnable, 1);
pause(0.1);

% Optical delay line
mdl = MDL002('COM6', 'Verbose', true);
mdl.setUnitsPS();

% Power meter for return-power measurements
PM100D_1550 = visa('ni', 'USB0::0x1313::0x8078::P0010849::INSTR');
fopen(PM100D_1550);

mDev.GoToDevicePosition(xs(1), ys(1), laserOnValue, 0);
pause(0.2);
mdl.setAbsoluteDelay(delayRange(1));
pause(0.3);

if measureNSR
    if ~exist('dev', 'var') || ~exist('voa', 'var')
        error('NSR measurement requires initialized laser (dev) and attenuator (voa).');
    end

    writeline(dev, 'SETLDCUR 1 0');
    writeline(voa, 'VSET1:8.000');
    pause(2);
    fprintf(PM100D_1550, 'SENS:CORR:COLL:ZERO:INIT');
    pause(2);
    writeline(dev, 'SETLDCUR 1 58');
    pause(0.1);

    input('Configure the background-off measurement, then press Enter.');
    P1 = PMreading(PM100D_1550);
    input('Configure the background-on measurement, then press Enter.');
    P2 = PMreading(PM100D_1550);
    NSR = 10 * log10((P2 - P1) / (eta * P1));
    fprintf('Measured NSR = %.3f dB\n', NSR);
end

% Measurement cube
meas3D = nan(Ny, Nx, nDelay);
Power_Before_SFG = nan(Ny, Nx);

nCompleted = 0;
totalPixels = Nx * Ny;
tStart = tic;

for row = 1:Ny
    for col = 1:Nx
        mDev.GoToDevicePosition(xs(col), ys(row), laserOnValue, 0);
        pause(mirrorSettle_s);

        Power_Before_SFG(row, col) = PMreading(PM100D_1550);
        trace = nan(1, nDelay);

        for k = 1:nDelay
            mdl.setAbsoluteDelay(delayRange(k));
            pause(odlSettle_s);
            fprintf(pm, 'MEAS:SCAL:FREQ?');
            trace(k) = str2double(fscanf(pm));
        end

        meas3D(row, col, :) = trace;
        mdl.setAbsoluteDelay(delayRange(1));
        pause(0.3);

        nCompleted = nCompleted + 1;
        fprintf('Acquired %d / %d pixels\n', nCompleted, totalPixels);
    end
end

elapsedScan_s = toc(tStart);
fileName = fullfile(targetDir, ...
    ['QOCT_DepthSweep_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.mat']);
save(fileName, 'meas3D', 'xs', 'ys', 'delayRange', ...
    'Power_Before_SFG', 'NSR', 'P1', 'P2', 'elapsedScan_s', 'eta', '-v7.3');

fprintf('Finished in %.1f s\nSaved: %s\n', elapsedScan_s, fileName);

mDev.ResetDevicePosition();
mDev.SetDeviceParam(MTIParam.MEMSDriverEnable, 0);
mDev.DisconnectDevice();
mdl.close();
fclose(pm);
delete(pm);
fclose(PM100D_1550);
delete(PM100D_1550);

function y = PMreading(powerMeter)
    fprintf(powerMeter, 'SENS:CORR:WAV %d', 1550);
    fprintf(powerMeter, 'SENS:POW:UNIT W');
    fprintf(powerMeter, 'CONF:POW');
    samples = zeros(1, 10);
    for n = 1:10
        samples(n) = str2double(query(powerMeter, 'READ?'));
    end
    y = mean(samples) * 1e6;
end
