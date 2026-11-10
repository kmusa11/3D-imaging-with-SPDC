%% Measure joint spectral correlations using scanned WaveShaper filters and SFG

% WaveShaper and SFG detector
url = "192.168.1.93";
port = 8001;
signal_port = 60;
idler_port = 51;

% Spectral scan settings
center_wavelength = 1559.67;
delta = 8;
step_nm = 0.1;
signal_width_nm = 0.3;
idler_width_nm = 0.3;
phaseCoeff = -50; % Coefficient for f in THz, not angular-frequency GDD
center_choice = "domain";

coincidenceFile = 'SPDC_JSI_CAR.mat';

if ~exist('pm', 'var') || isempty(pm)
    error('Initialize the SFG detector object (pm) before running the wavelength sweep.');
end

signal = center_wavelength - delta:step_nm:center_wavelength + delta;
idler = center_wavelength - delta:step_nm:center_wavelength + delta;
CAR = nan(numel(idler), numel(signal)); % Measured SFG counts

wl_nm = (1567.5:-0.02:1527.4).';
f_thz = nm2thz(wl_nm);
rect_nm = @(x, mu, width) abs(x - mu) <= width/2;

tStart = tic;
for i = 1:numel(idler)
    for s = 1:numel(signal)
        sig_mask = rect_nm(wl_nm, signal(s), signal_width_nm);
        idl_mask = rect_nm(wl_nm, idler(i), idler_width_nm);

        profile = [f_thz, repmat([50, 0, 0], numel(f_thz), 1)];
        profile(sig_mask, 2) = 0;
        profile(sig_mask, 4) = signal_port;
        profile(idl_mask, 2) = 0;
        profile(idl_mask, 4) = idler_port;
        profile(:, 3) = compute_phases(f_thz, signal(s), idler(i), ...
            phaseCoeff, center_choice);

        [~, order] = sort(profile(:, 1));
        set_profile(profile(order, :), url, port);
        pause(0.1);

        fprintf(pm, 'MEAS:SCAL:FREQ?');
        CAR(i, s) = str2double(fscanf(pm));
        pause(0.2);
    end
    fprintf('Completed idler scan %d / %d\n', i, numel(idler));
end

elapsed_s = toc(tStart);
fileName = ['JSI_SFG_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.mat'];
save(fileName, 'CAR', 'signal', 'idler', 'signal_width_nm', ...
    'idler_width_nm', 'phaseCoeff', 'elapsed_s');
fprintf('Saved %s\n', fileName);

figure('Color', 'w', 'Position', [150 120 1050 760]);
axMain = axes;
imagesc(axMain, signal, idler, CAR);
set(axMain, 'YDir', 'normal', 'TickLabelInterpreter', 'latex', ...
    'FontSize', 22, 'LineWidth', 1.3, 'Box', 'on');
xlabel(axMain, 'Signal wavelength $\lambda_s$ (nm)', 'Interpreter', 'latex');
ylabel(axMain, 'Idler wavelength $\lambda_i$ (nm)', 'Interpreter', 'latex');
cb = colorbar(axMain);
cb.Label.String = 'SFG counts';
cb.Label.Interpreter = 'latex';

if isfile(coincidenceFile)
    C = load(coincidenceFile, 'signal', 'idler', 'coincident_counts');
    if isfield(C, 'coincident_counts') && ...
            isfield(C, 'signal') && isfield(C, 'idler')
        mainPos = axMain.Position;
        insetPos = [mainPos(1) + mainPos(3) * 0.58, ...
            mainPos(2) + mainPos(4) * 0.58, ...
            mainPos(3) * 0.32, mainPos(4) * 0.32];
        axInset = axes('Position', insetPos);
        imagesc(axInset, C.signal, C.idler, C.coincident_counts);
        set(axInset, 'YDir', 'normal', 'FontSize', 12, ...
            'LineWidth', 1.0, 'Box', 'on', 'TickLabelInterpreter', 'latex');
        xlabel(axInset, '$\lambda_s$ (nm)', 'Interpreter', 'latex');
        ylabel(axInset, '$\lambda_i$ (nm)', 'Interpreter', 'latex');
        uistack(axInset, 'top');
    end
end

function thz = nm2thz(nm)
    thz = 299792.45 ./ nm;
end

function check_url(url, port)
    if ~(ischar(url) || isstring(url)), error('url must be char or string'); end
    if ~isscalar(port) || port < 1 || port > 65535, error('port must be 1..65535'); end
end

function result = set_profile(A, url, port)
    check_url(url, port);
    import matlab.net.*
    import matlab.net.http.*
    payload = np2bytes_be_c_row(A);
    body    = MessageBody(uint8(payload));
    hdrs = [ HeaderField('Content-Type','application/octet-stream'), ...
             HeaderField('Accept','*/*') ];
    uri = URI(sprintf('http://%s:%d/', char(url), port));
    req = RequestMessage('post', hdrs, body);
    resp = send(req, uri);
    resTxt = "";
    hdr = resp.getFields('Result');
    if ~isempty(hdr), resTxt = string(hdr(1).Value); end
    if resp.StatusCode ~= matlab.net.http.StatusCode.OK
        error('POST failed: %s | Result: %s', resp.StatusLine, resTxt);
    end
    fprintf("Server Result: %s\n", resTxt);
    result = resTxt;
end

function bytes = np2bytes_be_c_row(A)
    A   = double(A);
    shp = uint32(size(A));
    rnk = uint32(numel(shp));
    h_rank = typecast(swapbytes(rnk), 'uint8');   h_rank = h_rank(:).';
    h_dims = typecast(swapbytes(shp), 'uint8');   h_dims = h_dims(:).';
    header = [h_rank, h_dims];
    if ndims(A) == 2
        Ac = permute(A, [2 1]);
    else
        Ac = permute(A, ndims(A):-1:1);
    end
    payload = typecast(Ac(:), 'uint8'); payload = payload(:).';
    bytes = [header, payload];
end

function phi = compute_phases(f_thz, signal_mu_nm, idler_mu_nm, phaseCoeff, center_choice)
    switch lower(center_choice)
        case 'signal'
            f0 = nm2thz(signal_mu_nm);
        case 'idler'
            f0 = nm2thz(idler_mu_nm);
        otherwise
            f0 = nm2thz(1550);
    end
    phi = 0.5 * phaseCoeff * (f_thz - f0).^2;
    phi = phi(:);
end
