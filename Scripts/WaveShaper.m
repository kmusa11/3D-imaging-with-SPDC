%% Configure SPDC signal/idler spectral windows and dispersion compensation

% WaveShaper connection
url = "192.168.1.93";
port = 8001;
signal_port = 60;
idler_port = 51;

% Spectral windows
wl_nm = (1567.5:-0.02:1527.4).';
f_thz = nm2thz(wl_nm);
signal_mu_nm = 1555.48;
idler_mu_nm = 1563.914;
signal_width_nm = 7.5;
idler_width_nm = 7.5;

% Quadratic spectral phase coefficient for frequency in THz
phaseCoeff = -50; % Coefficient for f in THz, not angular-frequency GDD
center_choice = "domain";

rect_nm = @(x, mu, width) abs(x - mu) <= width/2;
sig_mask = rect_nm(wl_nm, signal_mu_nm, signal_width_nm);
idl_mask = rect_nm(wl_nm, idler_mu_nm, idler_width_nm);

profile = [f_thz, repmat([50, 0, 0], numel(f_thz), 1)];
profile(sig_mask, 2) = 0;
profile(sig_mask, 4) = signal_port;
profile(idl_mask, 2) = 0;
profile(idl_mask, 4) = idler_port;
profile(:, 3) = compute_phases(f_thz, signal_mu_nm, idler_mu_nm, ...
    phaseCoeff, center_choice);

[~, order] = sort(profile(:, 1));
profile = profile(order, :);
set_profile(profile, url, port);

receivedProfile = get_profile(url, port);
[~, order] = sort(receivedProfile(:, 1));
receivedProfile = receivedProfile(order, :);

figure('Color', 'w');
yyaxis left;
plot(thz2nm(receivedProfile(:, 1)), receivedProfile(:, 2), 'LineWidth', 1.5);
ylabel('Attenuation (dB)');
yyaxis right;
plot(thz2nm(receivedProfile(:, 1)), receivedProfile(:, 3), 'LineWidth', 1.5);
ylabel('Spectral phase (rad)');
xlabel('Wavelength (nm)');
grid on;

function thz = nm2thz(nm)
    thz = 299792.45 ./ nm;
end

function nm = thz2nm(thz)
    nm = 299792.45 ./ thz;
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

function A = get_profile(url, port)
    check_url(url, port);
    import matlab.net.*; import matlab.net.http.*;
    uri  = URI(sprintf('http://%s:%d/', char(url), port));
    resp = send(RequestMessage('get'), uri);
    hdr = resp.getFields('Result');
    if ~isempty(hdr), fprintf("Server Result: %s\n", string(hdr(1).Value)); end
    if resp.StatusCode ~= matlab.net.http.StatusCode.OK
        error('GET failed: %s', resp.StatusLine);
    end
    A = bytes2np_be_c_row(uint8(resp.Body.Data));
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

function A = bytes2np_be_c_row(bytes)
    assert(numel(bytes) >= 8, 'Buffer too small');
    rnk = swapbytes(typecast(uint8(bytes(1:4)), 'uint32')); r = double(rnk);
    dims = zeros(1, r, 'uint32');
    offs = 4;
    for i = 1:r
        lo = offs + 4*(i-1) + 1; hi = lo + 3;
        dims(i) = swapbytes(typecast(uint8(bytes(lo:hi)), 'uint32'));
    end
    dimsD = double(dims);
    dataStart = 4*(r + 1) + 1;
    vals = typecast(uint8(bytes(dataStart:end)), 'double');
    if r == 2
        Ac = reshape(vals, fliplr(dimsD));
        A  = ipermute(Ac, [2 1]);
    else
        Ac = reshape(vals, fliplr(dimsD));
        A  = ipermute(Ac, ndims(Ac):-1:1);
    end
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
