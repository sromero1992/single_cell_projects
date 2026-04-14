function progressbar(c)
%PROGRESSBAR  Single-line command-window progress bar with smoothed ETA.
%
%   progressbar('label')   – initialise (prints label then holds one line)
%   progressbar(pct)       – update progress 0–100
%   progressbar(100)       – finalise (appends newline)
%
%   ETA method
%     Tracks per-block DELTA times (not cumulative).  A 5-block sliding
%     window median drives the forecast, so one slow block doesn't poison
%     the estimate.  ETA is suppressed for the first MIN_SKIP blocks while
%     the estimate stabilises (shows '...' instead).
%
%   Overwrite method
%     Uses \r to return to start of line, then pads the remainder with
%     spaces to erase leftover characters from any longer previous line.
%     drawnow after each write forces the Command Window to refresh.
%
%   Author: Selim Romero (revised TimeSCape v0.2)

    persistent startTime prevC prevElapsed blockDeltaTimes prevLineLen;

    BAR_LEN   = 30;   % visible character width of the progress bar
    WIN_SIZE  =  5;   % sliding window for median delta-time
    MIN_SKIP  =  2;   % blocks before ETA is shown (suppress noisy early estimates)

    % ── Initialise ────────────────────────────────────────────────────────
    if ischar(c) || isstring(c)
        fprintf('%s\n', c);
        startTime      = tic;
        prevC          = 0;
        prevElapsed    = 0;
        blockDeltaTimes = [];
        prevLineLen    = 0;
        return;
    end

    if ~isnumeric(c)
        error('progressbar: input must be a string (init) or a number (update).');
    end

    pct = max(0, min(100, double(c)));

    % ── Finalise ──────────────────────────────────────────────────────────
    if pct >= 100
        total = toc(startTime);
        bar   = repmat('=', 1, BAR_LEN);
        line  = sprintf('100%% [%s] Done in %s', bar, sec2timestr(total));
        % Pad to overwrite any longer previous line, then end with newline
        pad  = max(0, prevLineLen - length(line));
        fprintf('\r%s%s\n', line, repmat(' ', 1, pad));
        drawnow;
        % Reset persistent state
        startTime       = [];
        blockDeltaTimes = [];
        prevC           = 0;
        prevElapsed     = 0;
        prevLineLen     = 0;
        return;
    end

    % ── Skip trivial changes ──────────────────────────────────────────────
    if pct - prevC < 1; return; end

    % ── Per-block delta time ──────────────────────────────────────────────
    % CRITICAL: store the TIME DELTA for this block, not the total elapsed
    % time. Cumulative storage causes wild early ETA because block 1 looks
    % cheap and the projection compounds the error.
    nowElapsed = toc(startTime);
    dt         = nowElapsed - prevElapsed;     % seconds this block took
    prevElapsed = nowElapsed;
    blockDeltaTimes(end+1) = max(dt, 0);      %#ok<AGROW>

    % ── ETA: sliding-window median ─────────────────────────────────────
    n_done       = length(blockDeltaTimes);
    n_total_est  = max(n_done, round(n_done / (pct / 100)));  % estimated total blocks
    n_left       = max(0, n_total_est - n_done);

    if n_done > MIN_SKIP && n_left > 0
        win      = blockDeltaTimes(max(1, end-WIN_SIZE+1) : end);
        avg_dt   = median(win);
        eta_sec  = avg_dt * n_left;
        eta_str  = sec2timestr(eta_sec);
    elseif n_done <= MIN_SKIP
        eta_str  = '...';
    else
        eta_str  = 'almost done';
    end

    % ── Draw bar ──────────────────────────────────────────────────────────
    filled  = round(BAR_LEN * pct / 100);
    bar     = [repmat('=', 1, filled), repmat('-', 1, BAR_LEN - filled)];
    line    = sprintf('%3.0f%% [%s] ETA: %s', pct, bar, eta_str);

    % Pad to clear leftover characters from a previous longer line
    pad     = max(0, prevLineLen - length(line));
    fprintf('\r%s%s', line, repmat(' ', 1, pad));
    drawnow;                        % force Command Window refresh
    prevLineLen = length(line);
    prevC       = pct;
end


% ── Local helper: seconds → human-readable string ─────────────────────────
function s = sec2timestr(sec)
    sec = max(0, sec);
    h   = floor(sec / 3600);
    m   = floor(mod(sec, 3600) / 60);
    sr  = mod(sec, 60);
    if h > 0
        s = sprintf('%dh%02dm', h, m);
    elseif m > 0
        s = sprintf('%dm%02ds', m, floor(sr));
    else
        s = sprintf('%ds', ceil(sr));
    end
end
