function ReachTask_twobeam(partNum)
%% Memory-guided reaching task (MATLAB R2014 + Psychtoolbox)
% IR beams: finger lift & return
% IR touch panel: touch coordinates only (TouchInput)

%% ---- Housekeeping ----
clc; close all; sca;
savePath = 'C:\Tahsin\Experiments\ReachTask\p1';
rng('shuffle');

if nargin < 1
    partNum = input('Participant Number: ');
end
sessionNum      = input('Session Number: ');
startBlockInput = input('Starting Block Number (1-40): ');
startBlock      = max(1, min(40, round(startBlockInput)));

%% --- Constants ---
nBlocks = 40;
nTrials = 48;

if ~exist(savePath,'dir'); mkdir(savePath); end

%% ---- Timing parameters ----
tTarget   = 0.150;
tMask     = 0.200;

fixMin   = 0.500;
fixRange = 0.200;

delayMin   = 1.000;
delayRange = 0.500;

ITI = 0.200;   % after finger returns

%% ---- Task parameters ----
gridN        = 5;
targetRadius = 25;
fixSize      = 10;
maskN        = 50;

WHITE = [255 255 255];
BLACK = [0 0 0];
RED   = [255 0 0];
GREEN = [0 255 0];

%% ---- Triggers ----
Trig.sessionID = TriggerInit(2);
Trig.trialStart = 194;
Trig.targetOn   = 195;
Trig.maskOn     = 212;
Trig.goCue      = 213;
Trig.moveStart  = 214;   % finger lift
Trig.touch      = 215;   % screen touch
Trig.moveStop   = 216;   % finger return

%% ---- LabJack / IR Beams ----
[~, ljParams] = LJInit();
pinNo.IRBeamA = 0;
pinNo.IRBeamB = 1;

%% ---- Psychtoolbox setup ----
HideCursor;
PsychDefaultSetup(2);
Screen('Preference','SkipSyncTests',1); ListenChar(2);

screenNumber = 0;
[win, winRect] = Screen('OpenWindow', screenNumber, BLACK);
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[xCenter, yCenter] = RectCenter(winRect);
ifi = Screen('GetFlipInterval', win);

locs = make25Locations(gridN, xCenter, yCenter, winRect);

KbName('UnifyKeyNames');
KeyESC   = KbName('ESCAPE');
KeySpace = KbName('SPACE');

%% ---- TOUCH INPUT OPEN ----
% Uses Psychtoolbox TouchInput API (HID touch frames)
useTouchInput = false;
tdev = [];
try
    tdev = TouchInput('Open');
    useTouchInput = true;
catch
    useTouchInput = false;
end

%% ---- BLOCK LOOP ----
for blockNum = startBlock:nBlocks

    %% ---- Block start ----
    drawFixation(win,xCenter,yCenter,WHITE,fixSize);
    DrawFormattedText(win, ...
        sprintf('Block %d of %d\n\nPress SPACE to start\nESC to quit',blockNum,nBlocks), ...
        'center', yCenter+80, WHITE);
    Screen('Flip', win);

    while 1
        [kd,~,kc] = KbCheck;
        if kd && kc(KeySpace), break; end
        if kd && kc(KeyESC), cleanup(win, useTouchInput, tdev); return; end
    end

    %% ---- Data struct ----
    theData.partNum = partNum;
    theData.sessionNum = sessionNum;
    theData.blockNum = blockNum;
    theData.nTrials = nTrials;

    theData.locIdx      = NaN(nTrials,1);
    theData.redPos      = NaN(nTrials,2);
    theData.greenPos    = NaN(nTrials,2);
    theData.goColor     = NaN(nTrials,1);

    theData.RT          = NaN(nTrials,1);
    theData.MT          = NaN(nTrials,1);
    theData.touchXY     = NaN(nTrials,2);
    theData.distTarget  = NaN(nTrials,1);
    theData.distDistr   = NaN(nTrials,1);
    theData.badTrial    = zeros(nTrials,1);

    theData.tTrialStart = NaN(nTrials,1);
    theData.tTargetOn   = NaN(nTrials,1);
    theData.tMaskOn     = NaN(nTrials,1);
    theData.tGoOn       = NaN(nTrials,1);
    theData.tTouch      = NaN(nTrials,1);
    
    %% ---- Enforced locations (24 locations) ----
    allLocs = setdiff(1:25,13);    % exclude fixation
    locList = repmat(allLocs,1,2);
    locList = locList(randperm(nTrials));

    %% ---- TRIAL LOOP ----
    for t = 1:nTrials

        % Fixation
        drawFixation(win,xCenter,yCenter,WHITE,fixSize);
        Screen('Flip', win);
        WaitSecs(fixMin + rand*fixRange);

         % Choose location
        idx = locList(t);
        posA = locs(idx,:);
        posB = [2*xCenter - posA(1), 2*yCenter - posA(2)];

        if rand < 0.5
            redPos = posA; greenPos = posB;
        else
            redPos = posB; greenPos = posA;
        end

        theData.tTrialStart(t) = GetSecs;
        SendTrigger(Trig.sessionID,Trig.trialStart);

        % Target flash
        Screen('FillRect',win,BLACK);
        drawFixation(win,xCenter,yCenter,WHITE,fixSize);
        drawCircle(win,redPos,targetRadius,RED);
        drawCircle(win,greenPos,targetRadius,GREEN);
        vbl = Screen('Flip',win);
        theData.tTargetOn(t) = vbl;
        SendTrigger(Trig.sessionID,Trig.targetOn);
        WaitSecs(tTarget - 0.5*ifi);

        % Mask
        Screen('FillRect',win,BLACK);
        drawFixation(win,xCenter,yCenter,WHITE,fixSize);
        drawMaskCircles(win,winRect,maskN,targetRadius,RED,GREEN);
        vbl = Screen('Flip',win);
        theData.tMaskOn(t) = vbl;
        SendTrigger(Trig.sessionID,Trig.maskOn);
        WaitSecs(tMask - 0.5*ifi);

        % Delay
        Screen('FillRect',win,BLACK);
        drawFixation(win,xCenter,yCenter,WHITE,fixSize);
        Screen('Flip',win);
        WaitSecs(delayMin + rand*delayRange);

        % GO cue
        if rand < 0.5
            goCol = RED;   goalPos = redPos;   theData.goColor(t) = 1;
        else
            goCol = GREEN; goalPos = greenPos; theData.goColor(t) = 2;
        end

        Screen('FillRect',win,BLACK);
        drawFixation(win,xCenter,yCenter,goCol,fixSize);
        vbl = Screen('Flip',win);
        theData.tGoOn(t) = vbl;
        SendTrigger(Trig.sessionID,Trig.goCue);

        % ----- IR BEAM: finger lift -----
        moveStartTime = LJWaitforRelease_TwoBeam(ljParams,pinNo.IRBeamA,pinNo.IRBeamB);
        SendTrigger(Trig.sessionID,Trig.moveStart);
        theData.RT(t) = moveStartTime - vbl;

        % ----- TOUCH PANEL: first touch -----
        touched = false;
        touchTimeout = 5.0; % seconds (safety so you never get stuck)
        tStart = GetSecs;

        if useTouchInput
            % Clear any old queued events
            while TouchEventAvail(tdev)
                TouchEventGet(tdev);
            end

            while ~touched && (GetSecs - tStart < touchTimeout)

                % Allow quit anytime
                [kd,~,kc] = KbCheck;
                if kd && kc(KeyESC)
                    saveBlock(savePath, partNum, sessionNum, blockNum, theData);
                    cleanup(win, useTouchInput, tdev);
                    return;
                end

                if TouchEventAvail(tdev)
                    evt = TouchEventGet(tdev);

                    % We only take the first BEGIN as the "touch-down"
                    if isfield(evt,'Phase') && strcmpi(evt.Phase,'Begin')
                        touched = true;
                        theData.tTouch(t) = evt.Time;

                        % evt.X / evt.Y are in window pixel coordinates
                        theData.touchXY(t,:) = [evt.X evt.Y];

                        SendTrigger(Trig.sessionID,Trig.touch);
                    end
                end

                WaitSecs(0.001);
            end
        end

        % Fallback to GetMouse if TouchInput didn't fire
        if ~touched
            while ~touched && (GetSecs - tStart < touchTimeout)
                [kd,~,kc] = KbCheck;
                if kd && kc(KeyESC)
                    saveBlock(savePath, partNum, sessionNum, blockNum, theData);
                    cleanup(win, useTouchInput, tdev);
                    return;
                end

                [mx,my,buttons] = GetMouse(win);
                if any(buttons)
                    touched = true;
                    theData.tTouch(t) = GetSecs;
                    theData.touchXY(t,:) = [mx my];
                    SendTrigger(Trig.sessionID,Trig.touch);
                end
                WaitSecs(0.001);
            end
        end

        % If they never touched, mark badTrial and continue 
        if isnan(theData.tTouch(t))
            theData.badTrial(t) = 1;
            theData.MT(t) = NaN;
            theData.distTarget(t) = NaN;
            theData.distDistr(t)  = NaN;
        else
            theData.MT(t) = theData.tTouch(t) - moveStartTime;

            % Distances
            if isequal(goalPos,redPos)
                distractPos = greenPos;
            else
                distractPos = redPos;
            end

            theData.distTarget(t) = norm(theData.touchXY(t,:) - goalPos);
            theData.distDistr(t)  = norm(theData.touchXY(t,:) - distractPos);
        end

        % ----- IR BEAM: finger return -----
        LJWaitforPass_TwoBeam(ljParams,pinNo.IRBeamA, pinNo.IRBeamB);
        SendTrigger(Trig.sessionID,Trig.moveStop);
        WaitSecs(ITI);

        theData.locIdx(t)     = idx;
        theData.redPos(t,:)   = redPos;
        theData.greenPos(t,:) = greenPos;
    end

    % Save block
    fname = sprintf('p%ds%db%d.mat',partNum,sessionNum,blockNum);
    save(fullfile(savePath,fname),'theData');
    fprintf('Saved block %d\n',blockNum);
end

cleanup(win, useTouchInput, tdev);
end

%% ================= HELPER FUNCTIONS =================

function locs = make25Locations(gridN, xC, yC, winRect) % center is excluded
[winW, winH] = RectSize(winRect);
maxOffsetX = winW * 0.30;
maxOffsetY = winH * 0.30;
xs = linspace(-maxOffsetX, maxOffsetX, gridN);
ys = linspace(-maxOffsetY, maxOffsetY, gridN);
locs = zeros(gridN*gridN, 2);
k = 1;
for iy = 1:gridN
    for ix = 1:gridN
        locs(k,:) = [xC + xs(ix), yC + ys(iy)];
        k = k + 1;
    end
end
end

function drawFixation(win, xC, yC, color, fixSize)
Screen('DrawLine', win, color, xC-fixSize, yC, xC+fixSize, yC, 3);
Screen('DrawLine', win, color, xC, yC-fixSize, xC, yC+fixSize, 3);
end

function drawCircle(win, posXY, radius, color)
baseRect = [0 0 radius*2 radius*2];
rect = CenterRectOnPointd(baseRect, posXY(1), posXY(2));
Screen('FillOval', win, color, rect);
end

function drawMaskCircles(win, winRect, n, radius, colA, colB)
[winW, winH] = RectSize(winRect);
margin = radius + 5;
for i = 1:n
    x = margin + rand*(winW - 2*margin);
    y = margin + rand*(winH - 2*margin);
    if rand < 0.5
        c = colA;
    else
        c = colB;
    end
    rect = CenterRectOnPointd([0 0 radius*2 radius*2], x, y);
    Screen('FillOval', win, c, rect);
end
end

function saveBlock(path,p,s,b,data)
fname = sprintf('p%ds%db%d.mat',p,s,b);
save(fullfile(path,fname),'data');
fprintf('Saved %s\n',fname);
end

function cleanup(win, useTouchInput, tdev)
if useTouchInput
    try
        TouchInput('Close', tdev);
    catch
    end
end
Screen('CloseAll');
sca;
ListenChar(0);
ShowCursor;
end
