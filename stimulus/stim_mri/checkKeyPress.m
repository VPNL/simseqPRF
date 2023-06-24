function [] = checkKeyPress(triggerKey)
% We use KbCheck in a while loop instead of KbWait, as KbWait is
% reported to be unreliable and would need a device input as well.
%
% We can either wait for any key press (if triggerKey = 'any') or we
% can wait for a specific trigger key (if triggerKey has any other
% value)
%
iwait = 0;
if strcmpi('any', triggerKey)
    while ~KbCheck(-1)
        WaitSecs(0.01);
    end
else  % wait for the specific key set by 'triggerKey'
    while ~iwait
        %             WaitSecs(0.01);
        [keyIsDown, ~, c] = KbCheck();
        if keyIsDown
            str = KbName(find(c));
            % we just compare the first element in str and triggerKey
            % because if KbName and KbCheck when used togegther can
            % return unwanted characters (for example, KbName(KbCheck)
            % return '5%' when only the '5' key is pressed).
            if strcmpi(str(1), triggerKey(1))
                iwait = 1;
            end
        end
    end
end

return