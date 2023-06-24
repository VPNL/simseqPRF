function [status, start_time] = newStartScan()
% This code will trigger the 3T GE scanner at CNI using the E-prime trigger
% cable. --Michael 09/30/2013

try
    s = serial('/dev/tty.usbmodem123451', 'BaudRate', 57600);
    fopen(s);
    fprintf(s, '[t]');
    fclose(s);
catch
    err = 1;
end

if exist('err','var') == 0
    start_time  = GetSecs;
    status = 0;
else
    start_time = GetSecs;
    status = 1;
end

end
