function time = sec2hms(t)
%Convert seconds into hour-minute-seconds   
    hours = floor(t / 3600);
    t = t - hours * 3600;
    mins = floor(t / 60);
    secs = floor(t - mins * 60);

    time = [num2str(hours), 'h ', num2str(mins), 'mins ', num2str(secs), 'secs '];
