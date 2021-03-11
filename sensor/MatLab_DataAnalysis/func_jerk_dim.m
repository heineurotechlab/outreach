
function [jerk_dim,jerk_dim_log] = func_jerk_dim(vx,vy,vz)
    speed = sqrt(vx.^2 + vy.^2 + vz.^2);
    Ts = 1/1000; % sample time, seconds
    Tt = length(vx) * Ts; % total time of the movement
    time = 0:Ts:Tt-Ts;
    vpeak = max(speed);
    acc = diff(spline(downsample(time,10),downsample(speed,10),time))./Ts;
    acc = [acc,0];
    jerk = diff(acc)./Ts;
    jerk = [jerk,0];
    %plot(time,speed,time,acc,time,jerk)
    normTerm = (vpeak^2 / Tt^3)^-1;
    jerk_dim = -normTerm * trapz(time,jerk.^2);
    jerk_dim_log = -log(normTerm * trapz(time,jerk.^2));
end
