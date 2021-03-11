function [shortenedMatrix] = peakTrim(peaks,IMUdata)
index = 1;
if (index < 2)
    shortenedMatrix = peaks;
else
while index <= length(peaks)
        if peaks(index,2) < 2 && peaks(index,2) > -2
            peaks(index,:) = [];
            index=index-1;
        end
        index=index+1;
end
shortenedMatrix = peaks;
end
return
end

%% With Standard deviations
% function [shortenedMatrix] = peakTrim(peaks,IMUdata)
% index = 1;
% while index <= length(peaks)
%         if peaks(index,2) > (mean(IMUdata) + 4.5*std(IMUdata)) | peaks(index,2) < (mean(IMUdata) - 4.5*std(IMUdata))
%             peaks(index,:) = [];
%             index=index-1;
%         end
%         index=index+1;
% end
% shortenedMatrix = peaks;
% return
% end