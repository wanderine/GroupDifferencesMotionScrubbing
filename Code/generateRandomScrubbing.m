function keep = generateRandomScrubbing(numberOfTimepoints,savedTimepoints)

% Generate more realistic motion scrubbing mask, 
% compared to randperm

throw = zeros(numberOfTimepoints,1);

currentTimepoint = randi(10);

while sum(throw) < (numberOfTimepoints - savedTimepoints)
    motionCluster = abs(round(5*randn));
    throw(currentTimepoint:(currentTimepoint + motionCluster-1)) = 1;
    currentTimepoint = currentTimepoint + motionCluster + abs(round(7*randn)) + 3;
end

throw = throw(1:numberOfTimepoints);
keep = 1 - throw;

if sum(keep) > savedTimepoints
    t = numberOfTimepoints;
    while sum(keep) > savedTimepoints
        keep(t) = 0;
        t = t - 1;
    end
elseif sum(keep) < savedTimepoints
    t = numberOfTimepoints;
    while sum(keep) < savedTimepoints
        keep(t) = 1;
        t = t - 1;
    end
end