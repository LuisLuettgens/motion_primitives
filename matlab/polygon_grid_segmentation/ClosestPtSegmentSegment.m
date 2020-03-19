function [dist_,s,t,c1,c2]=ClosestPtSegmentSegment(p1, q1, p2, q2)

EPSILON = 10e-6;
d1 = q1 - p1; % Direction vector of segment S1
d2 = q2 - p2; % Direction vector of segment S2
r = p1 - p2;
a = d1*d1'; % Squared length of segment S1, always nonnegative
e = d2*d2'; % Squared length of segment S2, always nonnegative
f = d2*r';
% Check if either or both segments degenerate into points
if a <= EPSILON && e <= EPSILON
    % Both segments degenerate into points
    s = 0.0;
    t = 0.0;
    c1 = p1;
    c2 = p2;
    dist_=(c1-c2)*(c1-c2)';
    return
end
if a <= EPSILON
    % First segment degenerates into a point
    s = 0.0;
    t = f / e; % s = 0 => t = (b*s + f) / e = f / e
    t = Clamp(t, 0.0, 1.0);
else
    c = d1*r';
    if e <= EPSILON
        % Second segment degenerates into a point
        t = 0.0;
        s = Clamp(-c / a, 0.0, 1.0); % t = 0 => s = (b*t - c) / a = -c / a
    else
        % The general nondegenerate case starts here
        b = d1*d2';
        denom = a*e-b*b; % Always nonnegative
        % If segments not parallel, compute closest point on L1 to L2 and
        % clamp to segment S1. Else pick arbitrary s (here 0)
        if denom ~= 0.0
            s = Clamp((b*f - c*e) / denom, 0.0, 1.0);
        else
            s = 0.0;
        end
        % Compute point on L2 closest to S1(s) using
        % t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e
        tnom = b*s + f;
        if tnom < 0.0
            t = 0.0;
            s = Clamp(-c / a, 0.0, 1.0);
        elseif tnom > e
            t = 1.0;
            s = Clamp((b - c) / a, 0.0, 1.0);
        else
            t = tnom / e;
        end
   end
end
c1 = p1 + d1 * s;
c2 = p2 + d2 * t;
dist_=(c1-c2)*(c1-c2)';
    
end

function out=Clamp(n, min_, max_)

out=n;
if n < min_
    out=min_;
end
if n > max_
    out=max_;
end

end