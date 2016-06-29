function [ Pmfs ] = mfs( P, nU, prate )

%function Pmfs = mfs( P, nU, prate )
%
% Initializes the strategy space according to a geometric distribution
%
%INPUTS
% Pmf = the strategy space og the game initialized with a uniform
% distribution
% nU = number of players
% prate = the parameter of the geometric distribution
%
%OUTPUTS
% Pmfs = the new strategy space

Pmfs = P ;

for p=1:nU
    [~,pi,~] = find(P(p,:)) ;
    lpi = length(pi) ;
    if lpi>1
        
        elem = 1:lpi;
        gpi = geopdf(elem,prate);
        gpi = gpi/sum(gpi);
        Pmfs(p,pi) = gpi ;
    end
end

end

