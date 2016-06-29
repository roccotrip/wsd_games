function [correct,given,unique,niter,tot_pos,unique_pos,correct_pos,given_pos,X,Y2] = wsdg(nU,W,P,Y,X,Wsim)


pos = {'n','v','a','r'} ;
nPos = length(pos);
WW = W;

% normalize the graph
D12 = diag(sum(WW,2).^-0.5) ;
D12(isinf(D12)) = 0;
S = D12*WW*D12 ;

% initializations
% initializations
p = P ;
niter = 0;
tot_pos = zeros(1,nPos) ;
correct_pos = zeros(1,nPos) ;
given_pos = zeros(1,nPos) ;
unique_pos = zeros(1,nPos) ;
C = 0;
Y2 = {};


% game dynamics

while true,
    q = zeros(length(S),length(p)) ;
    for i=1:nU
        [~,pi,si] = find(p(i,:)) ;                                          % i's classes
        if length(pi) > 1
            nn = find(S(i,:)) ;                                             % i's neighbours
            for n =1:length(nn)
                j = nn(n) ;
                sij = S(i,j) ;                                              % the similarity among i and j
                [~,pj,sj] = find(p(j,:)) ;                                  % j's classes
                sim = Wsim(pj,pi) ;                                         % similarity among i's classes and j's classes
                sim = sij * sim ;
                q(i,pi)=q(i,pi)+(si.*(sj*sim));
            end
        else
            q(i,pi)=1;
        end
    end
    
    %    q=S*p;
    C=max(-min(q(:))+1e-20,C); % In case of negative similarities
    dummy = p.*(q+C);
    dummySum = sum(dummy,2) ;
    pnew = bsxfun(@rdivide, dummy, dummySum);
    
    diff = norm(p(:)-pnew(:));
    p = pnew;
    niter = niter+1 ;
    
    if diff<10^-4 || niter==10^4
        break;
    end
end
%fprintf('\n');
% EVALUATION


correct = 0;
unique = 0 ;
given = 0 ;

for i=1:nU
    ind_pos = find(ismember(pos,X{i,2})) ;
    tot_pos(1,ind_pos) = tot_pos(1,ind_pos)+1 ;
    if size(find(P(i,:))) == 1
        unique = unique+1;
        unique_pos(1,ind_pos) = unique_pos(1,ind_pos)+1 ;
    end
    answ = Y(i,:);
    probmax = max(pnew(i,:));
    maxidx = find(pnew(i,:) == probmax);
    
    
    X{i,3}=0;
    Y2{i,1}=maxidx(1);
    y22 = Y(i,find(Y(i,:)));
    Y2{i,2}= y22 ;
    [~,y23] = find(P(i,:));
    Y2{i,3}= y23 ;
    Y2{i,4}= X{i,1} ;
    Y2{i,5}= X{i,2} ;
    Y2{i,6}= 0 ;
    Y2{i,7} = length(find(P(i,:)));
    if (probmax == 1 || probmax > 1/size(find(P(i,:)),2))
        given = given+1 ;
        given_pos(1,ind_pos)=given_pos(1,ind_pos)+1;
        if ismember(maxidx, answ) == 1
            correct=correct+1;
            correct_pos(1,ind_pos) = correct_pos(1,ind_pos)+1 ;
            X{i,3}=1;
            Y2{i,6}= 1 ;
        end
        %     else
        %         fprintf('not given:%d - %s - %s\n',i, X{i,1}, X{i,2});
    end
    
end



end


