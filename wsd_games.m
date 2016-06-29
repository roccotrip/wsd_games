close all
clear all

prates = [0,.4] ;
krates = 0:.05:.8;
k=5;
iter = 1; % number of iterations 1

XX={};
YY={};

semSim = {'WsimTFIDF-'} ;
distSim = {'NN','WS'};
datasets = {'S2','S3','S7','S7CG'};                                           % the datasets' directories
pos = {'n','v','a','r'} ;

ldist = length(distSim);
lsem = length(semSim);
lds = length(datasets);

EntireRes.SemMes=semSim;
EntireRes.DistMes=distSim;
EntireRes.Ks=k;
EntireRes.Prates=prates;

basePath = mfilename('fullpath');                                           % get the current file
basePath = fileparts(basePath);                                             % get the current directory

for nds = 1:lds
    Ris = cell(ldist,lsem);
    RisPos = cell(ldist,lsem);
    RisXXX = cell(ldist,lsem);
    RisvalXXX = cell(ldist,lsem);
    dataset = datasets{nds};
    for prate = prates
        for krate = krates
            YY = {'Lemma','Pos','Correct','# senses'};
            pStr = num2str(prate);
            prateName = strcat('prate_',pStr(3:end));
            fprintf('\n%s - k:%d - prate: %.2f - krate: %2f \n',dataset,k,prate,krate);
            for s = 1:lsem
                sem = semSim{s};
                for d = 1:ldist
                    dist = distSim{d};
                    fprintf('\nDataset: %s - sem: %s - dist: %s \n', dataset,sem,dist);
                    base = strcat(basePath(1:end-1),'/data/',dataset,'/') ;
                    if strcmp(dataset, 'S7CG') == 1 && prate > 0
                        Pcsv = dir(fullfile(base,strcat('P-',pStr,'*'))) ;                             % the strategy space files
                    else
                        Pcsv = dir(fullfile(base,'P-u*')) ;                             % the strategy space files
                    end
                    
                    Pcsv = dir(fullfile(base,'P-u*')) ;                             % the strategy space files
                    Wcsv = dir(fullfile(base, strcat(dist,'*'))) ;        % the similarity graph files
                    Ycsv = dir(fullfile(base,'Y-*')) ;                              % the ground truth files
                    Xcsv = dir(fullfile(base,'Words-*')) ;                          % the word-pos files
                    Wsimcsv = dir(fullfile(base, strcat(sem,'*'))) ;                % the sense similarity files
                    
                    
                    l = length(Pcsv) ;                                              % number of texts in the dataset
                    
                    % results initialization
                    XXX={'W','Lemma','Pos'};
                    valXXX=[0,0,0];
                    
                    % results initialization
                    totali = zeros(1,l);
                    totUni = zeros(1,l);
                    totCor = zeros(1,l);
                    totGiv = zeros(1,l);
                    niters = zeros(1,l);
                    totPos = zeros(l,4);
                    uniPos = zeros(l,4);
                    corPos = zeros(l,4);
                    posGiv = zeros(l,4);
                    
                    for x=1:l
                        Pname = Pcsv(x).name ;
                        Wname = Wcsv(x).name ;
                        Yname = Ycsv(x).name ;
                        Xname = Xcsv(x).name ;
                        Wsimname = Wsimcsv(x).name ;

                        P = dlmread(strcat(base,Pname)) ;                           % load the srategy space
                        X = csv2cell(strcat(base,Xname),'fromfile') ;               % load the word-pos
                        W = dlmread(strcat(base,Wname)) ;                           % load the graph
                        Y = dlmread(strcat(base,Yname)) ;                           % load the ground truth
                        Wsim = dlmread(strcat(base,Wsimname)) ;                        % load the sense similarity matrix
                        
                        k=5;
                        nU = length(W);
                        W(W<krate)=0;
                        if prate > 0 && strcmp(dataset, 'S7CG') == 0
                            for p=1:nU
                                [~,pi,si] = find(P(p,:)) ;
                                lpi = length(pi) ;
                                elem = 1:lpi;
                                gpi = geopdf(elem,prate);
                                gpi = gpi/sum(gpi);
                                P(p,pi) = gpi ;
                            end
                        end
                        
                        [correct,given,unique,niter,tot_pos,unique_pos,correct_pos,given_pos,X2o,Y2o]=wsdg(nU,W,P,Y,X,Wsim ) ;
                        YY = [YY;Y2o(:,4:7)];
                        X2o(:,5)={1};
                        for r = 1:nU
                            row = X2o(r,:);
                            lem = row{1};
                            p = row{2};
                            c = row{3};
                            if isempty(c)
                                c=0;
                            end
                            g = row{4};
                            if isempty(g)
                                g=0;
                            end
                            o = row{5};
                            w = strcat(lem,'-', p);
                            idx = find(strcmp(XXX(:,1), w));
                            if isempty(idx)
                                XXX=[XXX;{w,lem,p}];
                                valXXX = [valXXX;c,g,o];
                            else
                                valXXX(idx,:) = valXXX(idx,:)+[c,g,o];
                            end
                        end
                        nW = length(find(Y(:,1)));
                        precision = correct/given ;
                        recall = correct/nU ;
                        F1 = 2*((precision*recall)/(precision+recall));
                        minWsim = min(min(Wsim));
                        maxWsim = max(max(Wsim));
                        meanWsim = mean2(Wsim);
                        stdWsim = std2(Wsim);
                        nonzeroWsim = length(find(Wsim));
                        
                        fprintf('unique:%.3f \t given:%d \t space:%d \t pre: %.3f\t rec: %.3f\t F1: %.3f\tniter: %d\n',unique/nU,given,length(P),precision,recall,F1, niter);
                        totali(iter,x) = nW;
                        totUni(iter,x) = unique ;
                        totCor(iter,x) = correct ;
                        totGiv(iter,x) = given ;
                        niters(iter,x) = niter ;
                        totPos(x,:,iter) = tot_pos ;
                        uniPos(x,:,iter) = unique_pos ;
                        corPos(x,:,iter) = correct_pos ;
                        posGiv(x,:,iter) = given_pos ;
                    end % end texts
                    % calculate the precision, recall and F1 for the entire dataset
                    Ac = sum(totCor)/sum(totGiv) ;
                    Re = sum(totCor)/sum(totali) ;
                    FF1 =  2*((Ac*Re)/(Ac+Re)) ;
                    AcPos = sum(corPos)./sum(posGiv) ;
                    RePos = sum(corPos)./sum(totPos) ;
                    F1pos = 2*((AcPos.*RePos)./(AcPos+RePos)) ;
                    fprintf('Total precision = %.3f \n',Ac) ;
                    fprintf('Total recall = %.3f \n',Re) ;
                    fprintf('Total F1 = %.3f \n',FF1) ;
                    for xx = 1:length(pos)
                        fprintf('Total F1 on %s = %.3f \n',pos{xx}, F1pos(xx)) ;
                    end
                    Ris{s,d}= [Ac,Re,FF1];
                    RisPos{s,d} = [AcPos;RePos;F1pos] ;
                    RisXXX{s,d} = XXX ;
                    RisvalXXX{s,d} = valXXX;
                end % end disSim
            end % end semSim
            nkrate = num2str(krate);
            nkrate = strcat('k',nkrate(3:end));
            EntireRes.(dataset).(prateName).(nkrate).Ris=Ris;
            EntireRes.(dataset).(prateName).(nkrate).RisPos=RisPos;
            EntireRes.(dataset).(prateName).(nkrate).valXXX=RisvalXXX;
            ds = cell2dataset(YY);
            rnameY = strcat(base,'results/WNALL-Y-',dist,num2str(prate),dataset,'.csv') ;
            %export(ds,'file',rnameY,'delimiter',',');
        end % krates
    end % prates
    EntireRes.(dataset).XXX=RisXXX;
    %base = strcat(dsPath,ds,'/results/') ;
    rname = strcat(base,'results/WNALL-EntireRes',dist,'-',num2str(k),dataset,'.mat') ;
    %save(rname,'EntireRes')
end % end corpora
