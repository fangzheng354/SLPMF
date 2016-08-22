function [U,V,TraceU,TraceV]=PGSCD(incompleteNetwork, unknownsA, latentDim, coordGradMaxIter, gradMaxIter, gradeps, numberOfThreads)
%Inputs:
%   incompleteNetwork - incomplete matrix
%   unknownsA - list of unknown entries
%   latentDim - latent dimension
%   coordGradMaxIter - coordinate descent iterations
%   gradMaxIter - gradient descent iterations
%   gradeps - gradient descent epsilon
%   numberOfThreads - number of parallel threads
%
%Outputs:
%   U,V - latent factors
%   TraceU - Objective value Traces (U Iterations)
%   TraceV - Objective value Traces (V Iterations)
    N = size(incompleteNetwork,1);
    M = size(incompleteNetwork,2);
    time=tic;
    
    traceObjVal = 0;
    TraceU=[];
    TraceV=[];
    coordinateSelection = zeros(N+M,2);
    degreeSelection = zeros(N,2);
    
    gradTime = 0;
    tmpGradTime = 0;
    updateTime = 0;
    tmpUpdateTime = 0;
    
    %Initialization
    U = rand(N,latentDim);
    V = rand(M,latentDim);
    NeighsU = [];
    NeighsV = [];
    sumLatentU = 0;
    sumLatentV = 0;
    
    %Sum of all latent factors - one time calculation for all updates
    residualGradsU = zeros(N,latentDim);
    maxDegree = 0;
    for i=1:N
        sumLatentU = sumLatentU + U(i,:);
        outNeigh_i = find(incompleteNetwork(i,:) > 0);%Neighbour finding - NOT SCALABLE
        if(isempty(outNeigh_i)==1)
            outNeigh_i = [];
        end
        if(maxDegree < length(outNeigh_i))
            maxDegree = length(outNeigh_i);
        end
        
        for j=1:length(outNeigh_i)%Iterate over neighbours for objective
            traceObjVal = traceObjVal + incompleteNetwork(i,outNeigh_i(j))*log((incompleteNetwork(i,outNeigh_i(j))+1e-5)/(U(i,:)*V(outNeigh_i(j),:)'+1e-5)) - incompleteNetwork(i,outNeigh_i(j));
        end
    end
    NeighsU = zeros(N,maxDegree+1);
    for i=1:N
        outNeigh_i = find(incompleteNetwork(i,:) > 0);%Neighbour finding - NOT SCALABLE
        if(isempty(outNeigh_i)==1)
            outNeigh_i = [];
        end
        diff = size(NeighsU,2) - length(outNeigh_i)-1;
        NeighsU(i,:)=[length(outNeigh_i) outNeigh_i zeros(1,diff)];
    end
            
    maxDegree = 0;
    for i=1:M
        sumLatentV = sumLatentV + V(i,:);
        inNeigh_i = find(incompleteNetwork(:,i) > 0);%Neighbour finding - NOT SCALABLE
        if(isempty(inNeigh_i)==1)
            inNeigh_i = [];
        end
        if(maxDegree < length(inNeigh_i))
            maxDegree = length(inNeigh_i);
        end
        
    end
    NeighsV = zeros(M,maxDegree+1);
    for i=1:M
        inNeigh_i = find(incompleteNetwork(:,i) > 0);%Neighbour finding - NOT SCALABLE
        if(isempty(inNeigh_i)==1)
            inNeigh_i = [];
        end
        diff = size(NeighsV,2) - length(inNeigh_i)-1;
        NeighsV(i,:)=[length(inNeigh_i) inNeigh_i' zeros(1,diff)];
    end
    
    for i=1:N
        traceObjVal = traceObjVal + U(i,:)*sumLatentV';
    end
    
    traceTime2 = 0;
    if (matlabpool('size') == 0)
        matlabpool(numberOfThreads);
    end
    extime = 0;
    iterations = 0;
    for coordCount = 1:coordGradMaxIter%Coordinate Descent
        disp(['coordDes ' num2str(coordCount)]);
        
        itered = 0;%U
        while(true)
            iterations = iterations + 1;
            itered = itered + 1;
%             disp(num2str(itered));
            if(itered > 1000)
                break;
            end

            tmpextime = tic();
            FinalList = [];%Selecting the parallel coordinates
            targetSet = 1:N;
            neighSet = [];
            index=randi(length(targetSet));
            next=targetSet(itered);
            FinalList = [FinalList next];
            neighSet = union(neighSet,NeighsU(next,2:NeighsU(next,1)+1));
            tmpIter = 0;
            while(length(FinalList) < numberOfThreads && tmpIter < numberOfThreads)
                tmpIter = tmpIter + 1;
                index=randi(length(targetSet));
                next=targetSet(index);
                if(isempty(intersect(neighSet,NeighsU(next,2:NeighsU(next,1)+1)))==1)
                    FinalList = [FinalList next];
                    neighSet = union(neighSet,NeighsU(next,2:NeighsU(next,1)+1));
                end
            end
            
            tmptmps = U(FinalList,:);
            News = zeros(length(FinalList),latentDim);
            News = mat2cell(News,ones(1,length(FinalList)),latentDim);
            VNeighs = cell(length(FinalList),1);
            Us = mat2cell(U(FinalList,:),ones(1,length(FinalList)),latentDim);
            for Ui=1:length(FinalList)
                VNeighs{Ui,1} = V(NeighsU(FinalList(Ui),2:NeighsU(FinalList(Ui),1)+1),:);
            end
            parfor Ui=1:length(FinalList)
%                 [News{Ui,1}]=gradDes23(VNeighs{Ui,1},U(Ui,:),sumLatentV,gradMaxIter,gradeps);
                [tmpVec1]=gradDes2(VNeighs{Ui,1},Us{Ui,1},sumLatentV,gradMaxIter,gradeps);
                News{Ui,1} = tmpVec1;
            end
            extime = extime + toc(tmpextime);
            News = cell2mat(News);
            U(FinalList,:) = News;
            sumLatentU = sumLatentU - sum(tmptmps) + sum(News);
            tmpTraceTime = tic();
            descentAmount = 0;
            for Ui=1:length(FinalList)
                outNeigh_i = NeighsU(FinalList(Ui),2:NeighsU(FinalList(Ui),1)+1);
                tmptmp_i = tmptmps(Ui,:);
                for j=1:length(outNeigh_i)
                    descentAmount = descentAmount - incompleteNetwork(FinalList(Ui),outNeigh_i(j))*log((incompleteNetwork(FinalList(Ui),outNeigh_i(j))+1e-5)/(U(FinalList(Ui),:)*V(outNeigh_i(j),:)'+1e-5)) + incompleteNetwork(FinalList(Ui),outNeigh_i(j))*log((incompleteNetwork(FinalList(Ui),outNeigh_i(j))+1e-5)/(tmptmp_i*V(outNeigh_i(j),:)'+1e-5));
                    traceObjVal = traceObjVal + incompleteNetwork(FinalList(Ui),outNeigh_i(j))*log((incompleteNetwork(FinalList(Ui),outNeigh_i(j))+1e-5)/(U(FinalList(Ui),:)*V(outNeigh_i(j),:)'+1e-5)) - incompleteNetwork(FinalList(Ui),outNeigh_i(j))*log((incompleteNetwork(FinalList(Ui),outNeigh_i(j))+1e-5)/(tmptmp_i*V(outNeigh_i(j),:)'+1e-5));
                end
                descentAmount = descentAmount + tmptmp_i*sumLatentV' - U(FinalList(Ui),:)*sumLatentV';
                traceObjVal = traceObjVal - tmptmp_i*sumLatentV' + U(FinalList(Ui),:)*sumLatentV';
            end
            traceTime2 = traceTime2 + toc(tmpTraceTime);

            traceTime = toc(time);
            TraceU=[TraceU [traceTime-traceTime2;iterations;traceObjVal]];
        end
        
        itered = 0;%V
        while(true)
            iterations = iterations + 1;
            itered = itered + 1;
%             disp(num2str(itered));
            if(itered > 1000)
                break;
            end
            
            tmpextime = tic();
            FinalList = [];%Selecting the parallel coordinates
            targetSet = 1:M;
            neighSet = [];
            index=randi(length(targetSet));
            next=targetSet(itered);
            FinalList = [FinalList next];
            neighSet = union(neighSet,NeighsV(next,2:NeighsV(next,1)+1));
            tmpIter = 0;
            while(length(FinalList) < numberOfThreads && tmpIter < numberOfThreads)
                tmpIter = tmpIter + 1;
                index=randi(length(targetSet));
                next=targetSet(index);
                if(isempty(intersect(neighSet,NeighsV(next,2:NeighsV(next,1)+1)))==1)
                    FinalList = [FinalList next];
                    neighSet = union(neighSet,NeighsV(next,2:NeighsV(next,1)+1));
                end
            end
            
            tmptmps = V(FinalList,:);
            News = zeros(length(FinalList),latentDim);
            News = mat2cell(News,ones(1,length(FinalList)),latentDim);
            UNeighs = cell(length(FinalList),1);
            Vs = mat2cell(V(FinalList,:),ones(1,length(FinalList)),latentDim);
            for Vi=1:length(FinalList)
                UNeighs{Vi,1} = U(NeighsV(FinalList(Vi),2:NeighsV(FinalList(Vi),1)+1),:);
            end
            parfor Vi=1:length(FinalList)
                [tmpVec2]=gradDes2(UNeighs{Vi,1},Vs{Vi,1},sumLatentU,gradMaxIter,gradeps);
                News{Vi,1} = tmpVec2;
            end
            extime = extime + toc(tmpextime);
            News = cell2mat(News);
            V(FinalList,:) = News;
            sumLatentV = sumLatentV - sum(tmptmps) + sum(News);
            tmpTraceTime = tic();
            descentAmount = 0;
            for Vi=1:length(FinalList)
                inNeigh_i = NeighsV(FinalList(Vi),2:NeighsV(FinalList(Vi),1)+1);
                tmptmp_i = tmptmps(Vi,:);
                for j=1:length(inNeigh_i)
                    descentAmount = descentAmount - incompleteNetwork(inNeigh_i(j),FinalList(Vi))*log((incompleteNetwork(inNeigh_i(j),FinalList(Vi))+1e-5)/(V(FinalList(Vi),:)*U(inNeigh_i(j),:)'+1e-5)) + incompleteNetwork(inNeigh_i(j),FinalList(Vi))*log((incompleteNetwork(inNeigh_i(j),FinalList(Vi))+1e-5)/(tmptmp_i*U(inNeigh_i(j),:)'+1e-5));
                    traceObjVal = traceObjVal + incompleteNetwork(inNeigh_i(j),FinalList(Vi))*log((incompleteNetwork(inNeigh_i(j),FinalList(Vi))+1e-5)/(V(FinalList(Vi),:)*U(inNeigh_i(j),:)'+1e-5)) - incompleteNetwork(inNeigh_i(j),FinalList(Vi))*log((incompleteNetwork(inNeigh_i(j),FinalList(Vi))+1e-5)/(tmptmp_i*U(inNeigh_i(j),:)'+1e-5));
                end
                descentAmount = descentAmount + tmptmp_i*sumLatentU' - V(FinalList(Vi),:)*sumLatentU';
                traceObjVal = traceObjVal - tmptmp_i*sumLatentU' + V(FinalList(Vi),:)*sumLatentU';
            end
            traceTime2 = traceTime2 + toc(tmpTraceTime);

            traceTime = toc(time);
            TraceU=[TraceU [traceTime-traceTime2;iterations;traceObjVal]];
        end
    end
    matlabpool close;
    time2=tic;

    for i=1:size(unknownsA,1)
        unknownsA(i,5)=U(unknownsA(i,2),:)*V(unknownsA(i,3),:)'; 
    end;
    disp(['AUC: ' num2str(scoreAUC(unknownsA(:,4),unknownsA(:,5)))]);

    toc(time);
    disp(num2str(extime));

end

function hout=suptitle(str)%used to put title on subplots
    %SUPTITLE Puts a title above all subplots.
    %	SUPTITLE('text') adds text to the top of the figure
    %	above all subplots (a "super title"). Use this function
    %	after all subplot commands.

    % Drea Thomas 6/15/95 drea@mathworks.com

    % Warning: If the figure or axis units are non-default, this
    % will break.

    % Parameters used to position the supertitle.

    % Amount of the figure window devoted to subplots
    plotregion = .92;

    % Y position of title in normalized coordinates
    titleypos  = .95;

    % Fontsize for supertitle
    fs = get(gcf,'defaultaxesfontsize')+4;

    % Fudge factor to adjust y spacing between subplots
    fudge=1;

    haold = gca;
    figunits = get(gcf,'units');

    % Get the (approximate) difference between full height (plot + title
    % + xlabel) and bounding rectangle.

        if (~strcmp(figunits,'pixels')),
            set(gcf,'units','pixels');
            pos = get(gcf,'position');
            set(gcf,'units',figunits);
        else,
            pos = get(gcf,'position');
        end
        ff = (fs-4)*1.27*5/pos(4)*fudge;

            % The 5 here reflects about 3 characters of height below
            % an axis and 2 above. 1.27 is pixels per point.

    % Determine the bounding rectange for all the plots

    % h = findobj('Type','axes');   

    % findobj is a 4.2 thing.. if you don't have 4.2 comment out
    % the next line and uncomment the following block.

    h = findobj(gcf,'Type','axes');  % Change suggested by Stacy J. Hills

    % If you don't have 4.2, use this code instead
    %ch = get(gcf,'children');
    %h=[];
    %for i=1:length(ch),
    %  if strcmp(get(ch(i),'type'),'axes'),
    %    h=[h,ch(i)];
    %  end
    %end




    max_y=0;
    min_y=1;

    oldtitle =0;
    for i=1:length(h),
        if (~strcmp(get(h(i),'Tag'),'suptitle')),
            pos=get(h(i),'pos');
            if (pos(2) < min_y), min_y=pos(2)-ff/5*3;end;
            if (pos(4)+pos(2) > max_y), max_y=pos(4)+pos(2)+ff/5*2;end;
        else
            oldtitle = h(i);
        end
    end

    if max_y > plotregion,
        scale = (plotregion-min_y)/(max_y-min_y);
        for i=1:length(h),
            pos = get(h(i),'position');
            pos(2) = (pos(2)-min_y)*scale+min_y;
            pos(4) = pos(4)*scale-(1-scale)*ff/5*3;
            set(h(i),'position',pos);
        end
    end

    np = get(gcf,'nextplot');
    set(gcf,'nextplot','add');
    if (oldtitle),
        delete(oldtitle);
    end
    ha=axes('pos',[0 1 1 1],'visible','off','Tag','suptitle');
    ht=text(.5,titleypos-1,str);set(ht,'horizontalalignment','center','fontsize',fs);
    set(gcf,'nextplot',np);
    axes(haold);
    if nargout,
        hout=ht;
    end
end