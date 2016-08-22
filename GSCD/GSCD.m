function [U,V,TraceU,TraceV]=GSCD(incompleteNetwork, unknownsA, latentDim, coordGradMaxIter, gradMaxIter, gradeps,cyclic)
%Inputs:
%   incompleteNetwork - incomplete matrix
%   unknownsA - list of unknown entries
%   latentDim - latent dimension
%   coordGradMaxIter - coordinate descent iterations
%   gradMaxIter - gradient descent iterations
%   gradeps - gradient descent epsilon
%   cyclic - 1=cyclic optimization, 1=greedy optimization
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
    maxDegree = 0;
    
    gradTime = 0;
    updateTime = 0;
    
    %Initialization
    U = rand(N,latentDim);
    V = rand(M,latentDim);
    sumLatentU = 0;
    sumLatentV = 0;
    
    %Sum of all latent factors - one time calculation for all updates
    residualGradsU = zeros(N,latentDim);
    residualNormMaxHeapU = zeros(N,1);
    ResidualNormMaxHeapIndicesU = (1:N)';
    residualNormIndicesU = [(1:N)' (1:N)'];
    for i=1:N
        sumLatentU = sumLatentU + U(i,:);
        outNeigh_i = find(incompleteNetwork(i,:) > 0);%Neighbour finding
        if(isempty(outNeigh_i)==1)
            outNeigh_i = [];
        end
        if(maxDegree < length(outNeigh_i))
            maxDegree = length(outNeigh_i);
        end
        for j=1:length(outNeigh_i)%Iterate over neighbours for objective
            traceObjVal = traceObjVal + log((1+1e-5)/(U(i,:)*V(outNeigh_i(j),:)'+1e-5)) - 1;
        end
    end
            
    residualGradsV = zeros(M,latentDim);
    residualNormMaxHeapV = zeros(M,1);
    ResidualNormMaxHeapIndicesV = (1:M)';
    residualNormIndicesV = [(1:M)' (1:M)'];
    for i=1:M
        sumLatentV = sumLatentV + V(i,:);
    end
    
    for i=1:N
        traceObjVal = traceObjVal + U(i,:)*sumLatentV';
    end
    
    if(cyclic == 1)
        coordConst = coordGradMaxIter;%number of cyclic coordinate descents
    else
        coordConst = 0;%number of cyclic coordinate descents
    end
    emptyU = 0;
    emptyV = 0;
    for coordCount = 1:coordGradMaxIter%Coordinate Descent
        disp(['coordDes ' num2str(coordCount)]);
        if(isempty(ResidualNormMaxHeapIndicesU))
            emptyU = 1;
        else
            disp(['MaxNormU ' num2str(residualNormMaxHeapU(1,1))]);
        end
        if(isempty(ResidualNormMaxHeapIndicesV))
            emptyV = 1;
        else
            disp(['MaxNormV ' num2str(residualNormMaxHeapV(1,1))]);
        end
        if(emptyU == 1 && emptyV == 1)
            ResidualNormMaxHeapIndicesV = (1:M)';
            residualNormIndicesV = [(1:M)' (1:M)'];
            for i=1:M
                inNeigh_i = find(incompleteNetwork(:,i) > 0);%Neighbour finding
                if(isempty(inNeigh_i)==1)
                    inNeigh_i = [];
                end
                tmpResNorm_i = 0;
                for j=1:length(inNeigh_i)%Iterate over neighbours for objective
                    tmpResNorm_i = tmpResNorm_i - U(inNeigh_i(j),:)*(1/(V(i,:)*U(inNeigh_i(j),:)'));
                end
                residualGradsV(i,:) = sumLatentU + tmpResNorm_i;
                residualNormMaxHeapV(i,1) = norm(residualGradsV(i,:));
            end
            [residualNormMaxHeapV, ResidualNormMaxHeapIndicesV, residualNormIndicesV]=maxHeapify(residualNormMaxHeapV, ResidualNormMaxHeapIndicesV, residualNormIndicesV);
            ResidualNormMaxHeapIndicesU = (1:N)';
            residualNormIndicesU = [(1:N)' (1:N)'];
            for i=1:N
                outNeigh_i = find(incompleteNetwork(i,:) > 0);%Neighbour finding
                if(isempty(outNeigh_i)==1)
                    outNeigh_i = [];
                end
                tmpResNorm_i = 0;
                for j=1:length(outNeigh_i)%Iterate over neighbours for objective
                    tmpResNorm_i = tmpResNorm_i - V(outNeigh_i(j),:)*(1/(U(i,:)*V(outNeigh_i(j),:)'));
                end
                residualGradsU(i,:) = sumLatentV + tmpResNorm_i;
                residualNormMaxHeapU(i,1) = norm(residualGradsU(i,:));
            end
            [residualNormMaxHeapU, ResidualNormMaxHeapIndicesU, residualNormIndicesU]=maxHeapify(residualNormMaxHeapU, ResidualNormMaxHeapIndicesU, residualNormIndicesU);
            emptyU = 0;
            emptyV = 0;
        end
        if(emptyU == 0)
            for i=1:size(residualNormMaxHeapU,1)
                outNeigh_i = find(incompleteNetwork(ResidualNormMaxHeapIndicesU(i,1),:) > 0);%Neighbour finding
                if(isempty(outNeigh_i)==1)
                    outNeigh_i = [];
                end
                tmpResNorm_i = 0;
                for j=1:length(outNeigh_i)%Iterate over neighbours for objective
                    tmpResNorm_i = tmpResNorm_i - V(outNeigh_i(j),:)*(1/(U(ResidualNormMaxHeapIndicesU(i,1),:)*V(outNeigh_i(j),:)'));
                end
                residualGradsU(ResidualNormMaxHeapIndicesU(i,1),:) = sumLatentV + tmpResNorm_i;
                residualNormMaxHeapU(i,1) = norm(residualGradsU(ResidualNormMaxHeapIndicesU(i,1),:));
            end
            [residualNormMaxHeapU, ResidualNormMaxHeapIndicesU, residualNormIndicesU]=maxHeapify(residualNormMaxHeapU, ResidualNormMaxHeapIndicesU, residualNormIndicesU);
            i = 0;
            itered = 0;
            while(true)
                itered = itered + 1;
                if(itered > N)
                    break;
                end
                if(isempty(ResidualNormMaxHeapIndicesU))
                    emptyU = 1;
                end
                if(emptyU == 1)
                    break;
                end
                if(coordCount<=coordConst)%First cycle to initialize residual grads
                    i = i + 1;
                    if(i > N)
                        break;
                    end
                else
                    i = ResidualNormMaxHeapIndicesU(1,1);
                    if(residualNormMaxHeapU(1,1) < gradeps)
                        break;
                    end
                end
                outNeigh_i = find(incompleteNetwork(i,:) > 0);%Neighbour finding
                if(isempty(outNeigh_i)==1)
                    outNeigh_i = [];
                end
                tmptmp_i = U(i,:);
                tmpGradTime = tic();
                [tmpU,finalResNorm_i]=gradDes2(V(outNeigh_i,:),U(i,:),sumLatentV,gradMaxIter,gradeps);
                gradTime = gradTime + toc(tmpGradTime);
                if(norm(tmpU - U(i,:))/norm(U(i,:)) > 10^-1)
                    tisvalcount = 1;
                    U(i,:) = tmpU;
                else
                    tisvalcount = 0;
                end
                sumLatentU = sumLatentU - tmptmp_i + U(i,:);
                if(tisvalcount == 0)%No descent has occured, no need to continue coordinate descent
                    if(coordCount>coordConst)
                        if(residualNormIndicesU(i,2) ~= -1)
                            tmpUpdateTime = tic();
                            [residualNormMaxHeapU, ResidualNormMaxHeapIndicesU, residualNormIndicesU]=deleteMaxHeap(residualNormMaxHeapU, ResidualNormMaxHeapIndicesU, residualNormIndicesU, residualNormIndicesU(i,2));
                            updateTime = updateTime + toc(tmpUpdateTime);
                        end
                    end
                else
                    descentAmount = 0;
                    for j=1:length(outNeigh_i)
                        descentAmount = descentAmount - log((1+1e-5)/(U(i,:)*V(outNeigh_i(j),:)'+1e-5)) + log((1+1e-5)/(tmptmp_i*V(outNeigh_i(j),:)'+1e-5));
                        traceObjVal = traceObjVal + log((1+1e-5)/(U(i,:)*V(outNeigh_i(j),:)'+1e-5)) - log((1+1e-5)/(tmptmp_i*V(outNeigh_i(j),:)'+1e-5));
                    end
                    descentAmount = descentAmount + tmptmp_i*sumLatentV' - U(i,:)*sumLatentV';
                    traceObjVal = traceObjVal - tmptmp_i*sumLatentV' + U(i,:)*sumLatentV';
                    %Update the residual norms for i
                    if(finalResNorm_i ~= 0)
                        residualGradsU(i,:) = finalResNorm_i+sumLatentV;
                        if(coordCount<=coordConst)
                            residualNormMaxHeapU(i,1) = norm(residualGradsU(i,:));
                        else
                            %Update the residual heap for i
                            tmpUpdateTime = tic();
                            [residualNormMaxHeapU, ResidualNormMaxHeapIndicesU, residualNormIndicesU] = updateMaxHeap(residualNormMaxHeapU, ResidualNormMaxHeapIndicesU, residualNormIndicesU, residualNormIndicesU(i,2), norm(residualGradsU(i,:)));
                            updateTime = updateTime + toc(tmpUpdateTime);
                        end
                    end
                end
                if(isempty(ResidualNormMaxHeapIndicesU))
                    emptyU = 1;
                end
                if(isempty(ResidualNormMaxHeapIndicesV))
                    emptyV = 1;
                end
                traceTime = toc(time);

                if(emptyU == 0)
                    trmaxU = norm(residualGradsU(ResidualNormMaxHeapIndicesU(1,1),:));
                    TraceU=[TraceU [traceTime;traceObjVal;trmaxU]];
                end
                if(emptyV == 0)
                    trmaxV = norm(residualGradsV(ResidualNormMaxHeapIndicesV(1,1),:));
                    TraceV=[TraceV [traceTime;traceObjVal;trmaxV]];
                end
            end
        end
        if(emptyV == 0)
            for i=1:size(residualNormMaxHeapV,1)
                inNeigh_i = find(incompleteNetwork(:,ResidualNormMaxHeapIndicesV(i,1)) > 0);%Neighbour finding
                if(isempty(inNeigh_i)==1)
                    inNeigh_i = [];
                end
                tmpResNorm_i = 0;
                for j=1:length(inNeigh_i)%Iterate over neighbours for objective
                    tmpResNorm_i = tmpResNorm_i - U(inNeigh_i(j),:)*(1/(V(ResidualNormMaxHeapIndicesV(i,1),:)*U(inNeigh_i(j),:)'));
                end
                residualGradsV(ResidualNormMaxHeapIndicesV(i,1),:) = sumLatentU + tmpResNorm_i;
                residualNormMaxHeapV(i,1) = norm(residualGradsV(ResidualNormMaxHeapIndicesV(i,1),:));
            end
            [residualNormMaxHeapV, ResidualNormMaxHeapIndicesV, residualNormIndicesV]=maxHeapify(residualNormMaxHeapV, ResidualNormMaxHeapIndicesV, residualNormIndicesV);
            i = 0;
            itered = 0;
            while(true)
                itered = itered + 1;
                if(itered > M)
                    break;
                end
                if(isempty(ResidualNormMaxHeapIndicesV))
                    emptyV = 1;
                end
                if(emptyV == 1)
                    break;
                end
                if(coordCount<=coordConst)%First cycle to initialize residual grads
                    i = i + 1;
                    if(i > M)
                        break;
                    end
                else
                    i = ResidualNormMaxHeapIndicesV(1,1);
                    if(residualNormMaxHeapV(1,1) < gradeps)
                        break;
                    end
                end
                inNeigh_i = find(incompleteNetwork(:,i) > 0);%Neighbour finding
                if(isempty(inNeigh_i)==1)
                    inNeigh_i = [];
                end
                tmptmp_i = V(i,:);
                tmpGradTime = tic();
                [tmpV,finalResNorm_i]=gradDes2(U(inNeigh_i,:),V(i,:),sumLatentU,gradMaxIter,gradeps);
                gradTime = gradTime + toc(tmpGradTime);
                if(norm(tmpV - V(i,:))/norm(V(i,:)) > 10^-1)
                    tisvalcount = 1;
                    V(i,:) = tmpV;
                else
                    tisvalcount = 0;
                end
                sumLatentV = sumLatentV - tmptmp_i + V(i,:);
                if(tisvalcount == 0)%No descent has occured, no need to continue coordinate descent
                    if(coordCount>coordConst)
                        if(residualNormIndicesV(i,2) ~= -1 && residualNormIndicesV(i,1) ~= -1)
                            tmpUpdateTime = tic();
                            [residualNormMaxHeapV, ResidualNormMaxHeapIndicesV, residualNormIndicesV]=deleteMaxHeap(residualNormMaxHeapV, ResidualNormMaxHeapIndicesV, residualNormIndicesV, residualNormIndicesV(i,2));
                            updateTime = updateTime + toc(tmpUpdateTime);
                        end
                    end
                else
                    descentAmount = 0;
                    for j=1:length(inNeigh_i)
                        descentAmount = descentAmount - log((1+1e-5)/(V(i,:)*U(inNeigh_i(j),:)'+1e-5)) + log((1+1e-5)/(tmptmp_i*U(inNeigh_i(j),:)'+1e-5));
                        traceObjVal = traceObjVal + log((1+1e-5)/(V(i,:)*U(inNeigh_i(j),:)'+1e-5)) - log((1+1e-5)/(tmptmp_i*U(inNeigh_i(j),:)'+1e-5));
                    end
                    descentAmount = descentAmount + tmptmp_i*sumLatentU' - V(i,:)*sumLatentU';
                    traceObjVal = traceObjVal - tmptmp_i*sumLatentU' + V(i,:)*sumLatentU';
                    %Update the residual norms for i
                    if(finalResNorm_i ~= 0)
                        residualGradsV(i,:) = finalResNorm_i + sumLatentU;
                        if(coordCount<=coordConst)
                            residualNormMaxHeapV(i,1) = norm(residualGradsV(i,:));
                        else
                            %Update the residual heap for i
                            tmpUpdateTime = tic();
                            [residualNormMaxHeapV, ResidualNormMaxHeapIndicesV, residualNormIndicesV] = updateMaxHeap(residualNormMaxHeapV, ResidualNormMaxHeapIndicesV, residualNormIndicesV, residualNormIndicesV(i,2), norm(residualGradsV(i,:)));
                            updateTime = updateTime + toc(tmpUpdateTime);
                        end
                    end
                end
                if(isempty(ResidualNormMaxHeapIndicesU))
                    emptyU = 1;
                end
                if(isempty(ResidualNormMaxHeapIndicesV))
                    emptyV = 1;
                end
                traceTime = toc(time);

                if(emptyU == 0)
                    trmaxU = norm(residualGradsU(ResidualNormMaxHeapIndicesU(1,1),:));
                    TraceU=[TraceU [traceTime;traceObjVal;trmaxU]];
                end
                if(emptyV == 0)
                    trmaxV = norm(residualGradsV(ResidualNormMaxHeapIndicesV(1,1),:));
                    TraceV=[TraceV [traceTime;traceObjVal;trmaxV]];
                end
            end
        end
    end
%     time2=tic;
%     
%     for i=1:size(unknownsA,1)
%         unknownsA(i,5)=U(unknownsA(i,2),:)*V(unknownsA(i,3),:)'; 
%     end;
%     disp(['AUC: ' num2str(scoreAUC(unknownsA(:,4),unknownsA(:,5)))]);

    toc(time);
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