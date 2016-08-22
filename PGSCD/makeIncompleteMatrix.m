function [incompleteA, unknownsA]=makeIncompleteMatrix(incompleteA,unknownPercent,removeOnesOnly)
    count = 1;
    if(removeOnesOnly == 1)
        oneInds = find(incompleteA==1);
        unknownsA = zeros(floor((2*(unknownPercent/100) * (length(oneInds)))),4);
        for i=1:((unknownPercent/100) * (length(oneInds)))
            entry = floor(rand*(length(oneInds)))+1;
            if(incompleteA(oneInds(entry))~=-1)
                row = mod(oneInds(entry),size(incompleteA,1));
                col = floor(oneInds(entry)/size(incompleteA,1))+1;
                if(row == 0)
                    row = size(incompleteA,1); 
                    col = col - 1;
                end
                
                unknownsA(count,:) = [oneInds(entry) row col incompleteA(oneInds(entry))];
                incompleteA(oneInds(entry))=-1;
                count = count + 1;
                if(col ~= row)
                    unknownsA(count,:) = [((row-1)*size(incompleteA,1)+col) col row incompleteA(col,row)];
                    incompleteA(col,row)=-1;
                    count = count + 1;
                end
                if(mod(i,100000)==0)
                    disp(num2str(i));
                end
            end
        end
    else
        unknownsA = zeros(floor(2*((unknownPercent/100) * (size(incompleteA,1)*size(incompleteA,1)))),4);
        for i=1:((unknownPercent/100) * (size(incompleteA,1)*size(incompleteA,1)))
            entry = floor(rand*(size(incompleteA,1)*size(incompleteA,1)))+1;
            if(incompleteA(entry)~=-1)
                row = mod(entry,size(incompleteA,1));
                col = floor(entry/size(incompleteA,1))+1;
                if(row == 0) 
                    row = size(incompleteA,1); 
                    col = col - 1;
                end
                unknownsA(count,:) = [entry row col incompleteA(entry)];
                incompleteA(entry)=-1;
                count = count + 1;
                if(col ~= row)
                    unknownsA(count,:) = [((row-1)*size(incompleteA,1)+col) col row incompleteA(col,row)];
                    incompleteA(col,row)=-1;
                    count = count + 1;
                end
                if(mod(i,100000)==0)
                    disp(num2str(i));
                end
            end
        end
    end
    incompleteA(find(incompleteA==-1))=0;
    for i = size(unknownsA,1):-1:1
        if(unknownsA(i,1) ~= 0)
            unknownsA=unknownsA(1:i,:);
            break;
        end
    end
    
end