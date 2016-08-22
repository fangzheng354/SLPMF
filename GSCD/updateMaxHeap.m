function [heap, heapIndArray, indArray]=updateMaxHeap(heap, heapIndArray, indArray, Index, newVal)
    if (heap(Index) > newVal)
        heap(Index)=newVal;
        [heap, heapIndArray, indArray] = siftDown(heap, heapIndArray, indArray, Index-1, length(heap)-1);
    elseif(heap(Index) < newVal)
        heap(Index)=newVal;
        [heap, heapIndArray, indArray] = percolate(heap, heapIndArray, indArray, Index-1);
    end
end

function [a, b, c] = siftDown(a, b, c, start, endd)
     root = start;
     while (root * 2 + 1 <= endd)
         child = root * 2 + 1;
         swap = root;
         if (a(swap+1) < a(child+1))
             swap = child;
         end
         if (child+1 <= endd && a(swap+1) < a(child+2))
             swap = child + 1;
         end
         if (swap ~= root)
             tmp1=a(root+1);
             tmp2=b(root+1);
             a(root+1) = a(swap+1);
             b(root+1) = b(swap+1);
             a(swap+1) = tmp1;
             b(swap+1) = tmp2;
             c(b(root+1),2)=root+1;
             c(b(swap+1),2)=swap+1;
             root = swap;
         else
             return
         end
     end
end

function [a, b, c] = percolate(a, b, c, start)
    root = start;
    while(root > 0)
        if(a(floor((root-1)/2) + 1) < a(root+1))
            swap = floor((root-1)/2) + 1;
            tmp1 = a(swap);
            tmp2 = b(swap);
            a(swap) = a(root+1);
            b(swap) = b(root+1);
            a(root+1) = tmp1;
            b(root+1) = tmp2;
            c(b(root+1),2)=root+1;
            c(b(swap),2)=swap;
        else
            return;
        end
        root = floor((root-1)/2);
    end
end