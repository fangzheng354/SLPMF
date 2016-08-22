function [heap, heapIndArray, indArray]=maxHeapify(heap, heapIndArray, indArray)
    count = length(heap);
    start = floor((count - 1) / 2);
    while (start >= 0)
        [heap, heapIndArray, indArray] = siftDown(heap, heapIndArray, indArray, start, count-1);
        start = start - 1;
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