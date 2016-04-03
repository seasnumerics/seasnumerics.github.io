% this script is to find the points of each links
% l1(k) means the first  point at the k-th link
% l2(k) means the second point at the k-th link

function [l1,l2]=links(v,nv)

% assign the links for the 1st triangle, v(1,:)
%           l1(1) = v(1,1) = l1(2)
%                    /\
%                   /  \
%                  /    \
%          l(1)=> /      \<= l(2)
%                /        \
%               /          \
%   l1(1)=v(1,2)------------ v(1,3)=l2(2)
%        =l1(3)     l(3)           =l2(3)
l1=zeros(1,nv+length(v)-2);  % refine 4 times, 7680 edges, 5120 faces, 2562 vertices;
l2=zeros(1,nv+length(v)-2);  % Recall the formula : F + V - 2 = E
l1(1) = v(1,1);
l2(1) = v(1,2);
l1(2) = v(1,1);
l2(2) = v(1,3);
l1(3) = v(1,2);
l2(3) = v(1,3);
j=3;


for i=2:length(v)
    
    k=1;
    while ( k < j+1 )
        if (( v(i,1)==l1(k) && v(i,2)==l2(k) ) ||  ( v(i,1)==l2(k) && v(i,2)==l1(k) ))
            k=1000000;
        else
            k=k+1;
        end
        
        if (k == j+1)
            j=j+1;
            l1(j) = v(i,1);
            l2(j) = v(i,2);
            k=k+1;
            %l1(j)
            %l2(j)
        end
        
    end
    
    
    k=1;
    while ( k < j+1 )
        if (( v(i,1)==l1(k) && v(i,3)==l2(k) ) ||  ( v(i,1)==l2(k) && v(i,3)==l1(k) ))
            k=1000000;
        else
            k=k+1;
        end
        
        if (k == j+1)
            j=j+1;
            l1(j) = v(i,1);
            l2(j) = v(i,3);
            k=k+1;
        end
        
    end
    
    k=1;
    while ( k < j+1 )
        if (( v(i,2)==l1(k) && v(i,3)==l2(k) ) ||  ( v(i,2)==l2(k) && v(i,3)==l1(k) ))
            k=1000000;
        else
            k=k+1;
        end
        
        if (k == j+1)
            j=j+1;
            l1(j) = v(i,2);
            l2(j) = v(i,3);
            k=k+1;
        end
    end
    
end