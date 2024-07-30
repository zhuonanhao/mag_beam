clear all;
close all;
clc;

l1 = 1.0;
%l2 = 1.0;

n1 = 5;
%n2 = 11;

% delta = 0.001;

%nv1 = l1 / delta + 1;
%nv2 = l2 / delta + 1;

dl1 = l1 / (n1 - 1);
%dl2 = l2 / (n2 + 1);

bigNodes = zeros(3,2);

temp = 1;
for i = 1:n1
    for j = 1:n1
        
        %logo = 1;
        
%         if (i == 1 && j == 1)
%             logo = -1;
%         end
%         
%         if (i == 1 && j == n1)
%             logo = -1;
%         end
%         
%         if (i == n1 && j == 1)
%             logo = -1;
%         end
%         
%         if (i == n1 && j == n1)
%             logo = -1;
%         end
        
        %if (logo > 0)
            bigNodes(temp,1) = dl1 * (j-1) - l1 / 2;
            bigNodes(temp,2) = dl1 * (i-1) - l1 / 2;
            temp = temp + 1;
        %end
    end
end

step = 1

%figure(1)
%hold on;

%plot(bigNodes(:,1),bigNodes(:,2),'o');

maxEdge = abs(dl1 * n1 - l1 / 2);

[totalBigNode, ~] = size(bigNodes);

stretchIndex = zeros(3,2);
temp = 1;
for i = 1:totalBigNode
    for j = i+1:totalBigNode
        
        node1 = bigNodes(i,:);
        node2 = bigNodes(j,:);
        
        logo = 1;
        
%         if ( abs(node1(1)-maxEdge) < 1e-3 && abs(node2(1)-maxEdge) < 1e-3 )
%             logo = -1;
%         end
%         
%         if ( abs(node1(1)+maxEdge) < 1e-3 && abs(node2(1)+maxEdge) < 1e-3 )
%             logo = -1;
%         end
%         
%         if ( abs(node1(2)-maxEdge) < 1e-3 && abs(node2(2)-maxEdge) < 1e-3 )
%             logo = -1;
%         end
%         
%         if ( abs(node1(2)+maxEdge) < 1e-3 && abs(node2(2)+maxEdge) < 1e-3 )
%             logo = -1;
%         end
        
        if ( abs(norm(node1 - node2) - dl1) < 1e-3 && logo > 0)
            stretchIndex(temp,1) = i;
            stretchIndex(temp,2) = j;
            temp = temp + 1;
        end
       
    end
end

step = 2

[stretchNum, ~] = size(stretchIndex);

% for i = 1:stretchNum
%     
%     node1 = bigNodes(stretchIndex(i,1),:);
%     node2 = bigNodes(stretchIndex(i,2),:);
%     
%     plot([node1(1) node2(1)],[node1(2) node2(2)],'r-');
% 
% end

[totalBigNode, ~] = size(bigNodes);
temp = totalBigNode + 1;

insideNode = 10;
delta = dl1 / (insideNode + 1);

nodes = bigNodes;
for i = 1:stretchNum
    node1 = bigNodes(stretchIndex(i,1),:);
    node2 = bigNodes(stretchIndex(i,2),:);
    
    t = (node2 - node1) / norm(node2 - node1);
    
    for j = 1:insideNode
        localNode = node1 + t * j * delta;
        
        nodes(temp,1) = localNode(1);
        nodes(temp,2) = localNode(2);
        temp = temp + 1;
    end
end


step = 3

%plot(nodes(:,1),nodes(:,2),'s');
[nv,~] =size(nodes);

stretchIndex = zeros(3,2);
temp = 1;
for i = 1:nv
    for j = i+1:nv
        
        node1 = nodes(i,:);
        node2 = nodes(j,:);
       
        if ( abs(norm(node1 - node2) - delta) < 1e-3)
            stretchIndex(temp,1) = i;
            stretchIndex(temp,2) = j;
           
            temp = temp + 1;
        end
       
    end
end

step = 4

%figure(2)
%hold on;

[stretchNum, ~] = size(stretchIndex);

% for i = 1:stretchNum
%     
%     node1 = nodes(stretchIndex(i,1),:);
%     node2 = nodes(stretchIndex(i,2),:);
%     
%     plot([node1(1) node2(1)],[node1(2) node2(2)],'r-');
% 
% end


bendingIndex = zeros(3,3);
temp = 1;
for i = 1:stretchNum
    for j = i+1:stretchNum
        stretchIndex1 = stretchIndex(i,:);
        stretchIndex2 = stretchIndex(j,:);
        
        if (stretchIndex1(1) == stretchIndex2(1))
            bendingIndex(temp,1) = stretchIndex1(2);
            bendingIndex(temp,2) = stretchIndex1(1);
            bendingIndex(temp,3) = stretchIndex2(2);
            
            temp = temp + 1;
        end
        
        if (stretchIndex1(1) == stretchIndex2(2))
            bendingIndex(temp,1) = stretchIndex1(2);
            bendingIndex(temp,2) = stretchIndex1(1);
            bendingIndex(temp,3) = stretchIndex2(1);
            
            temp = temp + 1;
        end
        
        if (stretchIndex1(2) == stretchIndex2(1))
            bendingIndex(temp,1) = stretchIndex1(1);
            bendingIndex(temp,2) = stretchIndex1(2);
            bendingIndex(temp,3) = stretchIndex2(2);
            
            temp = temp + 1;
        end
        
        if (stretchIndex1(2) == stretchIndex2(2))
            bendingIndex(temp,1) = stretchIndex1(1);
            bendingIndex(temp,2) = stretchIndex1(2);
            bendingIndex(temp,3) = stretchIndex2(1);
           
            temp = temp + 1;
        end
        
    end
end

plot(nodes(:,1),nodes(:,2),'s');
hold on

[nv,~] = size(nodes);

fixIndex = zeros(3,1);
temp = 1;
for i = 1:nv
    
    xx = nodes(i,1);
    yy = nodes(i,2);
    
    if (yy < -0.4)
        fixIndex(temp) = i;
        temp = temp + 1;
    end
end

plot(nodes(fixIndex,1),nodes(fixIndex,2),'ro');


dlmwrite('nodesInput.txt',nodes,'delimiter',' ');
dlmwrite('edgeInput.txt',stretchIndex-1,'delimiter',' ');
dlmwrite('bendingInput.txt',bendingIndex-1,'delimiter',' ');
dlmwrite('constraint.txt',fixIndex-1,'delimiter',' ');
