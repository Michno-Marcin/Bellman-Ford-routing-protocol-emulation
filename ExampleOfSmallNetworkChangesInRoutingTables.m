% Emulation of a selected routing protocol: distance-vector 
% Application of Belman-Ford algorithm to find shortest routes

%% Preparation to run the written script

clear; close all; clc;
% Optional closing the session of the previous script
child_handles = allchild(0);
names = get(child_handles,'Name');
asd = find(strncmp('Biograph Viewer', names, 15));
close(child_handles(asd));

disp('DISTANCE VECTOR ROUTING PROTOCOL:');

%% Generating a random network graph (directed)

fprintf('\t');
AmountOfNodes = input('Enter the number of nodes: '); % Number of entered nodes, e.g. 15
fprintf('\t');
AmountOfLinks = input('Enter the number of random links: '); % Number of accepted channels, e.g. 60
% Assigning exmplary weights from the range to all network channels: 1-100
weights =[5 3 4 2 7];
% Assigning channels in a network (identifying sources and outlets)
source_nodes=[1 1 2 3 4];
target_nodes=[2 4 5 5 5];

% Assigning node names
names=["A" "B" "C" "D" "E"];
% names % <- checking channel names
 
%% Visualization of the created network graph

NetworkGraph = graph(source_nodes,target_nodes,weights,names);
% The thickness of the individual channels in the drawing corresponds to their weight
LWidths = 2*NetworkGraph.Edges.Weight/max(NetworkGraph.Edges.Weight);
x=randi(AmountOfNodes*5,1,AmountOfNodes);
y=randi(AmountOfNodes*5,1,AmountOfNodes);
GPlotGraph = plot(NetworkGraph,'XData',x,'YData',y,'EdgeLabel',...
    NetworkGraph.Edges.Weight,'LineWidth',LWidths); % Graph visualization
GPlotGraph.MarkerSize = 7; % Resizing nodes for clarity
GPlotGraph.NodeColor='black'; % Changing nodes color for drawing clarity
GPlotGraph.EdgeColor='blue'; % Changing channels color for drawing clarity
% It's also possible to show the ratio of distances between given nodes in the graph
layout(GPlotGraph,'force','WeightEffect','direct') ;

%% Sample matrices of a given network graph

% AdjandencyMatrix = full(adjacency(NetworkGraph)) % Adjacency matrix of given graph
AdWeMatrix = adjacency(NetworkGraph,'weighted');
% (Adjacency) weight matrix of given graph
AdjadencyWeightMatrix = full(AdWeMatrix); fprintf("\tAdjadencyWeightMatrix:\n"); 
disp(AdjadencyWeightMatrix); 
% IMatrix = incidence(NetworkGraph);
% IncidenceMatrix = full(IMatrix) % Incidence matrix of the given graph

%% Other network graph visualization (more transparent)

BioNetGraph = biograph(triu(AdjadencyWeightMatrix),names,'showarrows','off',...
    'ShowWeights','on','EdgeTextColor',[0 0 1]);
set(BioNetGraph.Nodes,'Color',[.5 .7 1]);
set(BioNetGraph.Edges,'LineColor',[0 0 0]);
view(BioNetGraph);

%% Example of routing tables 
% 
%   | from A | via A | via B | via C | via D | via E | ...
%    ------------------------------------------------
%   |  to A  |  inf  |  inf  |  inf  |  inf  |  inf  |
%   |  to B  |  inf  |       |       |       |       |
%   |  to C  |  inf  |       |       |       |       |  
%   |  to D  |  inf  |       |       |       |       |  
%   |  to E  |  inf  |       |       |       |       |  
%
% 
%   | from B | via A | via B | via C | via D | via E | ...
%    ------------------------------------------------
%   |  to A  |       |  inf  |       |       |       | 
%   |  to B  |  inf  |  inf  |  inf  |  inf  |  inf  | 
%   |  to C  |       |  inf  |       |       |       |  
%   |  to D  |       |  inf  |       |       |       | 
%   |  to E  |       |  inf  |       |       |       |
%
% etc. ...


%% Assign an index to each non-zero node in the matrix (IDNodeMatrix)
 NoOfNode=1;
 TriuMat=triu(AdjadencyWeightMatrix); IDNodeMatrix=int16.empty;
 for i=1:AmountOfNodes
     for j=1:AmountOfNodes
         if(TriuMat(i,j)~=0)
             IDNodeMatrix(i,j)=NoOfNode;
             IDNodeMatrix(j,i)=IDNodeMatrix(i,j);
             NoOfNode=NoOfNode+1;      
         end
     end
 end
% view(BioNetGraph)

%% First stage of the algorithm to determine routing tables for network nodes
% Taking into account the paths connecting neighbors for each node

RoutingTable=int16.empty;
for from=1:AmountOfNodes
    for via=1:AmountOfNodes
        for to=1:AmountOfNodes       
            % There is no point in looking for paths that return to the start node
            if(from~=via&&from~=to)
                % In the absence of "intermediates", the weight of the path 
                % is determined by the weighted neighbourhood matrix
                if(via==to&&AdjadencyWeightMatrix(from,to)~=0) 
                    RoutingTable(to,via,from)=AdjadencyWeightMatrix(from,to);
                else
                     % We will take a relatively large enough number as the path weight
                     % for the script to work correctly, meaning infinity
                     RoutingTable(to,via,from)=intmax/10; 
               end
            else
                    % As above, there is no weighted path with the same: source, outlet
                    RoutingTable(to,via,from)=inf;
            end           
        end
    end
end

%% Further steps to determine the final routing tables for nodes in the created network
% Application of the Belman-Ford algorithm
% The following loops execute as long as the routing tables are changed in them.
% That is, the shortest paths between nodes are updated as long as there are routing tables.
temp=int32.empty;
while(i>=0)    
    OldRoutingTable=RoutingTable;
    for from=1:AmountOfNodes
        for to=1:AmountOfNodes
            if(from~=to)   % There is no point in looking for paths that return to the start node
                if(AdjadencyWeightMatrix(from,to)~=0)  % Determination of adjacent nodes
                    for x=1:AmountOfNodes
                        for y=1:AmountOfNodes                    
                            % Updating routing tables by selecting better paths
                            temp(x,y)=AdjadencyWeightMatrix(from,to)+min(RoutingTable(y,:,to));
                            if(temp(x,y)<RoutingTable(y,to,from)&&RoutingTable(y,to,from)<inf)
                                RoutingTable(y,to,from)=temp(x,y);
                            end
                        end
                    end
                end
            end
        end
    end
    fprintf('Routing tables after %d iterations of Bellman - Ford algorithm: \n', i);
    disp(RoutingTable);
    i=i+1;
    % Once the routing tables are fixed, we can determine the optimal routing paths
    if(OldRoutingTable==RoutingTable) 
        fprintf('In order to determine the final routing tables %d iterations of the algorithm (Bellman - Ford) were performed. \n', i)
        break;
    end
end
% disp(path); <- Displaying routing tables; % e.g. "disp(path(:,:,2))" <- for node number 2

%% Determination of optimal paths between selected source and outlet
choice='y';
while(choice=='y')
fprintf('\t');
source=input('Enter the source node: ');
fprintf('\t');
destination=input('Enter the destination node: ');

% Checking whether there is a route between selected nodes or not
AreConnectedNodes = conncomp(NetworkGraph);
trace = int16.empty;
if ((destination > AmountOfNodes) || (destination < 1) || mod(destination,1))
    fprintf('\tNo destination nodes with these identifiers ! \n\tProper numbers are:');
    disp(1:AmountOfNodes) ;
    fprintf('\t');
    choice=input('Do you want to look for shortest paths again (y/n): ','s');
    continue;
    % Prompt about a possible wrong choice of target nodes  
else
    if ((source > AmountOfNodes) || (destination < 1) || mod(destination,1))
        fprintf('No destination nodes with these identifiers ! \nProper numbers are:');
        disp(1:AmountOfNodes) ;
        fprintf('\t');
        choice=input('Do you want to look for shortest paths again (y/n): ','s');
        continue;
        % Prompt about a possible wrong choice of starting nodes
    else
        if (AreConnectedNodes(source) == AreConnectedNodes(destination))
            trace(1)=source;
            j=2;
            % Selecting the best route between selected nodes
            while(source~=destination)
                [row,col] = find(RoutingTable(destination,:,source) ==...
                    min(RoutingTable(destination,:,source)));
                trace(j)=col;
                source=col;
                j=j+1;
            end
            % Displaying optimal path in comments and drawing
            disp('Best route between these nodes is:');
            disp(trace(1:j-1));
            BioNetGraph=biograph(triu(AdjadencyWeightMatrix),names,'showarrows','off',...
                'ShowWeights','on','EdgeTextColor',[0 0 1]);
            for i=1:j-1
                if (i<j-1)
                    if (i==1)
                        set(BioNetGraph.nodes(trace(i)),'color',[.5 .7 1]);
                        set(BioNetGraph.edges(IDNodeMatrix(trace(i),trace(i+1))),'linecolor',[1 0 0]);
                    else
                        set(BioNetGraph.nodes(trace(i)), 'color', [0 1 0]);
                        set(BioNetGraph.edges(IDNodeMatrix(trace(i),trace(i+1))),...
                            'linecolor',[1 0 0]);
                    end
                else
                    set(BioNetGraph.nodes(trace(i)),'color',[.5 .7 1]);
                end
            end
            view(BioNetGraph);
            choice=input('Do you want to look for shortest paths again (y/n): ','s');
        else 
            disp('Paths between these nodes does not exist.'); 
            % Prompt about missing connection between given nodes, if any
            choice=input('Do you want to look for shortest paths again (y/n): ','s');
        end
    end
end
end