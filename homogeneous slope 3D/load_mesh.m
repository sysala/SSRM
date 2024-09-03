function [coord,elem,surf,Q]=load_mesh()


% MATLAB script to load datasets from HDF5 file

% File path
file_path = 'slope_paper4.h5';

% Load datasets
boundary = h5read(file_path, '/boundary');
elem = h5read(file_path, '/elem') + 1;
face = h5read(file_path, '/face') + 1;
material = h5read(file_path, '/material');
node = h5read(file_path, '/node');

Q=true(size(node));

tmp = face(:,boundary==1);
Q(1,tmp(:)) = 0;
tmp = face(:,boundary==2);
Q(1,tmp(:)) = 0;
tmp = face(:,boundary==3);
Q(2,tmp(:)) = 0;
tmp = face(:,boundary==4);
Q(2,tmp(:)) = 0;
tmp = face(:,boundary==5);
Q(3,tmp(:)) = 0;

coord = double(node([1 3 2],:));
coord(1,:) = coord(1,:)+20;
Q = Q([1 3 2],:);
elem=double(elem([1 2 4 3],:));
surf = double(face);


% % Display the size of each dataset
% disp('Size of boundary:'), disp(size(boundary));
% disp('Size of elem:'), disp(size(elem));
% disp('Size of face:'), disp(size(face));
% disp('Size of material:'), disp(size(material));
% disp('Size of node:'), disp(size(node));

% 
% % Create a new figure
% figure;
% 
% % Plot each face
% hold on;
% for i = 1:size(face, 2)
%     if boundary(i) == 0
%         continue
%     end
%     % Extract the node indices for this face
%     node_indices = face(:, i);
% 
%     % Get the coordinates of the nodes
%     vertices = node(:,node_indices);
% 
%     % Plot the face as a polygon
%     fill3(vertices(1,:), vertices(2,:), vertices(3,:), 'cyan', 'FaceAlpha', 0.9);
% 
%     % Calculate the centroid of the face for placing the boundary number
%     %centroid = mean(vertices, 2);
% 
%     % Overlay the boundary number on the plot at the centroid
%     %text(centroid(1), centroid(2), centroid(3), num2str(boundary(i)), 'FontSize', 10, 'Color', 'red', 'HorizontalAlignment', 'center');
% end
% hold off;
% 
% % Set plot properties
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('3D Plot of All Faces');
% axis equal;
% grid on;
% view(3);
end