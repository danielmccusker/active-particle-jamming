noCells=512;
noFrames=141;

%% make one video of cell positions
A = importdata('positions.dat');
% C = mat2cell(A,noCells*ones(1,noFrames));
% 
% vidwriter = VideoWriter('positions1.avi'); 
% vidwriter.FrameRate = 10;
% vidwriter.Quality = 50;
% open(vidwriter);
% set(gca, 'Position',[0 0 1 1]);
% pbaspect([1 1 1]);
% for i=1:noFrames
%     %viscircles([C{i}(:,1) C{i}(:,2)],C{i}(:,3),'Color','black','LineStyle','-');
%     
%     circle(C{i}(:,1),C{i}(:,2),C{i}(:,3));
%         
%     writeVideo(vidwriter,getframe);
%     cla;
% end
% close(vidwriter);

%% get Fourier transform of this run
FC = C;
for i=1:noFrames
    temp = cell2mat(C(i));
    temp = abs(fftshift(fft2(temp)));
    FC(i) = mat2cell(temp,noCells);
    %FC(i) = abs(fftshift(FC(i)));
end

%% make a video of the Fourier transform
vidwriter = VideoWriter('FT1.avi'); 
vidwriter.FrameRate = 10;
open(vidwriter);
for i=1:noFrames
    surf(FC{i});
    currFrame = getframe;
    writeVideo(vidwriter,currFrame);
end
close(vidwriter);

% %% make a second position video
% A = importdata('positions2.dat');
% A = [A, ones(262144,1)];%define delta functions at each cell location
% C = mat2cell(A,noCells*ones(1,noFrames));
% vidwriter = VideoWriter('positions2.avi'); 
% vidwriter.FrameRate = 10;
% open(vidwriter);
% for i=1:256
%     scatter(C{i,1}(:,1),C{i,1}(:,2));
%     currFrame = getframe;
%     writeVideo(vidwriter,currFrame);
% end
% close(vidwriter);
% 
% %% get Fourier transform of this run
% FC = C;
% for i=1:256
%     temp = cell2mat(C(i));
%     temp = abs(fftshift(fft2(temp)));
%     FC(i) = mat2cell(temp,noCells);
%     %FC(i) = abs(fftshift(FC(i)));
% end
% 
% %% make a second video of the Fourier transform
% vidwriter = VideoWriter('FT2.avi'); 
% vidwriter.FrameRate = 10;
% open(vidwriter);
% for i=1:256
%     scatter(FC{i,1}(:,1),FC{i,1}(:,2));
%     currFrame = getframe;
%     writeVideo(vidwriter,currFrame);
% end
% close(vidwriter);

