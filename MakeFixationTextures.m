% function [FixTex] = MakeFixationTextures(p,ps,RDK)
% % create precue and cue fixation "crosses"
% %   default - fixation cross normal
% %   event 1 - fixation cross event shorter left
% %   event 2 - fixation cross event longer left
% %   event 3 - fixation cross event shorter right
% %   event 4 - fixation cross event longer right
% %   cue 1 - fixation cross attend left color 1
% %   cue 2 - fixation cross attend right color 1
% %   [Will extend if more than two colors are defined]
% 
% %% create precue fixation cross textures
% 
% % default - fixation cross normal
% tex_default = zeros(p.crs.size*2 +  p.crs.width-1, p.crs.size ); % square to draw arrow texture in; added event pixels (square must always be same size); added linewidth downards in y direction; 
% x_right = [1:1:p.crs.size p.crs.size:-1:1]; % x values for right pointing arrow; 
% x_left = [p.crs.size:-1:1 1:1:p.crs.size]; % x values for left pointing arrow; 
% 
% for y = [1 : p.crs.size*2]
%     tex_default(y:y+p.crs.width-1,x_right(y)) = 1; %draw x pixel(s) for every y (+ line width) for right pointing arrow
%     tex_default(y:y+p.crs.width-1,x_left(y)) = 1; %draw x pixel(s) for every y (+ line width) for left pointing arrow
% end
% 
% tex_col = repmat(tex_default,[1 1 4]) .* permute(p.crs.color,[3 1 2]);
% FixTex.default = Screen('MakeTexture', ps.window, tex_col);
% 
% % figure; imagesc(tex_default(:,:)) % test figure
% % imagesc(tex_col(:,:,1:3));
% % Screen('DrawTextures', ps.window, FixTex.default); Screen('Flip', ps.window, 0); % test drawing
% 
% 
% 
% 
% % figure; imagesc(tex_temp(:,:)) % test figure
% % Screen('DrawTextures', ps.window, FixTex.event(i_fixev)); Screen('Flip', ps.window, 0); % test drawing
% 
% %% create cue fixation cross textures
% i_fixcue = 0;
% 
% % cue left 
% tex_uncued = tex_default;
% %tex_cued = zeros(p.crs.size*2 + p.crs.width-1, p.crs.size);
% tex_cued = zeros(p.crs.size*2 + p.crs.width-1, p.crs.size);
% for y = [1 : p.crs.size*2]
%     tex_uncued(y:y+p.crs.width-1,x_left(y)) = 0;
%     tex_cued(y:y+p.crs.width-1,x_left(y)) = 1;    
% end
% 
% for i_col = 1:1 % for every color
%     tex_uncued_col = repmat(tex_uncued,[1 1 4]) .* permute(p.crs.color,[3 1 2]); % multiply with RGBA cross color %ToDo: use isoluminant colors?
%     tex_cued_col = repmat(tex_cued,[1 1 4]) .* permute(RDK.RDK(i_col).col(1,:),[3 1 2]); % multiply with RGBA RDK color
%     tex_col = tex_uncued_col + tex_cued_col;
%     i_fixcue = i_fixcue + 1;
%     FixTex.cue(i_fixcue) = Screen('MakeTexture', ps.window, tex_col);
% end
% 
% % figure; imagesc(tex_uncued(:,:)) % test figure
% % Screen('DrawTextures', ps.window, FixTex.cue(2)); Screen('Flip', ps.window, 0); % test drawing
% 
% % cue right
% tex_uncued = tex_default;
% %tex_cued = zeros(p.crs.size*2 + p.crs.width-1, p.crs.size);
% tex_cued = zeros(p.crs.size*2 + p.crs.width-1, p.crs.size);
% for y = [1 : p.crs.size*2]
%     tex_uncued(y:y+p.crs.width-1,x_right(y)) = 0;
%     tex_cued(y:y+p.crs.width-1,x_right(y)) = 1;    
% end
% 
% for i_col = 1:1 % for every color
%     tex_uncued_col = repmat(tex_uncued,[1 1 4]) .* permute(p.crs.color,[3 1 2]); % multiply with RGBA cross color %ToDo: use isoluminant colors?
%     tex_cued_col = repmat(tex_cued,[1 1 4]) .* permute(RDK.RDK(i_col).col(1,:),[3 1 2]); % multiply with RGBA RDK color
%     tex_col = tex_uncued_col + tex_cued_col;
%     i_fixcue = i_fixcue + 1;
%     FixTex.cue(i_fixcue) = Screen('MakeTexture', ps.window, tex_col);
% end
% 
% % figure; imagesc(tex_uncued(:,:)) % test figure
% % Screen('DrawTextures', ps.window, FixTex.cue(4)); Screen('Flip', ps.window, 0); % test drawing
% 
% end






function [FixTex] = MakeFixationTextures(p,ps,RDK)
% create precue and cue fixation "crosses"
%   default - fixation cross normal
%   event 1 - fixation cross event shorter left
%   event 2 - fixation cross event longer left
%   event 3 - fixation cross event shorter right
%   event 4 - fixation cross event longer right
%   cue 1 - fixation cross attend left color 1
%   cue 2 - fixation cross attend left color 2
%   cue 3 - fixation cross attend right color 1
%   cue 4 - fixation cross attend right color 2
%   [Will extend if more than two colors are defined]


%% create precue fixation cross textures


% default - fixation cross normal
tex_default = zeros(p.crs.size*2 + p.crs.event*2 + p.crs.width-1, p.crs.size + p.crs.event*2); % square to draw arrow texture in; added event pixels (square must always be same size); added linewidth downards in y direction; 
x_right = [1:1:p.crs.size+p.crs.event p.crs.size+p.crs.event:-1:1]; % x values for right pointing arrow; 
x_left = [p.crs.size+p.crs.event:-1:1 1:1:p.crs.size+p.crs.event]+p.crs.event; % x values for left pointing arrow; 


for y = [1 : p.crs.size*2]+p.crs.event
tex_default(y:y+p.crs.width-1,x_right(y)) = 1; %draw x pixel(s) for every y (+ line width) for right pointing arrow
tex_default(y:y+p.crs.width-1,x_left(y)) = 1; %draw x pixel(s) for every y (+ line width) for left pointing arrow
end


tex_col = repmat(tex_default,[1 1 4]) .* permute(p.crs.color,[3 1 2]);
FixTex.default = Screen('MakeTexture', ps.window, tex_col);


% figure; imagesc(tex_default(:,:)) % test figure
% imagesc(tex_col(:,:,1:3));
% Screen('FillRect', ps.window, [1 1 1 1], [10 10 20 20])
% Screen('DrawTextures', ps.window, FixTex.default); Screen('Flip', ps.window, 1); % test drawing
% Screen('DrawTextures', ps.window,FixTex.default, [],[10 10 47 44]); Screen('Flip', ps.window, 0); % test drawing

% % trouble shooting
% % create a white rectangle texture
% rectangletex = zeros(6,10);
% rectangletex(2:end-1,2:end-1) = ones(4,8)*255;
% % make it a texture
% rectTexture = Screen('MakeTexture', ps.window, rectangletex); 
% Screen('DrawTexture', ps.window, rectTexture, [], [], 0);



% event 1 - fixation cross event shorter left
tex_temp = tex_default;
for y = [[1: p.crs.event]+p.crs.event [p.crs.size*2-p.crs.event+1 : p.crs.size*2]+p.crs.event]
tex_temp(y:y+p.crs.width-1,x_right(y)) = 0; %delete pixels from default tex_mat
end


tex_col = repmat(tex_temp,[1 1 4]) .* permute(p.crs.color,[3 1 2]);
i_fixev = 1;
FixTex.event(i_fixev) = Screen('MakeTexture', ps.window, tex_col);


% figure; imagesc(tex_temp(:,:)) % test figure
% imagesc(tex_col(:,:,1:3));
% Screen('DrawTextures', ps.window, FixTex.event(i_fixev)); Screen('Flip', ps.window, 0); % test drawing




% event 2 - fixation cross event longer left
tex_temp = tex_default;
for y = [1:p.crs.event [p.crs.size*2+p.crs.event+1 : p.crs.size*2+p.crs.event*2]]
tex_temp(y:y+p.crs.width-1,x_right(y)) = 1; %delete pixels from default tex_mat
end


tex_col = repmat(tex_temp,[1 1 4]) .* permute(p.crs.color,[3 1 2]);
i_fixev = i_fixev + 1;
FixTex.event(i_fixev) = Screen('MakeTexture', ps.window, tex_col);


% figure; imagesc(tex_temp(:,:)) % test figure
% Screen('DrawTextures', ps.window, FixTex.event(i_fixev)); Screen('Flip', ps.window, 0); % test drawing




% event 3 -fixation cross event shorter right
tex_temp = tex_default;
for y = [[1: p.crs.event]+p.crs.event [p.crs.size*2-p.crs.event+1 : p.crs.size*2]+p.crs.event]
tex_temp(y:y+p.crs.width-1,x_left(y)) = 0; %delete pixels from default tex_mat
end


tex_col = repmat(tex_temp,[1 1 4]) .* permute(p.crs.color,[3 1 2]); % multiply with RGBA cross color
i_fixev = i_fixev + 1;
FixTex.event(i_fixev) = Screen('MakeTexture', ps.window, tex_col);




% figure; imagesc(tex_temp(:,:)) % test figure
% Screen('DrawTextures', ps.window, FixTex.event(i_fixev)); Screen('Flip', ps.window, 0); % test drawing


% event 4 - fixation cross event longer right
tex_temp = tex_default;
for y = [1:p.crs.event [p.crs.size*2+p.crs.event+1 : p.crs.size*2+p.crs.event*2]]
tex_temp(y:y+p.crs.width-1,x_left(y)) = 1; %delete pixels from default tex_mat
end


tex_col = repmat(tex_temp,[1 1 4]) .* permute(p.crs.color,[3 1 2]);
i_fixev = i_fixev + 1;
FixTex.event(i_fixev) = Screen('MakeTexture', ps.window, tex_col);


% figure; imagesc(tex_temp(:,:)) % test figure
% Screen('DrawTextures', ps.window, FixTex.event(i_fixev)); Screen('Flip', ps.window, 0); % test drawing


%% create cue fixation cross textures
i_fixcue = 0;


% cue left 
tex_uncued = tex_default;
%tex_cued = zeros(p.crs.size*2 + p.crs.width-1, p.crs.size);
tex_cued = zeros(p.crs.size*2 + p.crs.event*2 + p.crs.width-1, p.crs.size + p.crs.event*2);
for y = [1 : p.crs.size*2]+p.crs.event
    tex_uncued(y:y+p.crs.width-1,x_left(y)) = 0;
    tex_cued(y:y+p.crs.width-1,x_left(y)) = 1;    
end


for i_col = 1:1 % for every color
    tex_uncued_col = repmat(tex_uncued,[1 1 4]) .* permute(p.crs.color,[3 1 2]); % multiply with RGBA cross color %ToDo: use isoluminant colors?
    tex_cued_col = repmat(tex_cued,[1 1 4]) .* permute(RDK.RDK(i_col).col(1,:),[3 1 2]); % multiply with RGBA RDK color
    tex_col = tex_uncued_col + tex_cued_col;
    i_fixcue = i_fixcue + 1;
    FixTex.cue(i_fixcue) = Screen('MakeTexture', ps.window, tex_col);
end


% figure; imagesc(tex_uncued(:,:)) % test figure
% Screen('DrawTextures', ps.window, FixTex.cue(2)); Screen('Flip', ps.window, 0); % test drawing


% cue right
tex_uncued = tex_default;
%tex_cued = zeros(p.crs.size*2 + p.crs.width-1, p.crs.size);
tex_cued = zeros(p.crs.size*2 + p.crs.event*2 + p.crs.width-1, p.crs.size + p.crs.event*2);
for y = [1 : p.crs.size*2]+p.crs.event
    tex_uncued(y:y+p.crs.width-1,x_right(y)) = 0;
    tex_cued(y:y+p.crs.width-1,x_right(y)) = 1;    
end


for i_col = 1:1 % for every color
    tex_uncued_col = repmat(tex_uncued,[1 1 4]) .* permute(p.crs.color,[3 1 2]); % multiply with RGBA cross color %ToDo: use isoluminant colors?
    tex_cued_col = repmat(tex_cued,[1 1 4]) .* permute(RDK.RDK(i_col).col(1,:),[3 1 2]); % multiply with RGBA RDK color
    tex_col = tex_uncued_col + tex_cued_col;
    i_fixcue = i_fixcue + 1;
    FixTex.cue(i_fixcue) = Screen('MakeTexture', ps.window, tex_col);
end


% figure; imagesc(tex_uncued(:,:)) % test figure
% Screen('DrawTextures', ps.window, FixTex.cue(4)); Screen('Flip', ps.window, 0); % test drawing


end
