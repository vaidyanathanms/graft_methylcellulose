function [patchHndl,lineHndl] = myCustMarker(xData,yData,markerDataX,markerDataY,markerDataX2,markerDataY2,markerSize,mface1,medge1,mface2,medge2,inp_handle)

%% make sure the inputs are OK
if ~isvector(xData) || ~isvector(yData) || length(xData)~=length(yData)
    fprintf('Error! xData and yData mar must be vectors of the same length!\n');
    return
end
if ~isvector(markerDataX) || ~isvector(markerDataY) || length(markerDataX)~=length(markerDataY)
    fprintf('Error! markerDataX and markerDataY mar must be vectors of the same length!\n');
    return
end
% ------
xData = reshape(xData,length(xData),1) ;
yData = reshape(yData,length(yData),1) ;
markerDataX = markerSize * reshape(markerDataX,1,length(markerDataX)) ;
markerDataY = markerSize * reshape(markerDataY,1,length(markerDataY)) ;
markerDataX2 = markerSize * reshape(markerDataX2,1,length(markerDataX2)) ;
markerDataY2 = markerSize * reshape(markerDataY2,1,length(markerDataY2)) ;
% -------------------------------------------------------------


%% prepare and plot the patches
% markerEdgeColor = [0 0 0] ;
% markerFaceColor = [0 1 1] ;

lineStyle = 'None' ;
lineColor = [1 0 0] ;
% ------
vertX = repmat(markerDataX,length(xData),1) ;
vertX = vertX(:) ;
vertY = repmat(markerDataY,length(yData),1) ;
vertY = vertY(:) ;

vertX2 = repmat(markerDataX2,length(xData),1) ;
vertX2 = vertX2(:) ;
vertY2 = repmat(markerDataY2,length(yData),1) ;
vertY2 = vertY2(:) ;

% ------
vertX = repmat(xData,length(markerDataX),1) + vertX ;
vertY = repmat(yData,length(markerDataY),1) + vertY ;

vertX2 = repmat(xData,length(markerDataX2),1) + vertX2 ;
vertY2 = repmat(yData,length(markerDataY2),1) + vertY2 ;

% ------
faces = 0:length(xData):length(xData)*(length(markerDataY)-1) ;
faces = repmat(faces,length(xData),1) ;
faces = repmat((1:length(xData))',1,length(markerDataY)) + faces ;

faces2 = 0:length(xData):length(xData)*(length(markerDataY2)-1) ;
faces2 = repmat(faces2,length(xData),1) ;
faces2 = repmat((1:length(xData))',1,length(markerDataY2)) + faces2 ;


% ------
hold on;
figHndl = inp_handle ; box on ; 
lineHndl = plot(xData,yData,'lineStyle',lineStyle,'Color',lineColor) ;
patchHndl = patch('Faces',faces,'Vertices',[vertX vertY]);
patchHndl2 = patch('Faces',faces2,'Vertices',[vertX2 vertY2]);
set(patchHndl,'FaceColor',mface1,'EdgeColor',medge1) ;
set(patchHndl2,'FaceColor',mface2,'EdgeColor',medge2) ;
% -------------------------------------------------------------

