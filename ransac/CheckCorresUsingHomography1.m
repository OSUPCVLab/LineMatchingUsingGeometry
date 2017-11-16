

function [H, inliers,OrtErr,PtDist] = CheckCorresUsingHomography1(x1, x2, th1, th2,type,normlz)

    if ~all(size(x1)==size(x2))
        error('Data sets x1 and x2 must have the same dimension');
    end
    
    [rows,npts] = size(x1);
    
    % Based on (line end-points) or (just points)
    if rows~=4 && rows~=6 && rows~=2 && rows~=3 
        error('x1 and x2 must have 3 or 6 rows');
    end
    
    
    
    if npts < 4
        error('Must have at least 4 lines to fit homography');
    end
    
    if rows == 4    % Pad data with homogeneous scale factor of 1
        
        x1 = [x1(1:2,:); ones(1,npts);x1(3:4,:); ones(1,npts)];
        x2 = [x2(1:2,:); ones(1,npts);x2(3:4,:); ones(1,npts)];
    end
    
    if rows == 2    % Pad data with homogeneous scale factor of 1
        
        x1 = [x1(1:2,:); ones(1,npts)];
        x2 = [x2(1:2,:); ones(1,npts)];
    end
%% Based on Point Normalization
    % Normalise each set of points so that the origin is at centroid and
    % mean distance from origin is sqrt(2).  normalise2dpts also ensures the
    % scale parameter is 1.  
    
    if (strcmp('Points', normlz))
        [x1, T1,T1s] = normalise2dpts(x1);
        [x2, T2,T2s] = normalise2dpts(x2);
        T12=[T1;T2];
    
        [H] = NewHomography2d_1(x1, x2,T12,type);

        % Distance evaluation 
        [H, inliers,OrtErr,PtDist] = homogdist2d1(H, [x1;x2], th1, th2, T12,type);
    end
%% Based on Line normalization
    if (strcmp('Lines', normlz))
        
        Ls1 = cross(x1(1:3,:), x1(4:6,:));
        Ls2 = cross(x2(1:3,:), x2(4:6,:));
        
        [L1, T1,T1s] = normalise2dlines(Ls1);
        [L2, T2,T2s] = normalise2dlines(Ls2);
        T12=[T1;T2];
    
        [H] = NewHomography2d_2L(L1, L2,type);
    
        % Distance evaluation 
        [H, inliers,OrtErr,PtDist] = homogdist2d2L(H, [x1;x2], [L1;L2], th1, th2, T12,type);
    end
    

end
    