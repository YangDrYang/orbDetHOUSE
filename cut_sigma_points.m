addpath('../CUTpoints-master/')

for n = 1:20
    
    mu = zeros(n, 1);
    P = eye(n);
    
    % 4th order points
    [X4,w4]=conjugate_dir_gausspts_4thmoments(mu,P);
    writematrix(X4', sprintf('CUT/pts_4_%d.csv', n));
    writematrix(w4,  sprintf('CUT/wts_4_%d.csv', n));
    
    % 6th order points (dimension 2 to 9)
    if (n>=2 && n<=9)
        [X6,w6]=conjugate_dir_gausspts_6moment(mu,P);
        writematrix(X6', sprintf('CUT/pts_6_%d.csv', n));
        writematrix(w6,  sprintf('CUT/wts_6_%d.csv', n));
    end
    
    % 8th order points (dimension 2 to 6)
    if (n>=2 && n<=6)
        [X8,w8]=conjugate_dir_gausspts_8moment(mu,P);
        writematrix(X8', sprintf('CUT/pts_8_%d.csv', n));
        writematrix(w8,  sprintf('CUT/wts_8_%d.csv', n));
    end
    
end

rmpath('../CUTpoints-master/')
