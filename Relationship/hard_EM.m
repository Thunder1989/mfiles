function residual = hard_EM(Niter, data)
    N=10;    
    kf = StandardKalmanFilter(data_vav',8,N,'EWMA'); 
    diff1 = abs(data' - kf);
    diff1(isnan(diff1)) = 0;
    diff1 = diff1(16:end);

    mask = [zeros(1,15) diff1] >= mean(diff1) + 1*std(diff1);
    data_steady = data(1-mask);
    data_transition = data(mask);
    
    data_steady_kf = StandardKalmanFilter(data_steady,8,N,'EWMA'); 
    data_transition_kf = StandardKalmanFilter(data_transition,8,N,'EWMA'); 