function event = get_z_hmm(data)
    
    Niter = 20;
    Nburn = 0;
    Nplot = 0;

    data = abs([0 diff(data)']);
    data = reshape(data, 96, []);
    event_times = zeros(size(data));
    res = mmpp(data, [Niter,Nburn,Nplot], event_times, [3,3]);

    event = res.Z(:,:,end);
    event = event(:)';