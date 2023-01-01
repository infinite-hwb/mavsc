function [CA F P R nmi AR,missrate] = get_result(A,s)
for q = 1:20
    gnd=s;
    grp = SpectralClustering(A, max(s));
    grps = bestMap(s,grp);
    missrate = sum(s(:) ~= grps(:)) / length(s);
    CAq(q) = 1-compute_CE(grps, gnd); % clustering accuracy
    [Fq(q),Pq(q),Rq(q)] = compute_f(gnd,grps); % F1, precision, recall
    nmiq(q) = compute_nmi(gnd,grps);
    ARq(q) = rand_index(gnd,grps);
end
    CA(1)=mean(CAq); CA(2)=std(CAq);
    F(1)=mean(Fq); F(2)=std(Fq);
    P(1)=mean(Pq); P(2)=std(Pq);
    R(1)=mean(Rq); R(2)=std(Rq);
    nmi(1)=mean(nmiq); nmi(2)=std(nmiq);
    AR(1)=mean(ARq); AR(2)=std(ARq);
end