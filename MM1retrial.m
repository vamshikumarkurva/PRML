clock =0 ;
N = input('enter the number of Telephone lines: ');
AT = input('enter mean time between two succesive arrivals(Arrival time): ');
ST = input('enter mean length of call duration(service time): ');
RT = input('enter the retrial time: ');
LINKS = input('enter the number of servers(Links): ');
end_time = input('enter the time for simulation: ');
processed = 0;
busy = 0;
blocked = 0;
completed = 0;
links_in_use = 0;
status_lines = zeros(N,1);
status_links = zeros(LINKS,1);
lines_links = zeros(N,1);
CAT = round(exprnd(AT));  
service_time = zeros(LINKS);
retrial = inf(N,1);
waittime = zeros(N,1);
SCT = inf(N,1);

while (clock < end_time)
    not_busy_lines = find(status_lines==0);
    not_busy_links = find(status_links==0);
    [min_time,min_ind] = min(SCT);
    if ~all(retrial==inf)
        MIN = min(retrial);
        clock1 = min(CAT,min_time);
        clock1 = min(clock1,MIN);
    else
        clock1 = min(CAT,min_time);
    end
    
    if (clock1==CAT || clock1==MIN)
       if clock1==CAT
           CAT = CAT + round(exprnd(AT))+1;
           if size(not_busy_lines,1)>2
               call = randsample(not_busy_lines,2);
           else
               clock = clock1;
               continue;
           end
       elseif clock1==MIN
           call = find(retrial==MIN);
       end
           service_time = round(exprnd(ST));
       if (status_lines(call(1))==1 || status_lines(call(2))==1)
           busy = busy+1;
           processed = processed+1;
           retrial(call(1)) = clock1+RT;
           retrial(call(2)) = clock1+RT;
           waittime(call(1)) = waittime(call(1))+RT;
           waittime(call(2)) = waittime(call(2))+RT;
           status_lines(call(1)) = 2;
           status_lines(call(2)) = 2;
       elseif all(status_links==1)
           blocked = blocked+1;
           processed = processed+1;
           retrial(call(1)) = clock1+RT;
           retrial(call(2)) = clock1+RT;
           waittime(call(1)) = waittime(call(1))+RT;
           waittime(call(2)) = waittime(call(2))+RT;           
           status_lines(call(1)) = 2;
           status_lines(call(2)) = 2;
       else
           lines_links(call(1)) = randsample(not_busy_links,1);
           lines_links(call(2)) = lines_links(call(1));
           status_lines(call(1)) = 1;
           status_lines(call(2)) = 1; 
           status_links(lines_links(call(1)))=1;
           SCT(call(1)) = clock1+service_time;
           SCT(call(2)) = clock1+service_time;
           retrial(call(1)) = inf;
           retrial(call(2)) = inf;
       end
    elseif clock1==min_time
        call = find(SCT==min_time);
        status_lines(call) = 0;
        status_links(lines_links(call(1)))=0;
        lines_links(call(1)) = 0;
        lines_links(call(2)) = 0;
        completed = completed + 1;
        processed = processed + 1;
        SCT(call) = inf;
    end 
    clock = clock1;
end

avg_wait_time = sum(waittime)/(2*processed); 
disp('Processed calls: '), disp(processed);
disp('Blocked calls: '), disp(blocked);
%disp('Busy calls: '), disp(busy);
disp('completed calls: '), disp(completed);
disp('average waiting time'), disp(avg_wait_time);
